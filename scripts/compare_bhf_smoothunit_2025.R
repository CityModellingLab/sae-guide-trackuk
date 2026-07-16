# Quick check: project unit-level area-summary code vs SUMMER::smoothUnit().

pkgs <- c("dplyr", "readr", "tidyr", "tibble", "sf", "spdep", "survey", "INLA", "SUMMER", "ggplot2", "scales")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
invisible(lapply(pkgs, library, character.only = TRUE))

year <- 2025
n_draw <- 300
var_tol <- 1e-8

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(script_path)) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) normalizePath(sub("^--file=", "", file_arg[1])) else NA_character_
}

root <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
data_dir <- file.path(root, "data")
out_dir <- file.path(root, "outputs", "compare_bhf_smoothunit_2025")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_limit <- function(x) c(0, max(0.1, quantile(x[is.finite(x)], 0.99, na.rm = TRUE, names = FALSE)))
pct <- label_percent(accuracy = 0.1)

# Inputs.
domain <- read_csv(file.path(data_dir, "clean", "domain_msoayear.csv"), show_col_types = FALSE) |>
  filter(time_id == year) |>
  select(area_id, time_id, area_inla_id, msoa_name) |>
  arrange(area_inla_id)

survey_2025 <- read_csv(file.path(data_dir, "clean", "surveyind_freqcyc.csv"), show_col_types = FALSE) |>
  filter(time_id == year) |>
  transmute(
    unit_id,
    area_id,
    time_id,
    strata,
    weight,
    freq_cyclist,
    age,
    female = as.integer(sex == 2),
    caraccess
  ) |>
  drop_na(area_id, strata, weight, freq_cyclist, age, female, caraccess)

area_x <- read_csv(file.path(data_dir, "clean", "areapred_freqcyc.csv"), show_col_types = FALSE) |>
  transmute(
    area_id,
    age = age_adult_mean,
    female = female_adult_share,
    caraccess = caraccess_adult_share
  )

agg <- domain |>
  left_join(area_x, by = "area_id") |>
  drop_na(age, female, caraccess) |>
  arrange(area_inla_id)

boundaries <- st_read(file.path(data_dir, "boundaries-msoa.geojson"), quiet = TRUE) |>
  rename(area_id = all_of("MSOA21CD")) |>
  filter(area_id %in% agg$area_id) |>
  left_join(agg |> select(area_id, area_inla_id), by = "area_id") |>
  arrange(area_inla_id)

nb <- poly2nb(boundaries, queen = TRUE, row.names = boundaries$area_id)
graph_file <- tempfile(fileext = ".adj")
nb2INLA(file = graph_file, nb = nb)

adj_mat <- nb2mat(nb, style = "B", zero.policy = TRUE)
rownames(adj_mat) <- boundaries$area_id
colnames(adj_mat) <- boundaries$area_id

# Direct-estimate context.
design <- svydesign(ids = ~1, strata = ~strata, weights = ~weight, data = survey_2025)

direct_est <- svyby(~freq_cyclist, ~area_id, design, svymean, vartype = "se", na.rm = TRUE, keep.names = FALSE) |>
  as_tibble() |>
  transmute(area_id, direct_estimate = freq_cyclist, variance = se^2)

direct <- domain |>
  left_join(direct_est, by = "area_id") |>
  left_join(count(survey_2025, area_id, name = "n"), by = "area_id") |>
  mutate(n = replace_na(n, 0L))

max_var <- max(direct$variance[direct$variance > 0], na.rm = TRUE)
if (!is.finite(max_var)) stop("No positive direct variances found.", call. = FALSE)

direct <- direct |>
  mutate(
    adj_variance = case_when(
      is.na(variance) ~ max_var,
      variance <= var_tol ~ max_var,
      TRUE ~ variance
    )
  )

# INLA fit.
model_data <- survey_2025 |>
  left_join(domain |> select(area_id, area_inla_id), by = "area_id") |>
  drop_na(area_inla_id) |>
  arrange(area_inla_id)

area_ids <- sort(unique(agg$area_inla_id))

bhf_fit <- inla(
  freq_cyclist ~
    f(area_inla_id, model = "bym2", graph = graph_file, scale.model = TRUE, constr = TRUE) +
    age + female + caraccess,
  family = "binomial",
  Ntrials = 1,
  data = model_data,
  control.predictor = list(compute = TRUE),
  control.compute = list(config = TRUE)
)

set.seed(100)
draws <- inla.posterior.sample(n = n_draw, result = bhf_fit)
fixed_terms <- rownames(bhf_fit$summary.fixed)

fixed_draws <- sapply(draws, \(draw) {
  setNames(as.numeric(draw$latent[paste0(fixed_terms, ":1"), 1]), fixed_terms)
})

area_draws <- sapply(draws, \(draw) {
  setNames(as.numeric(draw$latent[paste0("area_inla_id:", area_ids), 1]), area_ids)
})

X <- model.matrix(~ age + female + caraccess, data = agg)[, fixed_terms, drop = FALSE]
inla_pred <- sapply(seq_along(draws), \(i) {
  plogis(as.numeric(X %*% fixed_draws[, i]) + area_draws[as.character(agg$area_inla_id), i])
})

# SUMMER fit.
summer_fit <- smoothUnit(
  formula = freq_cyclist ~ age + female + caraccess,
  family = "binomial",
  domain = ~area_id,
  design = design,
  X.pop = agg |> select(area_id, age, female, caraccess),
  adj.mat = adj_mat,
  n.sample = n_draw,
  return.samples = TRUE
)

summer_est <- as_tibble(summer_fit$bym2.model.est) |>
  transmute(area_id = domain, summer = mean)

# Outputs.
final_results <- agg |>
  transmute(
    area_id,
    time_id,
    msoa_name,
    inla = rowMeans(inla_pred)
  ) |>
  left_join(
    direct |> select(area_id, direct_estimate, variance, adj_variance, n),
    by = "area_id"
  ) |>
  left_join(summer_est, by = "area_id") |>
  mutate(diff = inla - summer)

write_csv(final_results, file.path(out_dir, "final_results.csv"))

limits <- plot_limit(c(final_results$inla, final_results$summer))

parity_plot <- final_results |>
  ggplot(aes(summer, inla)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed") +
  geom_point(alpha = 0.75, size = 1.8) +
  coord_equal(xlim = limits, ylim = limits) +
  scale_x_continuous(labels = pct) +
  scale_y_continuous(labels = pct) +
  labs(title = "Unit summaries: INLA vs smoothUnit", x = "SUMMER", y = "INLA") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "parity_plot.png"), parity_plot, width = 7, height = 6, dpi = 300)

message("Wrote: ", out_dir)
