# Quick check: project FH code vs SUMMER::smoothArea().

pkgs <- c("dplyr", "readr", "tidyr", "tibble", "sf", "spdep", "survey", "INLA", "SUMMER", "ggplot2", "scales")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
invisible(lapply(pkgs, library, character.only = TRUE))

year <- 2025
n_draw <- 300
var_tol <- 1e-8
pseudo_high_var <- 1e8

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(script_path)) {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script_path <- if (length(file_arg) > 0) normalizePath(sub("^--file=", "", file_arg[1])) else NA_character_
}

root <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
data_dir <- file.path(root, "data")
out_dir <- file.path(root, "outputs", "compare_fh_smootharea_2025")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_limit <- function(x) c(0, max(0.1, quantile(x[is.finite(x)], 0.99, na.rm = TRUE, names = FALSE)))
pct <- label_percent(accuracy = 0.1)

# Inputs.
domain <- read_csv(file.path(data_dir, "clean", "domain_msoayear.csv"), show_col_types = FALSE) |>
  filter(time_id == year) |>
  select(area_id, time_id, area_inla_id, msoa_name)

survey_2025 <- read_csv(file.path(data_dir, "clean", "surveyind_freqcyc.csv"), show_col_types = FALSE) |>
  filter(time_id == year) |>
  drop_na(freq_cyclist, weight, area_id, strata)

area_x <- read_csv(file.path(data_dir, "clean", "areapred_freqcyc.csv"), show_col_types = FALSE) |>
  select(area_id, age_adult_mean, caraccess_adult_share, pop_density) |>
  mutate(
    across(
      c(age_adult_mean, caraccess_adult_share, pop_density),
      ~ as.numeric(scale(.x)),
      .names = "z_{.col}"
    )
  )

boundaries <- st_read(file.path(data_dir, "boundaries-msoa.geojson"), quiet = TRUE) |>
  rename(area_id = MSOA21CD) |>
  filter(area_id %in% domain$area_id) |>
  left_join(domain |> select(area_id, area_inla_id), by = "area_id") |>
  arrange(area_inla_id)

nb <- poly2nb(boundaries, queen = TRUE, row.names = boundaries$area_id)
graph_file <- tempfile(fileext = ".adj")
nb2INLA(file = graph_file, nb = nb)

adj_mat <- nb2mat(nb, style = "B", zero.policy = TRUE)
rownames(adj_mat) <- boundaries$area_id
colnames(adj_mat) <- boundaries$area_id

# Direct estimates.
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
    has_direct = !is.na(direct_estimate),
    adj_variance = case_when(
      is.na(variance) ~ max_var,
      variance <= var_tol ~ max_var,
      TRUE ~ variance
    )
  )

# INLA fit.
model_data <- direct |>
  left_join(area_x, by = "area_id") |>
  arrange(area_inla_id)

fh_fit <- inla(
  direct_estimate ~
    f(area_inla_id, model = "bym2", graph = graph_file, scale.model = TRUE, constr = TRUE) +
    z_age_adult_mean + z_pop_density + z_caraccess_adult_share,
  family = "gaussian",
  data = model_data,
  scale = 1 / model_data$adj_variance,
  control.family = list(hyper = list(prec = list(initial = 0, fixed = TRUE))),
  control.predictor = list(compute = TRUE),
  control.compute = list(config = TRUE)
)

set.seed(100)
inla_draws <- inla.posterior.sample(n = n_draw, result = fh_fit)
inla_pred <- sapply(inla_draws, \(draw) {
  as.numeric(draw$latent[startsWith(rownames(draw$latent), "Predictor"), 1])
})

# SUMMER fit.
smooth_direct <- direct |>
  transmute(
    area_id,
    direct_estimate = if_else(has_direct, direct_estimate, weighted.mean(survey_2025$freq_cyclist, survey_2025$weight)),
    variance = if_else(has_direct, adj_variance, pseudo_high_var)
  )

summer_fit <- smoothArea(
  formula = freq_cyclist ~ z_age_adult_mean + z_pop_density + z_caraccess_adult_share,
  domain = ~area_id,
  direct.est = smooth_direct,
  adj.mat = adj_mat,
  X.domain = area_x,
  transform = "identity",
  n.sample = n_draw,
  var.tol = var_tol
)

summer_est <- as_tibble(summer_fit$bym2.model.est) |>
  transmute(area_id = domain, summer = mean)

# Outputs.
final_results <- model_data |>
  transmute(
    area_id,
    time_id,
    msoa_name,
    n,
    direct_estimate,
    variance,
    adj_variance,
    inla = rowMeans(inla_pred)
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
  labs(title = "FH: INLA vs smoothArea", x = "SUMMER", y = "INLA") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "parity_plot.png"), parity_plot, width = 7, height = 6, dpi = 300)

message("Wrote: ", out_dir)
