# Quick check: project unit-level area-summary code vs SUMMER::smoothUnit().

pkgs <- c("dplyr", "readr", "tidyr", "tibble", "sf", "spdep", "survey", "INLA", "SUMMER", "ggplot2", "scales")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
invisible(lapply(pkgs, library, character.only = TRUE))

year <- 2025
n_draw <- 300

root <- here::here("tutorial-site")
data_dir <- file.path(root, "data")
out_dir <- file.path(root, "outputs", "compare_smoothUnit_2025")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_limit <- function(x) c(0, max(0.1, quantile(x[is.finite(x)], 0.99, na.rm = TRUE, names = FALSE)))
pct <- label_percent(accuracy = 0.1)

# Inputs.
spatial_graph_file <- file.path(data_dir, "clean", "msoa_queen.adj")

domain <- read_csv(file.path(data_dir, "clean", "domain_msoayear.csv"), show_col_types = FALSE) |>
  filter(time_id == year) |>
  select(area_id, time_id, area_inla_id, msoa_name) |>
  arrange(area_inla_id)

survey_2025 <- read_csv(file.path(data_dir, "clean", "surveyind_freqcyc.csv"), show_col_types = FALSE) |>
  filter(time_id == year) |>
  transmute(
    area_id,
    strata,
    weight,
    freq_cyclist,
    age,
    female = as.integer(sex == 2),
    caraccess
  ) |>
  drop_na(area_id, strata, weight, freq_cyclist, age, female, caraccess)

area_x <- read_csv(file.path(data_dir, "clean", "areapred_freqcyc.csv"), show_col_types = FALSE) |>
  select(
    area_id,
    age = age_adult_mean,
    female = female_adult_share,
    caraccess = caraccess_adult_share
  )

aggregation_data <- domain |>
  left_join(area_x, by = "area_id") |>
  drop_na(age, female, caraccess) |>
  select(area_id, msoa_name, time_id, age, female, caraccess, area_inla_id) |>
  arrange(area_inla_id)

boundaries <- st_read(file.path(data_dir, "boundaries-msoa.geojson"), quiet = TRUE) |>
  rename(area_id = MSOA21CD) |>
  filter(area_id %in% aggregation_data$area_id) |>
  left_join(aggregation_data |> select(area_id, area_inla_id), by = "area_id") |>
  arrange(area_inla_id)

nb <- poly2nb(boundaries, queen = TRUE, row.names = boundaries$area_id)

adj_mat <- nb2mat(nb, style = "B", zero.policy = TRUE)
rownames(adj_mat) <- boundaries$area_id
colnames(adj_mat) <- boundaries$area_id

# Survey design for SUMMER.
design <- svydesign(ids = ~1, strata = ~strata, weights = ~weight, data = survey_2025)

# INLA fit.
model_data <- survey_2025 |>
  select(area_id, freq_cyclist, age, female, caraccess) |>
  left_join(domain |> select(area_id, area_inla_id), by = "area_id") |>
  arrange(area_inla_id)

bhf_fit <- inla(
  freq_cyclist ~
    1 +
    f(area_inla_id, model = "bym2", graph = spatial_graph_file, scale.model = TRUE) +
    age + female + caraccess,
  family = "binomial",
  data = model_data,
  control.compute = list(config = TRUE)
)

set.seed(100)
posterior_draws <- inla.posterior.sample(
  n = n_draw,
  result = bhf_fit
)

extract_draws <- function(posterior_draws, nodes) {
  matrix(
    sapply(posterior_draws, \(draw) draw$latent[nodes, 1]),
    nrow = length(nodes)
  )
}

predict_from_draws <- function(posterior_draws, prediction_data, fixed_term_formula) {
  prediction_matrix <- model.matrix(
    fixed_term_formula,
    prediction_data
  )

  fixed_terms <- colnames(prediction_matrix)
  area_ids <- sort(unique(prediction_data$area_inla_id))

  fixed_draws <- extract_draws(posterior_draws, paste0(fixed_terms, ":1"))
  area_draws <- extract_draws(posterior_draws, paste0("area_inla_id:", area_ids))

  linear_predictor <- prediction_matrix %*% fixed_draws +
    area_draws[match(prediction_data$area_inla_id, area_ids), ]

  plogis(linear_predictor)
}

prediction_draws <- predict_from_draws(
  posterior_draws,
  aggregation_data,
  ~ age + female + caraccess
)

# SUMMER fit.
summer_fit <- smoothUnit(
  formula = freq_cyclist ~ age + female + caraccess,
  family = "binomial",
  domain = ~area_id,
  design = design,
  X.pop = aggregation_data |> select(area_id, age, female, caraccess),
  adj.mat = adj_mat,
  n.sample = n_draw,
  return.samples = TRUE
)

summer_est <- as_tibble(summer_fit$bym2.model.est) |>
  transmute(area_id = domain, summer = mean)

# Outputs.
final_results <- aggregation_data |>
  select(area_id, time_id, msoa_name) |>
  bind_cols(tibble(inla = rowMeans(prediction_draws))) |>
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