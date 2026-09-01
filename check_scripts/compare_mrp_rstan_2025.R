# Compare 2025 MRP implementations: INLA vs RSTAN 
# Codes adaptd from https://bookdown.org/jl5522/MRP-case-studies/introduction-to-mister-p.html

pkgs <- c("dplyr", "readr", "tidyr", "tibble", "INLA", "rstanarm", "ggplot2", "scales")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
invisible(lapply(pkgs, library, character.only = TRUE))

year <- 2025
n_draw <- 300
set.seed(100)

root <- here::here("tutorial-site")
data_dir <- file.path(root, "data")
out_dir <- file.path(root, "outputs", "compare_mrp_rstan_2025")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

age_band_levels <- c("16_24", "25_34", "35_49", "50_64", "65_plus")

domain <- read_csv(file.path(data_dir, "clean", "domain_msoayear.csv"), show_col_types = FALSE) |>
  filter(time_id == year) |>
  select(area_id, time_id, area_inla_id, msoa_name) |>
  arrange(area_inla_id)

model_data <- read_csv(file.path(data_dir, "clean", "surveyind_freqcyc.csv"), show_col_types = FALSE) |>
  filter(time_id == year) |>
  mutate(
    age_band = cut(
      age,
      breaks = c(16, 25, 35, 50, 65, Inf),
      labels = age_band_levels,
      right = FALSE
    ),
    sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female"))
  ) |>
  select(freq_cyclist, area_id, age_band, sex, caraccess) |>
  drop_na() |>
  left_join(domain |> select(area_id, area_inla_id), by = "area_id")

poststrat_data <- read_csv(
  file.path(data_dir, "clean", "poststrat_freqcyc_age_sex_car_21.csv"),
  show_col_types = FALSE
) |>
  filter(time_id == year, pop_count > 0) |>
  mutate(
    age_band = factor(age_band, levels = age_band_levels),
    sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female"))
  ) |>
  select(area_id, age_band, sex, caraccess, pop_count) |>
  left_join(domain |> select(area_id, area_inla_id), by = "area_id") |>
  arrange(area_inla_id, age_band, sex, caraccess)

inla_fit <- inla(
  freq_cyclist ~ age_band + sex + caraccess + f(area_inla_id, model = "iid"),
  family = "binomial",
  data = model_data,
  control.compute = list(config = TRUE)
)

stan_fit <- stan_glmer(
  freq_cyclist ~ age_band + sex + caraccess + (1 | area_id),
  family = binomial(link = "logit"),
  data = model_data,
  prior = normal(0, 1, autoscale = TRUE),
  prior_covariance = decov(scale = 0.5),
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  adapt_delta = 0.99,
  refresh = 100,
  seed = 100
)

posterior_draws <- inla.posterior.sample(n_draw, inla_fit)

extract_draws <- function(nodes) {
  sapply(posterior_draws, \(draw) draw$latent[nodes, 1])
}

prediction_matrix <- model.matrix(~ age_band + sex + caraccess, poststrat_data)
area_ids <- sort(unique(poststrat_data$area_inla_id))

inla_cell_draws <- plogis(
  prediction_matrix %*% extract_draws(paste0(colnames(prediction_matrix), ":1")) +
    extract_draws(paste0("area_inla_id:", area_ids))[
      match(poststrat_data$area_inla_id, area_ids),
    ]
)

stan_cell_draws <- t(
  posterior_epred(
    stan_fit,
    newdata = poststrat_data,
    draws = n_draw
  )
)

poststratify <- function(cell_draws) {
  rowsum(cell_draws * poststrat_data$pop_count, poststrat_data$area_id) /
    as.numeric(rowsum(poststrat_data$pop_count, poststrat_data$area_id))
}

inla_prediction_draws <- poststratify(inla_cell_draws)
stan_prediction_draws <- poststratify(stan_cell_draws)

inla_est <- tibble(area_id = rownames(inla_prediction_draws), inla = rowMeans(inla_prediction_draws))
stan_est <- tibble(area_id = rownames(stan_prediction_draws), rstan = rowMeans(stan_prediction_draws))

final_results <- domain |>
  select(area_id, time_id, msoa_name) |>
  left_join(inla_est, by = "area_id") |>
  left_join(stan_est, by = "area_id") |>
  mutate(diff = inla - rstan)

limits <- range(c(final_results$inla, final_results$rstan), na.rm = TRUE)

parity_plot <- ggplot(final_results, aes(rstan, inla)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "dashed") +
  geom_point(alpha = 0.65, colour = "#1f4e5f") +
  coord_equal(xlim = limits, ylim = limits) +
  scale_x_continuous(labels = label_percent(accuracy = 0.1)) +
  scale_y_continuous(labels = label_percent(accuracy = 0.1)) +
  labs(
    title = "MRP: INLA vs RSTAN",
    x = "RSTAN",
    y = "INLA"
  ) +
  theme_minimal(base_size = 12)

write_csv(final_results, file.path(out_dir, "final_results.csv"))
ggsave(file.path(out_dir, "parity_plot.png"), parity_plot, width = 7, height = 7, dpi = 300)

print(final_results)
print(parity_plot)
