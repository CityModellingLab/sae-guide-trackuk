# Bayesian MRP with population cells
#
# Input objects:
# - survey_unit_data: unit_id, domain_id, indicator, weight, optional stratum,
#   plus unit-level covariate_* columns used in the model
# - population_cell_data: cell_id, domain_id, pop_count, and the same
#   unit-level covariate_* columns used in the model

library(dplyr)
library(INLA)

predictor_vars <- c("covariate_1", "covariate_2", "covariate_3")
domain_levels <- sort(unique(population_cell_data$domain_id))

model_levels <- lapply(
  predictor_vars,
  function(var) sort(unique(population_cell_data[[var]]))
)
names(model_levels) <- predictor_vars

add_model_ids <- function(data) {
  for (var in predictor_vars) {
    data[[paste0(var, "_id")]] <- match(data[[var]], model_levels[[var]])
  }

  data$domain_re_id <- match(data$domain_id, domain_levels)
  data
}

survey_model <- survey_unit_data |>
  transmute(
    row_type = "survey",
    indicator,
    domain_id,
    pop_count = NA_real_,
    across(all_of(predictor_vars))
  ) |>
  add_model_ids()

poststrat_model <- population_cell_data |>
  transmute(
    row_type = "poststrat",
    indicator = NA_real_,
    domain_id,
    pop_count,
    across(all_of(predictor_vars))
  ) |>
  add_model_ids() |>
  filter(pop_count > 0)

model_data <- bind_rows(survey_model, poststrat_model) |>
  mutate(model_row = row_number())

mrp_formula <- indicator ~ 1 +
  f(covariate_1_id, model = "iid") +
  f(covariate_2_id, model = "iid") +
  f(covariate_3_id, model = "iid") +
  f(domain_re_id, model = "iid")

mrp_fit <- INLA::inla(
  formula = mrp_formula,
  family = "gaussian",  # change family to match the outcome
  data = model_data,
  control.predictor = list(compute = TRUE),
  control.compute = list(config = TRUE),
  verbose = FALSE
)

poststrat_rows <- model_data |>
  filter(row_type == "poststrat") |>
  select(model_row, domain_id, pop_count)

posterior_samples <- INLA::inla.posterior.sample(
  n = 1000,
  result = mrp_fit
)

predictor_names <- paste0("Predictor:", poststrat_rows$model_row)

cell_prediction_draws <- vapply(
  posterior_samples,
  function(sample) {
    as.numeric(sample$latent[predictor_names, 1])
  },
  numeric(nrow(poststrat_rows))
)

domain_index <- split(seq_len(nrow(poststrat_rows)), poststrat_rows$domain_id)

# Standardise results for diagnostics
unit_level_results <- lapply(names(domain_index), function(area) {
  row_ids <- domain_index[[area]]
  weights <- poststrat_rows$pop_count[row_ids] / sum(poststrat_rows$pop_count[row_ids])
  area_draws <- colSums(cell_prediction_draws[row_ids, , drop = FALSE] * weights)

  tibble(
    domain_id = area,
    method = "MRP",
    estimate = mean(area_draws),
    se = sd(area_draws),
    lower = as.numeric(quantile(area_draws, 0.025)),
    upper = as.numeric(quantile(area_draws, 0.975)),
    uncertainty_type = "credible_interval"
  )
}) |>
  bind_rows()
