library(dplyr)
library(INLA)

# File import
# Provide domain_frame with: area_id, time_id.
# Provide direct_estimates with: area_id, time_id, estimate, variance, optional n.
# Provide area_predictors with: area_id and raw areapred_* columns.
# Provide spatial_graph_file as the INLA .adj file ordered by area_inla_id.

areapred_cols <- c("areapred_1", "areapred_2")

check_cols <- function(data, cols, name) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(name, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

z_score <- function(x) {
  x <- as.numeric(x)
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

check_cols(domain_frame, c("area_id", "time_id"), "domain_frame")
check_cols(direct_estimates, c("area_id", "time_id", "estimate", "variance"), "direct_estimates")
check_cols(area_predictors, c("area_id", areapred_cols), "area_predictors")

area_index <- domain_frame |>
  distinct(area_id) |>
  arrange(area_id) |>
  mutate(area_inla_id = row_number())

time_index <- domain_frame |>
  distinct(time_id) |>
  arrange(time_id) |>
  mutate(time_inla_id = row_number())

z_cols <- paste0("z_", areapred_cols)

area_predictors <- area_predictors |>
  mutate(across(all_of(areapred_cols), z_score, .names = "z_{.col}"))

positive_variance <- direct_estimates$variance[direct_estimates$variance > 0]
max_variance <- max(positive_variance, na.rm = TRUE)

model_data <- domain_frame |>
  left_join(direct_estimates, by = c("area_id", "time_id")) |>
  left_join(area_predictors, by = "area_id") |>
  left_join(area_index, by = "area_id") |>
  left_join(time_index, by = "time_id") |>
  mutate(
    sampling_variance = case_when(
      is.na(estimate) ~ max_variance,
      is.na(variance) ~ max_variance,
      variance <= 1e-8 ~ max_variance,
      TRUE ~ variance
    )
  )

formula <- as.formula(paste(
  "estimate ~",
  paste(
    c(
      z_cols,
      "f(time_inla_id, model = 'rw1', scale.model = TRUE, constr = TRUE)",
      "f(area_inla_id, model = 'bym2', graph = spatial_graph_file, scale.model = TRUE, constr = TRUE)"
    ),
    collapse = " + "
  )
))

fit <- inla(
  formula,
  family = "gaussian",
  data = model_data,
  scale = 1 / model_data$sampling_variance,
  control.family = list(hyper = list(prec = list(initial = 0, fixed = TRUE))),
  control.predictor = list(compute = TRUE),
  control.compute = list(waic = TRUE)
)

fh_estimates <- model_data |>
  transmute(
    area_id,
    time_id,
    estimate = fit$summary.fitted.values$mean,
    variance = fit$summary.fitted.values$sd^2,
    lower = fit$summary.fitted.values$`0.025quant`,
    upper = fit$summary.fitted.values$`0.975quant`,
    ci_width = upper - lower,
    n = NA_integer_,
    remarks = "Area-level SAE"
  )

# File export
# Export fh_estimates using the common output schema:
# area_id, time_id, estimate, variance, lower, upper, ci_width, n, remarks.
