library(dplyr)
library(INLA)

# File import
# Provide survey_data with: area_id, time_id, outcome, and unitpred_* columns.
# Provide domain_frame with: area_id, time_id.
# Provide area_summaries with: area_id and matching *_mean / *_share columns.
# Provide spatial_graph_file as the INLA .adj file ordered by area_inla_id.

family <- "binomial"
unitpred_cols <- c("unitpred_1", "unitpred_2")
area_summary_map <- c(
  unitpred_1 = "unitpred_1_mean",
  unitpred_2 = "unitpred_2_share"
)

check_cols <- function(data, cols, name) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(name, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

rename_by_map <- function(data, map) {
  for (target in names(map)) {
    names(data)[names(data) == unname(map[target])] <- target
  }
  data
}

area_summaries <- area_summaries |>
  rename_by_map(area_summary_map)

check_cols(survey_data, c("area_id", "time_id", "outcome", unitpred_cols), "survey_data")
check_cols(domain_frame, c("area_id", "time_id"), "domain_frame")
check_cols(area_summaries, c("area_id", unitpred_cols), "area_summaries")

area_index <- domain_frame |>
  distinct(area_id) |>
  arrange(area_id) |>
  mutate(area_inla_id = row_number())

time_index <- domain_frame |>
  distinct(time_id) |>
  arrange(time_id) |>
  mutate(time_inla_id = row_number())

model_rows <- survey_data |>
  filter(!is.na(area_id), !is.na(time_id), !is.na(outcome)) |>
  filter(if_all(all_of(unitpred_cols), ~ !is.na(.x))) |>
  left_join(area_index, by = "area_id") |>
  left_join(time_index, by = "time_id") |>
  mutate(row_type = "model")

prediction_rows <- domain_frame |>
  left_join(area_summaries, by = "area_id") |>
  left_join(area_index, by = "area_id") |>
  left_join(time_index, by = "time_id") |>
  mutate(row_type = "prediction", outcome = NA_real_)

inla_data <- bind_rows(
  select(model_rows, row_type, outcome, all_of(unitpred_cols), area_inla_id, time_inla_id),
  select(prediction_rows, row_type, outcome, all_of(unitpred_cols), area_inla_id, time_inla_id)
)

prediction_ids <- which(inla_data$row_type == "prediction")

formula <- as.formula(paste(
  "outcome ~",
  paste(
    c(
      unitpred_cols,
      "f(time_inla_id, model = 'rw1', constr = TRUE)",
      "f(area_inla_id, model = 'bym2', graph = spatial_graph_file, scale.model = TRUE, constr = TRUE)"
    ),
    collapse = " + "
  )
))

fit_args <- list(
  formula = formula,
  family = family,
  data = inla_data,
  control.predictor = list(compute = TRUE),
  control.compute = list(waic = TRUE)
)

if (family == "binomial") {
  fit_args$Ntrials <- 1
}

fit <- do.call(inla, fit_args)

unit_sae_estimates <- prediction_rows |>
  transmute(
    area_id,
    time_id,
    estimate = fit$summary.fitted.values$mean[prediction_ids],
    variance = fit$summary.fitted.values$sd[prediction_ids]^2,
    lower = fit$summary.fitted.values$`0.025quant`[prediction_ids],
    upper = fit$summary.fitted.values$`0.975quant`[prediction_ids],
    ci_width = upper - lower,
    n = NA_integer_,
    remarks = "Unit-level SAE with area summaries"
  )

# File export
# Export unit_sae_estimates using the common output schema:
# area_id, time_id, estimate, variance, lower, upper, ci_width, n, remarks.
