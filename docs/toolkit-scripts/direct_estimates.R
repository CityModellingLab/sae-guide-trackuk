library(dplyr)
library(survey)
library(tidyr)

# File import
# Provide survey_data with: area_id, time_id, outcome, weight, optional strata.
# Provide domain_frame with: area_id, time_id for every target area-time cell.

check_cols <- function(data, cols, name) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop(name, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

check_cols(survey_data, c("area_id", "time_id", "outcome", "weight"), "survey_data")
check_cols(domain_frame, c("area_id", "time_id"), "domain_frame")

survey_data <- survey_data |>
  filter(!is.na(area_id), !is.na(time_id), !is.na(outcome), !is.na(weight))

design_args <- list(
  ids = ~1,
  weights = ~weight,
  data = survey_data
)

if ("strata" %in% names(survey_data)) {
  survey_data <- filter(survey_data, !is.na(strata))
  design_args$data <- survey_data
  design_args$strata <- ~strata
}

survey_design <- do.call(svydesign, design_args)

direct <- svyby(
  ~outcome,
  ~area_id + time_id,
  survey_design,
  svymean,
  vartype = "se",
  na.rm = TRUE,
  keep.names = FALSE
) |>
  as_tibble() |>
  transmute(
    area_id,
    time_id,
    estimate = outcome,
    variance = se^2,
    lower = estimate - qnorm(0.975) * se,
    upper = estimate + qnorm(0.975) * se,
    ci_width = upper - lower
  )

counts <- count(survey_data, area_id, time_id, name = "n")

direct_estimates <- domain_frame |>
  select(area_id, time_id) |>
  left_join(direct, by = c("area_id", "time_id")) |>
  left_join(counts, by = c("area_id", "time_id")) |>
  mutate(
    n = replace_na(n, 0L),
    remarks = "Direct estimation"
  )

# File export
# Export direct_estimates using the common output schema:
# area_id, time_id, estimate, variance, lower, upper, ci_width, n, remarks.
