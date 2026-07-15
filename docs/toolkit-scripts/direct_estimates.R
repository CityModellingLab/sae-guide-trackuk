# Direct area-time survey estimates
#
# Edit the config block, then run this script from the project root.
# Required inputs:
# - survey file: one row per respondent or respondent-period
# - domain frame: one row per area-time cell to estimate

library(dplyr)
library(readr)
library(survey)
library(tidyr)
library(tibble)

config <- list(
  survey_file = "data/clean/survey_indicator.csv",
  domain_file = "data/clean/domain_frame.csv",
  output_dir = "outputs",
  output_file = "direct_estimates.csv",

  area_col = "area_id",
  time_col = "time_id",
  outcome_col = "indicator",
  weight_col = "weight",
  strata_col = "strata"
)

rename_column <- function(data, source, target) {
  if (is.null(source) || source == target) {
    return(data)
  }

  if (!source %in% names(data)) {
    stop("Column not found: ", source, call. = FALSE)
  }

  names(data)[names(data) == source] <- target
  data
}

check_columns <- function(data, required, data_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(
      data_name,
      " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

survey_data <- read_csv(config$survey_file, show_col_types = FALSE) |>
  rename_column(config$area_col, "area_id") |>
  rename_column(config$time_col, "time_id") |>
  rename_column(config$outcome_col, "indicator") |>
  rename_column(config$weight_col, "weight") |>
  rename_column(config$strata_col, "strata")

domain_frame <- read_csv(config$domain_file, show_col_types = FALSE) |>
  rename_column(config$area_col, "area_id") |>
  rename_column(config$time_col, "time_id")

required_survey_columns <- c("area_id", "time_id", "indicator", "weight")
if (!is.null(config$strata_col)) {
  required_survey_columns <- c(required_survey_columns, "strata")
}

check_columns(survey_data, required_survey_columns, "survey_data")
check_columns(domain_frame, c("area_id", "time_id"), "domain_frame")

survey_data <- survey_data |>
  filter(
    !is.na(area_id),
    !is.na(time_id),
    !is.na(indicator),
    !is.na(weight)
  )

if ("strata" %in% names(survey_data)) {
  survey_data <- survey_data |>
    filter(!is.na(strata))
}

design_args <- list(
  ids = as.formula("~1"),
  weights = as.formula("~weight"),
  data = survey_data
)

if ("strata" %in% names(survey_data)) {
  design_args$strata <- as.formula("~strata")
}

survey_design <- do.call(svydesign, design_args)

direct_estimates <- svyby(
  formula = as.formula("~indicator"),
  by = as.formula("~area_id + time_id"),
  design = survey_design,
  FUN = svymean,
  vartype = "se",
  na.rm = TRUE,
  keep.names = FALSE
) |>
  as_tibble() |>
  transmute(
    area_id,
    time_id,
    direct_estimate = indicator,
    standard_error = se,
    variance = standard_error^2
  )

sample_counts <- survey_data |>
  count(area_id, time_id, name = "n")

results <- domain_frame |>
  left_join(direct_estimates, by = c("area_id", "time_id")) |>
  left_join(sample_counts, by = c("area_id", "time_id")) |>
  mutate(n = replace_na(n, 0L), method = "direct")

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(
  results |>
    select(
      area_id,
      time_id,
      direct_estimate,
      standard_error,
      variance,
      n,
      method
    ),
  file.path(config$output_dir, config$output_file)
)
