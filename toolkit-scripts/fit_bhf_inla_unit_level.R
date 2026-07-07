# Spatio-temporal BHF-style unit-level SAE model with R-INLA
#
# Edit the config block, then run this script from the project root.
# Required inputs:
# - unit file: respondent records with outcome and predictors
# - aggregation file: one row per area-time cell with target-level predictor values
# - area boundaries: polygons used to build the spatial adjacency graph

library(dplyr)
library(readr)
library(sf)
library(spdep)
library(INLA)
library(tibble)

config <- list(
  unit_file = "data/clean/unit_model_data.csv",
  aggregation_file = "data/clean/unit_aggregation_frame.csv",
  boundary_file = "data/boundaries.geojson",
  output_dir = "outputs",

  area_col = "area_id",
  time_col = "time_id",
  outcome_col = "indicator",
  predictor_cols = c("no_car"),

  family = "gaussian"
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

make_formula <- function(predictor_cols) {
  fixed_terms <- if (length(predictor_cols) == 0) {
    "1"
  } else {
    paste(predictor_cols, collapse = " + ")
  }

  as.formula(
    paste(
      "indicator ~",
      fixed_terms,
      "+ f(time_inla_id, model = 'rw1', constr = TRUE)",
      "+ f(area_spatial_id, model = 'bym2', graph = spatial_graph_file,",
      "scale.model = TRUE, constr = TRUE)"
    )
  )
}

unit_data <- read_csv(config$unit_file, show_col_types = FALSE) |>
  rename_column(config$area_col, "area_id") |>
  rename_column(config$time_col, "time_id") |>
  rename_column(config$outcome_col, "indicator")

aggregation_frame <- read_csv(config$aggregation_file, show_col_types = FALSE) |>
  rename_column(config$area_col, "area_id") |>
  rename_column(config$time_col, "time_id")

boundaries <- st_read(config$boundary_file, quiet = TRUE) |>
  rename_column(config$area_col, "area_id")

check_columns(
  unit_data,
  c("area_id", "time_id", "indicator", config$predictor_cols),
  "unit_data"
)
check_columns(
  aggregation_frame,
  c("area_id", "time_id", config$predictor_cols),
  "aggregation_frame"
)
check_columns(boundaries, "area_id", "boundaries")

area_index <- aggregation_frame |>
  distinct(area_id) |>
  arrange(area_id) |>
  mutate(area_inla_id = row_number())

time_index <- aggregation_frame |>
  distinct(time_id) |>
  arrange(time_id) |>
  mutate(time_inla_id = row_number())

boundaries_sorted <- boundaries |>
  inner_join(area_index, by = "area_id") |>
  arrange(area_inla_id)

area_neighbours <- poly2nb(
  boundaries_sorted,
  queen = TRUE,
  row.names = boundaries_sorted$area_id
)

spatial_graph_file <- tempfile(fileext = ".adj")
nb2INLA(file = spatial_graph_file, nb = area_neighbours)

model_rows <- unit_data |>
  filter(
    !is.na(area_id),
    !is.na(time_id),
    !is.na(indicator)
  ) |>
  filter(if_all(all_of(config$predictor_cols), \(x) !is.na(x))) |>
  left_join(area_index, by = "area_id") |>
  left_join(time_index, by = "time_id") |>
  mutate(row_type = "model", area_spatial_id = area_inla_id) |>
  select(
    row_type,
    indicator,
    all_of(config$predictor_cols),
    time_inla_id,
    area_spatial_id
  )

aggregation_rows <- aggregation_frame |>
  left_join(area_index, by = "area_id") |>
  left_join(time_index, by = "time_id") |>
  mutate(row_type = "aggregation", indicator = NA_real_, area_spatial_id = area_inla_id) |>
  select(
    row_type,
    indicator,
    all_of(config$predictor_cols),
    time_inla_id,
    area_spatial_id
  )

inla_input <- bind_rows(model_rows, aggregation_rows)
aggregation_row_ids <- which(inla_input$row_type == "aggregation")

bhf_formula <- make_formula(config$predictor_cols)

fit <- inla(
  formula = bhf_formula,
  family = config$family,
  data = inla_input,
  control.predictor = list(compute = TRUE),
  control.compute = list(dic = TRUE, waic = TRUE),
  verbose = FALSE
)

model_comparison <- tibble(
  model = "BHF temporal-spatial",
  waic = fit$waic$waic,
  effective_parameters = fit$waic$p.eff,
  dic = fit$dic$dic
)

estimates <- aggregation_frame |>
  transmute(
    area_id,
    time_id,
    estimate = fit$summary.fitted.values$mean[aggregation_row_ids],
    lower = fit$summary.fitted.values$`0.025quant`[aggregation_row_ids],
    upper = fit$summary.fitted.values$`0.975quant`[aggregation_row_ids],
    method = "BHF temporal-spatial"
  )

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(model_comparison, file.path(config$output_dir, "bhf_inla_model_comparison.csv"))

write_csv(estimates, file.path(config$output_dir, "bhf_inla_estimates.csv"))
