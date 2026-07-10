# Spatio-temporal Fay-Herriot SAE model with R-INLA
#
# Edit the config block, then run this script from the project root.
# Required inputs:
# - domain frame: one row per area-time cell to estimate
# - direct estimates: direct estimate and sampling variance for observed cells
# - area predictors: one row per area
# - area boundaries: polygons used to build the spatial adjacency graph

library(dplyr)
library(readr)
library(sf)
library(spdep)
library(INLA)
library(purrr)
library(tibble)

config <- list(
  domain_file = "data/clean/domain_frame.csv",
  direct_file = "data/clean/direct_estimates.csv",
  predictor_file = "data/clean/area_predictors.csv",
  boundary_file = "data/boundaries.geojson",
  output_dir = "outputs",

  area_col = "area_id",
  time_col = "time_id",
  direct_estimate_col = "direct_estimate",
  variance_col = "variance",

  core_covariates = c(
    "scaled_population_density",
    "scaled_share_no_car_households"
  ),
  extra_covariates = c(
    "scaled_share_work_from_home",
    "scaled_share_students"
  )
)

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

make_formula <- function(
  covariates,
  include_temporal = TRUE,
  include_spatial = FALSE
) {
  fixed_terms <- if (length(covariates) == 0) {
    "1"
  } else {
    paste(covariates, collapse = " + ")
  }

  random_terms <- "f(area_iid_id, model = 'iid')"

  if (include_temporal) {
    random_terms <- c(
      random_terms,
      "f(time_inla_id, model = 'rw1', constr = TRUE)"
    )
  }

  if (include_spatial) {
    random_terms <- paste0(
      "f(area_spatial_id, model = 'bym2', graph = spatial_graph_file, ",
      "scale.model = TRUE, constr = TRUE)"
    )

    if (include_temporal) {
      random_terms <- c(
        random_terms,
        "f(time_inla_id, model = 'rw1', constr = TRUE)"
      )
    }
  }

  as.formula(
    paste(
      "direct_estimate ~",
      paste(c(fixed_terms, random_terms), collapse = " + ")
    )
  )
}

domain_frame <- read_csv(config$domain_file, show_col_types = FALSE) |>
  rename(
    area_id = all_of(config$area_col),
    time_id = all_of(config$time_col)
  )

direct_estimates <- read_csv(config$direct_file, show_col_types = FALSE) |>
  rename(
    area_id = all_of(config$area_col),
    time_id = all_of(config$time_col),
    direct_estimate = all_of(config$direct_estimate_col),
    variance = all_of(config$variance_col)
  )

area_predictors <- read_csv(config$predictor_file, show_col_types = FALSE) |>
  rename(area_id = all_of(config$area_col))

boundaries <- st_read(config$boundary_file, quiet = TRUE) |>
  rename(area_id = all_of(config$area_col))

all_covariates <- c(config$core_covariates, config$extra_covariates)

check_columns(domain_frame, c("area_id", "time_id"), "domain_frame")
check_columns(
  direct_estimates,
  c("area_id", "time_id", "direct_estimate", "variance"),
  "direct_estimates"
)
check_columns(area_predictors, c("area_id", all_covariates), "area_predictors")
check_columns(boundaries, "area_id", "boundaries")

area_index <- domain_frame |>
  distinct(area_id) |>
  arrange(area_id) |>
  mutate(area_inla_id = row_number())

time_index <- domain_frame |>
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

max_variance <- max(direct_estimates$variance[direct_estimates$variance > 0], na.rm = TRUE)

model_data <- domain_frame |>
  left_join(direct_estimates, by = c("area_id", "time_id")) |>
  left_join(area_predictors, by = "area_id") |>
  left_join(area_index, by = "area_id") |>
  left_join(time_index, by = "time_id") |>
  mutate(
    area_iid_id = area_inla_id,
    area_spatial_id = area_inla_id,
    sampling_variance = case_when(
      is.na(direct_estimate) ~ max_variance,
      is.na(variance) ~ max_variance,
      variance <= 1e-8 ~ max_variance,
      TRUE ~ variance
    )
  ) |>
  arrange(time_id, area_inla_id)

model_formulas <- list(
  "M0 Area only" = make_formula(character(0), include_temporal = FALSE),
  "M1 Temporal" = make_formula(character(0)),
  "M2 Core predictors" = make_formula(config$core_covariates),
  "M3 All predictors" = make_formula(all_covariates),
  "M4 Spatial" = make_formula(all_covariates, include_spatial = TRUE)
)

fit_fh_model <- function(formula) {
  inla(
    formula = formula,
    family = "gaussian",
    data = model_data,
    scale = 1 / model_data$sampling_variance,
    control.family = list(
      hyper = list(
        prec = list(initial = 0, fixed = TRUE)
      )
    ),
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, waic = TRUE),
    verbose = FALSE
  )
}

fits <- lapply(model_formulas, fit_fh_model)

model_comparison <- imap_dfr(
  fits,
  \(fit, model_name) {
    tibble(
      model = model_name,
      waic = fit$waic$waic,
      effective_parameters = fit$waic$p.eff
    )
  }
)

best_model_name <- model_comparison |>
  filter(waic == min(waic, na.rm = TRUE)) |>
  slice(1) |>
  pull(model)

best_fit <- fits[[best_model_name]]

estimates <- model_data |>
  transmute(
    area_id,
    time_id,
    direct_estimate,
    variance,
    estimate = best_fit$summary.fitted.values$mean,
    lower = best_fit$summary.fitted.values$`0.025quant`,
    upper = best_fit$summary.fitted.values$`0.975quant`,
    method = best_model_name
  )

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(model_comparison, file.path(config$output_dir, "fh_inla_model_comparison.csv"))

write_csv(estimates, file.path(config$output_dir, "fh_inla_estimates.csv"))
