# Clean data

This folder contains the cleaned inputs used by the analysis. Fitted estimates and posterior outputs are stored in `outputs/`.

## Geography

- `msoa_lad_lookup.csv`: Wales 2021 MSOAs and their current Local Authorities.
- `dvla_lad_lookup.csv`: National Survey for Wales `DvLA` codes matched to current Local Authority codes.
- `msoa_queen.adj`: INLA adjacency graph for Wales MSOAs. Two MSOAs are neighbours when they share a boundary or corner.
- `domain_msoayear.csv`: complete grid of 408 MSOAs and the three analysis periods, coded as 2021, 2023 and 2025.

## Survey data

- `survey_microdata_msoayear.csv`: one row per respondent-wave, with harmonised survey fields, survey weights, respondent predictors and a synthetic MSOA assignment.
- `surveyind_freqcyc.csv`: adult respondent-wave records used for the cycling models. `freq_cyclist` equals 1 for adults who cycle at least several times a week and 0 for other valid responses. Missing outcomes and weights are retained for downstream filtering.

The cycling file contains `unit_id`, `time_id`, `area_id`, `strata`, `freq_cyclist`, `weight`, `sex`, `age` and `caraccess`.

## Census and population data

- `areapred_freqcyc.csv`: Census 2021 MSOA predictors: population density, mean adult age, female share and car-access share.
- `poststrat_freqcyc_age_sex_212325.csv`: age-by-sex population cells for all three periods. The population sources are mid-2020, mid-2022 and mid-2024 respectively.
- `poststrat_freqcyc_age_sex_car_21.csv`: Census 2021 age-by-sex-by-car-access cells. The same 2021 counts are used for all three periods because equivalent joint counts are not available for each period.
- `validation_bicycle_to_work_share.csv`: Census 2021 bicycle-to-work share for each MSOA. This is a related benchmark, not the same outcome as frequent cycling.

## Data limitations

- The survey does not identify respondents' MSOAs. Each respondent is assigned at random to an MSOA within their Local Authority using a fixed seed. The resulting MSOA estimates are analysis examples, not observed local estimates.
- The `DvLA` lookup maps code 7 to Powys (`W06000023`) and code 17 to Merthyr Tydfil (`W06000024`).
- The 2024–25 survey does not contain exact age. Missing ages are generated within `DvAgeGrp7` bands using the observed 2020–21 age distribution and a fixed seed. `age_source` records whether age is observed or generated.
- The travel weight is used in waves where it is supplied. The adult weight is used when a wave does not contain the travel-weight field.
- The 2025 age-by-sex cells use mid-2024 population estimates. The source's youngest band is age 15–19, so it is used as the nearest available match to the survey population aged 16 and over.

## Build sources

- `scripts/build_clean_survey_microdata.py`: survey microdata, cycling analysis file and domain grid.
- `scripts/build_cycling_area_predictors.py`: MSOA predictors.
- `scripts/build_cycling_poststrat_shallow3y.py`: age-by-sex cells.
- `scripts/build_cycling_poststrat_rich1y.py`: age-by-sex-by-car-access cells.
- `scripts/build_validation_bicycle_to_work.py`: bicycle-to-work benchmark.
- `dataprep.qmd`: Queen-contiguity adjacency graph.
