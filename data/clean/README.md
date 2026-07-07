# Clean Data

Tutorial-ready files derived from the raw survey, lookup, boundary, Census, and population sources.

## Lookups

- `msoa_lad_lookup.csv`: Wales MSOA-to-local-authority lookup used by the clean-data scripts.
- `dvla_lad_lookup.csv`: National Survey for Wales `DvLA` code to current local-authority code bridge.

Notes:

- `DvLA == 7` maps to Powys, `W06000023`.
- `DvLA == 17` maps to Merthyr Tydfil, `W06000024`.

## Survey Files

- `survey_microdata_msoayear.csv`: one row per respondent-wave. Includes harmonised survey fields, synthetic MSOA assignment, walking indicator, weights, and respondent predictors.
- `surveyind_freqwalk.csv`: one row per respondent-wave for the walking estimand. Includes `unit_id`, `time_id`, `area_id`, `strata`, `indicator`, `weight`, and unit predictors `sex`, `age`, and `no_car`.
- `domain_msoayear.csv`: complete MSOA-year estimation grid, including cells with no survey respondents.

## Walking Auxiliary Files

- `areapred_freqwalk.csv`: one row per MSOA with raw and z-scored area predictors for daily walking models, including Census 2021 mean age among usual residents aged 16+ and proportion female from Nomis bulk tables.
- `poststrat_freqwalk.csv`: one row per MSOA-year-age-sex cell with `pop_count` for poststratification.

The 2025 survey wave uses the latest supplied mid-2024 population denominator and records it in `source_population_year`.

## Build Scripts

- `scripts/build_clean_survey_microdata.py`: builds `survey_microdata_msoayear.csv`, `surveyind_freqwalk.csv`, and `domain_msoayear.csv`.
- `scripts/build_walking_area_predictors.py`: builds `areapred_freqwalk.csv`.
- `scripts/build_walking_poststrat.py`: builds `poststrat_freqwalk.csv`.
