# Clean Data

Tutorial-ready files derived from the raw survey, lookup, boundary, Census, and population sources.

## Lookups

- `msoa_lad_lookup.csv`: Wales MSOA-to-local-authority lookup used by the clean-data scripts.
- `dvla_lad_lookup.csv`: National Survey for Wales `DvLA` code to current local-authority code bridge.

## Survey Files

- `survey_microdata_msoayear.csv`: one row per respondent-wave. Includes harmonised survey fields, synthetic MSOA assignment, `strata`, cycling indicator, weights, and respondent predictors.
- `surveyind_freqcyc.csv`: one row per respondent-wave for the cycling estimand. Includes `unit_id`, `time_id`, `area_id`, `strata`, `freq_cyclist`, `weight`, and unit predictors `sex`, `age`, and `caraccess`.
- `domain_msoayear.csv`: complete MSOA-year estimation grid, including cells with no survey respondents.

## Cycling Auxiliary and Validation Files

- `areapred_freqcyc.csv`: one row per MSOA, keyed by `area_id`, with raw area predictors for frequent cycling models. Scaled model terms use the `z_*` prefix when needed.
- `poststrat_freqcyc_age_sex_212325.csv`: one row per MSOA-year-age-sex cell with `pop_count` for shallow poststratification. Age bands are grouped to match the rich poststratification table.
- `poststrat_freqcyc_age_sex_car_21.csv`: one row per MSOA-year-age-sex-car-access cell with `pop_count` for richer poststratification.
- `validation_bicycle_to_work_share.csv`: one row per Wales MSOA with a simple 2021 bicycle-to-work benchmark. The numerator is the TS061 bicycle-to-work count and the denominator is the TS061 total usual residents aged 16+ in employment.

## Notes

- Updated MSOA codes with mismatch from survey microdata.
    - `DvLA == 7` maps to Powys, `W06000023`.
    - `DvLA == 17` maps to Merthyr Tydfil, `W06000024`.
- The 2024-25 raw survey file does not contain exact respondent `Age`. The clean-data builder derives a synthetic numeric `age` from `DvAgeGrp7` by sampling from the 2021 observed age mix within each age-group interval, with deterministic random jitter and a cap of 95 for the open-ended 75+ group. The `survey_microdata_msoayear.csv` file records this in `age_source`.
- The shallow 2025 poststratification table uses the latest supplied mid-2024 population denominator, as the mid-2025 estimates have not been published as of July 2026.
- In the shallow poststratification table, the youngest available ONS age band is `15 to 19`, so age 15+ is used as the closest available match to the NSW adult definition of age 16+.

## Build Scripts

- `scripts/build_clean_survey_microdata.py`: builds `survey_microdata_msoayear.csv`, `surveyind_freqcyc.csv`, and `domain_msoayear.csv`.
- `scripts/build_cycling_area_predictors.py`: builds `areapred_freqcyc.csv`.
- `scripts/build_cycling_poststrat_shallow3y.py`: builds `poststrat_freqcyc_age_sex_212325.csv`.
- `scripts/build_cycling_poststrat_rich1y.py`: builds `poststrat_freqcyc_age_sex_car_21.csv`.
- `scripts/build_validation_bicycle_to_work.py`: builds `validation_bicycle_to_work_share.csv`.
