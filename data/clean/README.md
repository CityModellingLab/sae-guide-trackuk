# Clean Data

Tutorial-ready files derived from the raw survey, lookup, boundary, Census, and population sources.

## Lookups

- `msoa_lad_lookup.csv`: Wales MSOA-to-local-authority lookup used by the clean-data scripts.
- `dvla_lad_lookup.csv`: National Survey for Wales `DvLA` code to current local-authority code bridge.

## Survey Files

- `survey_microdata_msoayear.csv`: one row per respondent-wave. Includes harmonised survey fields, synthetic MSOA assignment, walking indicator, weights, and respondent predictors.
- `surveyind_freqwalk.csv`: one row per respondent-wave for the walking estimand. Includes `unit_id`, `time_id`, `area_id`, `strata`, `indicator`, `weight`, and unit predictors `sex`, `age`, and `car_access`.
- `domain_msoayear.csv`: complete MSOA-year estimation grid, including cells with no survey respondents.

## Walking Auxiliary and Validation Files

- `areapred_freqwalk.csv`: one row per MSOA with raw and z-scored area predictors for daily walking models, including Census 2021 mean age among usual residents aged 16+ and proportion female from Nomis bulk tables.
- `poststrat_freqwalk_age_sex_212325.csv`: one row per MSOA-year-age-sex cell with `pop_count` for shallow poststratification. Age bands are grouped to match the rich poststratification table.
- `poststrat_freqwalk_age_sex_car_21.csv`: one row per MSOA-year-age-sex-car-access cell with `pop_count` for richer poststratification.
- `validation_adult_walk_share.csv`: one row per Wales MSOA with a simple 2021 walk-to-work benchmark. The numerator is the TS061 walk-to-work count and the denominator is the TS061 total usual residents aged 16+ in employment.
- `validation_ts061_walk_to_work.csv`: one row per Wales MSOA with Census 2021 TS061 walking-to-work, cycling-to-work, and work-from-home counts and shares. Use `ts061_walk_to_work_share` as an external benchmark for spatial pattern checks only: it covers usual residents aged 16+ in employment on Census Day 2021, not all adults walking for transport daily.

## Notes

- Updated MSOA codes with mismatch from survey microdata.
    - `DvLA == 7` maps to Powys, `W06000023`.
    - `DvLA == 17` maps to Merthyr Tydfil, `W06000024`.
- The 2024-25 raw survey file does not contain exact respondent `Age`. The clean-data builder derives a synthetic numeric `age` from `DvAgeGrp7` by sampling from the 2021 observed age mix within each age-group interval, with deterministic random jitter and a cap of 95 for the open-ended 75+ group. The `survey_microdata_msoayear.csv` file records this in `age_source`.
- The shallow 2025 poststratification table uses the latest supplied mid-2024 population denominator, as the mid-2025 estimates have not been published as of July 2026.
- In the shallow poststratification table, the youngest available ONS age band is `15 to 19`, so age 15+ is used as the closest available match to the NSW adult definition of age 16+.

## Build Scripts

- `scripts/build_clean_survey_microdata.py`: builds `survey_microdata_msoayear.csv`, `surveyind_freqwalk.csv`, and `domain_msoayear.csv`.
- `scripts/build_walking_area_predictors.py`: builds `areapred_freqwalk.csv`.
- `scripts/build_walking_poststrat_shallow3y.py`: builds `poststrat_freqwalk_age_sex_212325.csv`.
- `scripts/build_walking_poststrat_rich1y.py`: builds `poststrat_freqwalk_age_sex_car_21.csv`.
- `scripts/build_validation_adult_walk_share.py`: builds `validation_adult_walk_share.csv`.
