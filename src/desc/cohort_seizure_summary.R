# Cohort-level seizure summary.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

source("src/data_corrections.R")

PANEL_INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
SURVEY_INPUT_PATH <- "data/prospective_surveys.csv"
REGISTRY_INPUT_PATH <- "data/registry.csv"
CLUSTER_ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
SEIZURE_FEATURES_INPUT_PATH <- "output/tabs/clustering/seizure_freq_features.csv"
MANUAL_PATIENT_MAP_INPUT_PATH <- "data/manual_patient_map.csv"
OUTPUT_TAB_DIR <- "output/tabs"
OUTPUT_PATH <- file.path(OUTPUT_TAB_DIR, "summary_table.csv")
DAYS_PER_MONTH <- STANDARD_MONTH_DAYS

US_STATE_NAMES <- c(
  "alabama", "alaska", "arizona", "arkansas", "california", "colorado",
  "connecticut", "delaware", "florida", "georgia", "hawaii", "idaho",
  "illinois", "indiana", "iowa", "kansas", "kentucky", "louisiana",
  "maine", "maryland", "massachusetts", "michigan", "minnesota",
  "mississippi", "missouri", "montana", "nebraska", "nevada",
  "new hampshire", "new jersey", "new mexico", "new york",
  "north carolina", "north dakota", "ohio", "oklahoma", "oregon",
  "pennsylvania", "rhode island", "south carolina", "south dakota",
  "tennessee", "texas", "utah", "vermont", "virginia", "washington",
  "west virginia", "wisconsin", "wyoming", "district of columbia"
)

is_us_location <- function(location) {
  location_clean <- stringr::str_to_lower(stringr::str_squish(as.character(location)))
  location_clean <- dplyr::na_if(location_clean, "")
  state_pattern <- paste0(
    "\\b(",
    paste(stringr::str_replace_all(US_STATE_NAMES, " ", "\\\\s+"), collapse = "|"),
    ")\\b"
  )

  dplyr::coalesce(
    stringr::str_detect(location_clean, "\\b(united states|usa|u\\.s\\.a\\.|u\\.s\\.|us)\\b") |
      stringr::str_detect(location_clean, state_pattern),
    FALSE
  )
}

dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PANEL_INPUT_PATH)) {
  stop("Input file not found: ", PANEL_INPUT_PATH, call. = FALSE)
}

if (!file.exists(SURVEY_INPUT_PATH)) {
  stop("Input file not found: ", SURVEY_INPUT_PATH, call. = FALSE)
}

if (!file.exists(REGISTRY_INPUT_PATH)) {
  stop("Input file not found: ", REGISTRY_INPUT_PATH, call. = FALSE)
}

if (!file.exists(CLUSTER_ASSIGNMENTS_INPUT_PATH)) {
  stop("Input file not found: ", CLUSTER_ASSIGNMENTS_INPUT_PATH, call. = FALSE)
}

if (!file.exists(SEIZURE_FEATURES_INPUT_PATH)) {
  stop("Input file not found: ", SEIZURE_FEATURES_INPUT_PATH, call. = FALSE)
}

if (!file.exists(MANUAL_PATIENT_MAP_INPUT_PATH)) {
  stop("Input file not found: ", MANUAL_PATIENT_MAP_INPUT_PATH, call. = FALSE)
}

patient_month_panel <- readr::read_csv(PANEL_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    study_start_date = as.Date(.data$study_start_date),
    study_end_date = pmin(as.Date(.data$study_end_date), ANALYSIS_CUTOFF_DATE),
    month = as.Date(.data$month),
    seizure_count = suppressWarnings(as.numeric(.data$seizure_count)),
    month_end = as.Date(lubridate::ceiling_date(.data$month, "month") - lubridate::days(1)),
    observed_start = as.Date(pmax(.data$study_start_date, .data$month), origin = "1970-01-01"),
    observed_end = as.Date(pmin(.data$study_end_date, .data$month_end), origin = "1970-01-01"),
    observed_days = pmax(as.integer(.data$observed_end - .data$observed_start + 1L), 0L),
    observed_patient_months_28d = .data$observed_days / DAYS_PER_MONTH,
    seizures_per_28d_month = dplyr::if_else(
      .data$observed_days > 0,
      .data$seizure_count / .data$observed_days * DAYS_PER_MONTH,
      NA_real_
    )
  ) %>%
  dplyr::filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    !is.na(.data$patient_id),
    !is.na(.data$seizure_count),
    !is.na(.data$seizures_per_28d_month),
    .data$observed_days > 0
  )

if (dplyr::n_distinct(patient_month_panel$patient_id) != length(TARGET_PATIENT_IDS)) {
  stop(
    "Expected ",
    length(TARGET_PATIENT_IDS),
    " analysis patients from ",
    PANEL_INPUT_PATH,
    " but found ",
    dplyr::n_distinct(patient_month_panel$patient_id),
    ".",
    call. = FALSE
  )
}

completed_prospective_survey_patient_ids <- read_prospective_surveys_corrected(SURVEY_INPUT_PATH) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    prospective_study_complete = suppressWarnings(as.numeric(.data$prospective_study_complete))
  ) %>%
  dplyr::filter(
    .data$patient_id %in% unique(patient_month_panel$patient_id),
    .data$prospective_study_complete == 2
  ) %>%
  dplyr::distinct(.data$patient_id) %>%
  dplyr::pull("patient_id")

cluster_assignments <- readr::read_csv(CLUSTER_ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    pam_cluster = suppressWarnings(as.integer(.data$pam_k3))
  ) %>%
  dplyr::distinct(.data$patient_id, .keep_all = TRUE)

manual_patient_ids <- readr::read_csv(MANUAL_PATIENT_MAP_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::transmute(patient_id = as.character(.data$patient_id)) %>%
  dplyr::filter(.data$patient_id %in% unique(patient_month_panel$patient_id)) %>%
  dplyr::distinct() %>%
  dplyr::pull("patient_id")

cohort_levels <- c("Cluster 1", "Cluster 2", "Cluster 3", "Manual patients", "Total")

cluster_cohorts <- tibble::tibble(patient_id = unique(patient_month_panel$patient_id)) %>%
  dplyr::left_join(cluster_assignments, by = "patient_id") %>%
  dplyr::mutate(
    cohort = dplyr::if_else(
      .data$pam_cluster %in% 1:3,
      paste("Cluster", .data$pam_cluster),
      NA_character_
    )
  )

unassigned_patients <- cluster_cohorts %>%
  dplyr::filter(is.na(.data$cohort)) %>%
  dplyr::pull("patient_id")

if (length(unassigned_patients) > 0) {
  stop(
    "Cohort labels could not be assigned for patient(s): ",
    paste(unassigned_patients, collapse = ", "),
    ".",
    call. = FALSE
  )
}

patient_cohorts <- dplyr::bind_rows(
  cluster_cohorts %>% dplyr::select(patient_id, cohort),
  tibble::tibble(patient_id = manual_patient_ids, cohort = "Manual patients")
)

total_patient_cohort <- tibble::tibble(
  patient_id = unique(patient_month_panel$patient_id),
  cohort = "Total"
)

patient_month_cohort_panel <- patient_month_panel %>%
  dplyr::inner_join(patient_cohorts, by = "patient_id", relationship = "many-to-many")

registry_raw <- readr::read_csv(REGISTRY_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::mutate(patient_id = as.character(.data$patient_id))

required_registry_columns <- c("patient_id", "sex", "age_years", "location")
missing_registry_columns <- setdiff(required_registry_columns, names(registry_raw))
if (length(missing_registry_columns) > 0) {
  stop(
    "Required registry column(s) missing from ",
    REGISTRY_INPUT_PATH,
    ": ",
    paste(missing_registry_columns, collapse = ", "),
    call. = FALSE
  )
}

missing_registry_patients <- setdiff(unique(patient_month_panel$patient_id), registry_raw$patient_id)
if (length(missing_registry_patients) > 0) {
  warning(
    "Registry rows are missing for ",
    length(missing_registry_patients),
    " cohort patient(s): ",
    paste(missing_registry_patients, collapse = ", "),
    ". Registry summaries use available registry rows only.",
    call. = FALSE
  )
}

registry_patient_data <- registry_raw %>%
  dplyr::filter(.data$patient_id %in% unique(patient_month_panel$patient_id)) %>%
  dplyr::distinct(.data$patient_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    sex_clean = stringr::str_to_lower(stringr::str_squish(as.character(.data$sex))),
    age_years = suppressWarnings(as.numeric(.data$age_years)),
    location_clean = stringr::str_squish(as.character(.data$location)),
    location_clean = dplyr::na_if(.data$location_clean, ""),
    international_patient = !is_us_location(.data$location_clean) & !is.na(.data$location_clean)
  )

registry_summary <- dplyr::bind_rows(
  registry_patient_data %>% dplyr::inner_join(patient_cohorts, by = "patient_id"),
  registry_patient_data %>% dplyr::inner_join(total_patient_cohort, by = "patient_id")
) %>%
  dplyr::group_by(.data$cohort) %>%
  dplyr::summarise(
    sex_male = as.integer(sum(.data$sex_clean == "male", na.rm = TRUE)),
    international_patients = as.integer(sum(.data$international_patient, na.rm = TRUE)),
    mean_age = mean(.data$age_years, na.rm = TRUE),
    .groups = "drop"
  )

survey_completion_summary <- dplyr::bind_rows(
  patient_cohorts,
  total_patient_cohort
) %>%
  dplyr::group_by(.data$cohort) %>%
  dplyr::summarise(
    patients_completed_prospective_survey = as.integer(sum(
      .data$patient_id %in% completed_prospective_survey_patient_ids
    )),
    .groups = "drop"
  )

patient_follow_up <- dplyr::bind_rows(
  patient_month_cohort_panel,
  patient_month_panel %>% dplyr::inner_join(total_patient_cohort, by = "patient_id")
) %>%
  dplyr::group_by(.data$cohort, .data$patient_id) %>%
  dplyr::summarise(
    patient_months_28d = sum(.data$observed_patient_months_28d),
    .groups = "drop"
  )

follow_up_summary <- patient_follow_up %>%
  dplyr::group_by(.data$cohort) %>%
  dplyr::summarise(
    min_patient_months_28d = min(.data$patient_months_28d),
    max_patient_months_28d = max(.data$patient_months_28d),
    .groups = "drop"
  )

patient_rate_data <- readr::read_csv(SEIZURE_FEATURES_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    patient_mean_seizures_per_28d_month = suppressWarnings(as.numeric(.data$mean_monthly_seizure_rate)),
    patient_zero_month_proportion = suppressWarnings(as.numeric(.data$proportion_zero_seizure_months))
  ) %>%
  dplyr::filter(
    .data$patient_id %in% unique(patient_month_panel$patient_id),
    is.finite(.data$patient_mean_seizures_per_28d_month),
    is.finite(.data$patient_zero_month_proportion)
  )

missing_rate_patients <- setdiff(unique(patient_month_panel$patient_id), patient_rate_data$patient_id)
if (length(missing_rate_patients) > 0) {
  stop(
    "Patient-level seizure-rate features are missing for patient(s): ",
    paste(missing_rate_patients, collapse = ", "),
    ".",
    call. = FALSE
  )
}

patient_rate_summary <- dplyr::bind_rows(
  patient_rate_data %>% dplyr::inner_join(patient_cohorts, by = "patient_id"),
  patient_rate_data %>% dplyr::inner_join(total_patient_cohort, by = "patient_id")
) %>%
  dplyr::group_by(.data$cohort) %>%
  dplyr::summarise(
    median_patient_mean_seizures_per_28d_month = stats::median(.data$patient_mean_seizures_per_28d_month),
    q1_patient_mean_seizures_per_28d_month = as.numeric(stats::quantile(
      .data$patient_mean_seizures_per_28d_month,
      probs = 0.25,
      names = FALSE
    )),
    q3_patient_mean_seizures_per_28d_month = as.numeric(stats::quantile(
      .data$patient_mean_seizures_per_28d_month,
      probs = 0.75,
      names = FALSE
    )),
    median_patient_zero_month_proportion = stats::median(.data$patient_zero_month_proportion),
    .groups = "drop"
  )

summary_table <- dplyr::bind_rows(
  patient_month_cohort_panel,
  patient_month_panel %>% dplyr::inner_join(total_patient_cohort, by = "patient_id")
) %>%
  dplyr::group_by(.data$cohort) %>%
  dplyr::summarise(
    n = as.integer(dplyr::n_distinct(.data$patient_id)),
    observed_patient_months_28d = sum(.data$observed_patient_months_28d),
    seizure_events = as.integer(sum(.data$seizure_count)),
    .groups = "drop"
  ) %>%
  dplyr::left_join(registry_summary, by = "cohort") %>%
  dplyr::left_join(survey_completion_summary, by = "cohort") %>%
  dplyr::left_join(follow_up_summary, by = "cohort") %>%
  dplyr::left_join(patient_rate_summary, by = "cohort") %>%
  dplyr::mutate(
    cohort = factor(.data$cohort, levels = cohort_levels),
    sex_male_display = sprintf(
      "%d (%.1f%%)",
      .data$sex_male,
      100 * .data$sex_male / .data$n
    ),
    international_patients_display = sprintf(
      "%d (%.1f%%)",
      .data$international_patients,
      100 * .data$international_patients / .data$n
    ),
    mean_age = round(.data$mean_age, 1),
    observed_patient_months_28d = round(.data$observed_patient_months_28d, 1),
    prospective_survey_display = sprintf(
      "%d (%.1f%%)",
      .data$patients_completed_prospective_survey,
      100 * .data$patients_completed_prospective_survey / .data$n
    ),
    seizures_per_28d_month_display = sprintf(
      "%.1f (%.1f–%.1f)",
      .data$median_patient_mean_seizures_per_28d_month,
      .data$q1_patient_mean_seizures_per_28d_month,
      .data$q3_patient_mean_seizures_per_28d_month
    ),
    zero_months_display = sprintf(
      "%.1f%%",
      100 * round(.data$median_patient_zero_month_proportion, 2)
    ),
    patient_observed_months_range = sprintf(
      "%.1f–%.1f",
      .data$min_patient_months_28d,
      .data$max_patient_months_28d
    )
  ) %>%
  dplyr::arrange(.data$cohort) %>%
  dplyr::select(
    cohort,
    n,
    sex_male_display,
    international_patients_display,
    mean_age,
    observed_patient_months_28d,
    seizure_events,
    prospective_survey_display,
    seizures_per_28d_month_display,
    zero_months_display,
    patient_observed_months_range
  ) %>%
  dplyr::rename(
    Cohort = cohort,
    `Sex (male), n (%)` = sex_male_display,
    `International patients, n (%)` = international_patients_display,
    `Mean age, years` = mean_age,
    `Observed patient-months (28-day)` = observed_patient_months_28d,
    `Seizure events` = seizure_events,
    `Completed prospective survey, n (%)` = prospective_survey_display,
    `Median patient mean seizures per 28-day month (IQR)` = seizures_per_28d_month_display,
    `Median patient zero-seizure months, %` = zero_months_display,
    `Patient observed months (28-day), min–max` = patient_observed_months_range
  )

readr::write_csv(summary_table, OUTPUT_PATH)
