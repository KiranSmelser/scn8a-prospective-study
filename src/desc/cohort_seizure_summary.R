# Cohort-level seizure summary.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

source("src/analysis_config.R")

PANEL_INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
SURVEY_INPUT_PATH <- "data/prospective_surveys.csv"
OUTPUT_TAB_DIR <- "output/tabs/seizure_patterns"
OUTPUT_PATH <- file.path(OUTPUT_TAB_DIR, "cohort_seizure_summary.csv")
DAYS_PER_MONTH <- STANDARD_MONTH_DAYS

dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PANEL_INPUT_PATH)) {
  stop("Input file not found: ", PANEL_INPUT_PATH, call. = FALSE)
}

if (!file.exists(SURVEY_INPUT_PATH)) {
  stop("Input file not found: ", SURVEY_INPUT_PATH, call. = FALSE)
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

completed_prospective_survey_patients <- readr::read_csv(SURVEY_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    prospective_study_complete = suppressWarnings(as.numeric(.data$prospective_study_complete))
  ) %>%
  dplyr::filter(
    .data$patient_id %in% unique(patient_month_panel$patient_id),
    .data$prospective_study_complete == 2
  ) %>%
  dplyr::distinct(.data$patient_id) %>%
  nrow()

patient_follow_up <- patient_month_panel %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(
    patient_months_28d = sum(.data$observed_patient_months_28d),
    .groups = "drop"
  )

summary_table <- patient_month_panel %>%
  dplyr::summarise(
    n = as.integer(dplyr::n_distinct(.data$patient_id)),
    observed_patient_months_28d = sum(.data$observed_patient_months_28d),
    seizure_events = as.integer(sum(.data$seizure_count)),
    patients_completed_prospective_survey = as.integer(completed_prospective_survey_patients),
    mean_seizures_per_28d_month = mean(.data$seizures_per_28d_month),
    median_seizures_per_28d_month = stats::median(.data$seizures_per_28d_month),
    sd_seizures_per_28d_month = stats::sd(.data$seizures_per_28d_month),
    min_seizures_per_28d_month = min(.data$seizures_per_28d_month),
    max_seizures_per_28d_month = max(.data$seizures_per_28d_month),
    q1_seizures_per_28d_month = as.numeric(stats::quantile(.data$seizures_per_28d_month, probs = 0.25, names = FALSE)),
    q3_seizures_per_28d_month = as.numeric(stats::quantile(.data$seizures_per_28d_month, probs = 0.75, names = FALSE)),
    zero_months = as.integer(sum(.data$seizure_count == 0)),
    min_patient_months_28d = min(patient_follow_up$patient_months_28d),
    median_patient_months_28d = stats::median(patient_follow_up$patient_months_28d),
    max_patient_months_28d = max(patient_follow_up$patient_months_28d)
  )

readr::write_csv(summary_table, OUTPUT_PATH)
