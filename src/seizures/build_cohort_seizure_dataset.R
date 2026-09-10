# Build a cleaned seizure-event dataset for patients with high compliance.

suppressPackageStartupMessages({
  library(tidyverse)
})

source("src/desc/seizure_type_standardization.R")
source("src/analysis_config.R")

OUTPUT_TAB_DIR <- "output/tabs/seizures"
OUTPUT_PATH <- file.path(OUTPUT_TAB_DIR, "seizures.csv")

dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

sanitize_text <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_squish() %>%
    dplyr::na_if("") %>%
    dplyr::if_else(
      stringr::str_to_lower(dplyr::coalesce(., "")) %in% c("na", "n/a", "nan", "null"),
      NA_character_,
      .
    )
}

patient_metadata <- readr::read_csv("data/whatsapp_status.csv", show_col_types = FALSE) %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    first_name = sanitize_text(.data$first_name),
    last_name = sanitize_text(.data$last_name),
    variant_p = sanitize_text(.data$variant_p),
    variant_c = sanitize_text(.data$variant_c),
    run_timestamp = as.character(.data$run_timestamp)
  ) %>%
  dplyr::filter(!is.na(.data$patient_id), .data$patient_id != "") %>%
  dplyr::arrange(.data$patient_id, dplyr::desc(.data$run_timestamp)) %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(
    first_name = dplyr::first(.data$first_name),
    last_name = dplyr::first(.data$last_name),
    variant_p = dplyr::first(.data$variant_p),
    variant_c = dplyr::first(.data$variant_c),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    patient_name = stringr::str_squish(stringr::str_trim(paste(.data$first_name, .data$last_name))),
    patient_name = dplyr::na_if(.data$patient_name, "")
  )

cohort_events <- read_events_corrected() %>%
  standardize_seizure_events(filter_to_seizure = TRUE) %>%
  dplyr::mutate(patient_id = as.character(.data$patient_id)) %>%
  dplyr::filter(.data$patient_id %in% TARGET_PATIENT_IDS) %>%
  dplyr::left_join(PATIENT_START_DATES, by = "patient_id") %>%
  dplyr::filter(
    !is.na(.data$event_date),
    .data$event_date <= ANALYSIS_CUTOFF_DATE,
    is.na(.data$min_event_date) | .data$event_date >= .data$min_event_date
  ) %>%
  dplyr::left_join(patient_metadata, by = "patient_id") %>%
  dplyr::arrange(.data$patient_id, .data$event_date, .data$date, .data$event_id) %>%
  dplyr::transmute(
    event_id = .data$event_id,
    patient_id = .data$patient_id,
    type = .data$seizure_group,
    date = .data$date,
    duration = .data$duration,
    triggers = .data$triggers,
    during_sleep = .data$during_sleep
  )

readr::write_csv(cohort_events, OUTPUT_PATH)

missing_patients <- setdiff(TARGET_PATIENT_IDS, unique(cohort_events$patient_id))
