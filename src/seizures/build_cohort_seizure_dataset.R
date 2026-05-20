# Build a cleaned seizure-event dataset for patients with high compliance.

suppressPackageStartupMessages({
  library(tidyverse)
})

source("src/desc/seizure_type_standardization.R")

TARGET_PATIENT_IDS <- c(
  "67b6126de8efcd6464eeb8c4",
  "67b59b49e8efcd6464eeb595",
  "690a3a45b2c2dee67c5423d6",
  "68d58414b783c7baf7292504",
  "691b83cb7db39ba79465d9cf",
  "6808969d5413b5941d3d58d8",
  "67fd6a03773a9ee7d191cbea",
  "67e82887773a9ee7d190fb0a",
  "5dcb3f4d559dc95b43af4872",
  "67ead18b773a9ee7d1911a39",
  "67d181092099f5b1ffcbc0db",
  "67b63d09e8efcd6464eebcba",
  "6836626b1191a9613d19ac21",
  "67d0c5202099f5b1ffcbbbc9",
  "67a8e3f7890c36004c181796",
  "67b529afe8efcd6464eeb287",
  "67e5f605773a9ee7d190e931",
  "6954f548ce6a468b362fce7b",
  "67b60d2ee8efcd6464eeb882",
  "688de8e1e20e57ca0b5335eb",
  "681876911191a9613d188b77",
  "67b562cae8efcd6464eeb345",
  "693e0b36ce6a468b362ee35d"
)

PATIENT_START_DATES <- tibble::tibble(
  patient_id = c(
    "67b59b49e8efcd6464eeb595",
    "690a3a45b2c2dee67c5423d6",
    "68d58414b783c7baf7292504",
    "691b83cb7db39ba79465d9cf",
    "6808969d5413b5941d3d58d8"
  ),
  min_event_date = as.Date(c(
    "2025-08-28",
    "2025-12-26",
    "2025-10-04",
    "2025-12-14",
    "2025-05-13"
  ))
)

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

cohort_events <- readr::read_csv("data/events.csv", show_col_types = FALSE) %>%
  standardize_seizure_events(filter_to_seizure = TRUE) %>%
  dplyr::mutate(patient_id = as.character(.data$patient_id)) %>%
  dplyr::filter(.data$patient_id %in% TARGET_PATIENT_IDS) %>%
  dplyr::left_join(PATIENT_START_DATES, by = "patient_id") %>%
  dplyr::filter(is.na(.data$min_event_date) | (!is.na(.data$event_date) & .data$event_date >= .data$min_event_date)) %>%
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
