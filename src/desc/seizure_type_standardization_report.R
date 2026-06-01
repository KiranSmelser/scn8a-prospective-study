# Build standardized seizure-type tables for downstream analysis.

suppressPackageStartupMessages({
  library(tidyverse)
})

source("src/desc/seizure_type_standardization.R")

output_dir <- "output/tabs/seizure_standardization"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_date <- analysis_end_date()

events_raw <- readr::read_csv("data/events.csv", show_col_types = FALSE)

seizure_events_standardized <- standardize_seizure_events(events_raw) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    analysis_date = analysis_date
  ) %>%
  dplyr::filter(.data$patient_id %in% TARGET_PATIENT_IDS)

seizure_type_reference <- seizure_events_standardized %>%
  dplyr::group_by(
    .data$seizure_type_raw,
    .data$seizure_type_standardized,
    .data$seizure_type_primary,
    .data$seizure_group,
    .data$standardization_status,
    .data$requires_review,
    .data$review_reason
  ) %>%
  dplyr::summarise(
    n_events = dplyr::n(),
    n_patients = dplyr::n_distinct(.data$patient_id),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(.data$n_events), desc(.data$n_patients), .data$seizure_type_raw)

unsure_seizure_types <- seizure_type_reference %>%
  dplyr::filter(.data$requires_review) %>%
  dplyr::arrange(desc(.data$n_events), desc(.data$n_patients), .data$seizure_type_raw)

readr::write_csv(
  seizure_events_standardized,
  file.path(output_dir, "seizure_events_standardized.csv")
)
readr::write_csv(
  seizure_type_reference,
  file.path(output_dir, "seizure_type_mapping_reference.csv")
)
readr::write_csv(
  unsure_seizure_types,
  file.path(output_dir, "seizure_types_unsure.csv")
)
