# Build standardized non-rescue epilepsy medication tables for downstream analysis.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

source("src/desc/medication_standardization.R")
source("src/data_corrections.R")

output_dir <- "output/tabs/medication_standardization"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

analysis_date <- Sys.Date()

medications_raw <- read_medications_corrected()
medications_standardized <- standardize_non_rescue_epilepsy_medications(
  medications_raw,
  include_non_drug = FALSE
) %>%
  mutate(analysis_date = analysis_date)

med_intervals <- read_medication_intervals_corrected(
  medications_for_analysis = medications_standardized
)

active_medication_ids <- med_intervals %>%
  filter(
    start_date <= analysis_date,
    is.na(end_date) | end_date >= analysis_date
  ) %>%
  distinct(patient_id, medication_id) %>%
  mutate(is_active_on_analysis_date = TRUE)

medications_standardized <- medications_standardized %>%
  left_join(active_medication_ids, by = c("patient_id", "medication_id")) %>%
  mutate(is_active_on_analysis_date = tidyr::replace_na(.data$is_active_on_analysis_date, FALSE))

active_medications <- medications_standardized %>%
  filter(is_active_on_analysis_date) %>%
  select(
    analysis_date, patient_id, medication_id, name, name_standardized,
    standardized_components, standardization_status, requires_review, review_reason
  ) %>%
  arrange(patient_id, name_standardized, medication_id)

patient_active_medication_summary <- active_medications %>%
  group_by(patient_id) %>%
  summarise(
    n_active_medication_records = n_distinct(medication_id),
    n_active_standardized_names = n_distinct(name_standardized),
    active_medication_names_standardized = paste(sort(unique(name_standardized)), collapse = "; "),
    active_records_needing_review = sum(requires_review),
    .groups = "drop"
  ) %>%
  arrange(desc(n_active_standardized_names), patient_id)

unsure_medications <- medications_standardized %>%
  filter(requires_review) %>%
  group_by(name, name_standardized, standardized_components, standardization_status, review_reason) %>%
  summarise(
    n_records = n(),
    n_patients = n_distinct(patient_id),
    n_active_records = sum(is_active_on_analysis_date),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients), desc(n_records), name)

readr::write_csv(
  medications_standardized,
  file.path(output_dir, "non_rescue_epilepsy_medications_standardized.csv")
)
readr::write_csv(
  active_medications,
  file.path(output_dir, "non_rescue_epilepsy_active_medications_standardized.csv")
)
readr::write_csv(
  patient_active_medication_summary,
  file.path(output_dir, "non_rescue_epilepsy_active_medication_summary.csv")
)
readr::write_csv(
  unsure_medications,
  file.path(output_dir, "non_rescue_epilepsy_medications_unsure.csv")
)
