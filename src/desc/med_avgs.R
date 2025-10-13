# Average medications per month per patient

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

target_year <- 2025L
current_date <- Sys.Date()
target_start <- as.Date(sprintf("%d-01-01", target_year))
target_end_cap <- as.Date(sprintf("%d-12-31", target_year))
analysis_end <- min(target_end_cap, current_date)

# Load medication data
medications <- readr::read_csv("data/medications.csv", show_col_types = FALSE) %>%
  dplyr::select(medication_id, med_patient_id = .data$patient_id, medication_name = .data$name)

# Load dosage data
med_dosages_raw <- readr::read_csv("data/med_dosages.csv", show_col_types = FALSE) %>%
  dplyr::mutate(
    start_date = as.Date(.data$from),
    end_date = as.Date(.data$to)
  ) %>%
  dplyr::mutate(
    end_date = dplyr::if_else(is.na(.data$end_date), analysis_end, .data$end_date)
  ) %>%
  dplyr::filter(!is.na(.data$start_date))

# Derive distinct medication intervals
dosage_intervals <- med_dosages_raw %>%
  dplyr::left_join(medications, by = "medication_id") %>%
  dplyr::filter(is.na(.data$med_patient_id) | .data$med_patient_id == .data$patient_id) %>%
  dplyr::mutate(patient_id = dplyr::coalesce(.data$patient_id, .data$med_patient_id)) %>%
  dplyr::filter(!is.na(.data$patient_id)) %>%
  dplyr::select(.data$patient_id, .data$medication_id, .data$start_date, .data$end_date) %>%
  dplyr::filter(.data$start_date <= analysis_end, .data$end_date >= target_start) %>%
  dplyr::mutate(
    start_date = dplyr::if_else(.data$start_date < target_start, target_start, .data$start_date),
    end_date = dplyr::if_else(.data$end_date > analysis_end, analysis_end, .data$end_date)
  ) %>%
  dplyr::filter(.data$start_date <= .data$end_date) %>%
  dplyr::distinct()

# Create month grid up to analysis end date
if (analysis_end < target_start) {
  months <- tibble::tibble(
    month = as.Date(character()),
    month_end = as.Date(character())
  )
} else {
  month_seq <- seq.Date(from = target_start, to = floor_date(analysis_end, "month"), by = "month")
  months <- tibble::tibble(month = month_seq) %>%
    dplyr::mutate(
      month_end = pmin((.data$month %m+% months(1)) - days(1), analysis_end)
    )
}

# Determine active medications
patient_month_med <- dosage_intervals %>%
  tidyr::expand_grid(months) %>%
  dplyr::filter(.data$start_date <= .data$month_end, .data$end_date >= .data$month) %>%
  dplyr::select(.data$patient_id, .data$medication_id, .data$month) %>%
  dplyr::distinct()

monthly_med_counts <- patient_month_med %>%
  dplyr::count(.data$patient_id, .data$month, name = "active_medications")

# Fill in months with zero medications
monthly_med_counts_complete <- monthly_med_counts %>%
  dplyr::group_by(.data$patient_id) %>%
  tidyr::complete(month = months$month, fill = list(active_medications = 0)) %>%
  dplyr::ungroup()

# Summarize
unique_medications <- patient_month_med %>%
  dplyr::distinct(.data$patient_id, .data$medication_id) %>%
  dplyr::count(.data$patient_id, name = "unique_medications")

med_avgs <- monthly_med_counts_complete %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(
    avg_medications_per_month = mean(.data$active_medications),
    max_concurrent_medications = max(.data$active_medications),
    months_with_medications = sum(.data$active_medications > 0),
    months_recorded = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::left_join(unique_medications, by = "patient_id") %>%
  dplyr::mutate(unique_medications = dplyr::coalesce(.data$unique_medications, 0L)) %>%
  dplyr::arrange(dplyr::desc(.data$avg_medications_per_month))

if (!dir.exists("output/tabs")) {
  dir.create("output/tabs", recursive = TRUE, showWarnings = FALSE)
}

write.csv(med_avgs, "output/tabs/med_avgs.csv", row.names = FALSE)
