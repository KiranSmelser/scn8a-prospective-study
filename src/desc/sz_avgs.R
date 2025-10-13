# Average seizures per month per patient

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

target_year <- 2025L
analysis_start <- as.Date(sprintf("%d-01-01", target_year))
analysis_end_cap <- as.Date(sprintf("%d-12-31", target_year))
current_date <- Sys.Date()
analysis_end <- min(analysis_end_cap, current_date)
analysis_months <- if (analysis_end < analysis_start) {
  as.Date(character())
} else {
  seq.Date(from = analysis_start, to = floor_date(analysis_end, "month"), by = "month")
}

events <- readr::read_csv("data/events.csv", show_col_types = FALSE) %>%
  dplyr::filter(tolower(.data$type) == "seizure") %>%
  dplyr::mutate(
    date_time = lubridate::ymd_hms(.data$date, quiet = TRUE, tz = "UTC"),
    event_date = as.Date(.data$date_time),
    year = lubridate::year(.data$date_time),
    month = as.Date(lubridate::floor_date(.data$date_time, "month"))
  ) %>%
  dplyr::filter(
    .data$year == target_year,
    !is.na(.data$patient_id),
    !is.na(.data$month),
    .data$event_date <= analysis_end
  )

monthly_counts <- events %>%
  dplyr::group_by(.data$patient_id, .data$month) %>%
  dplyr::summarise(seizure_count = dplyr::n(), .groups = "drop")

if (length(analysis_months) > 0) {
  monthly_counts <- monthly_counts %>%
    dplyr::filter(.data$month %in% analysis_months)
}

monthly_counts_complete <- monthly_counts %>%
  dplyr::group_by(.data$patient_id) %>%
  tidyr::complete(month = analysis_months, fill = list(seizure_count = 0)) %>%
  dplyr::ungroup()

if (length(analysis_months) == 0) {
  sz_avgs <- tibble::tibble(
    patient_id = character(),
    avg_seizures_per_month = numeric(),
    total_seizures = integer(),
    months_with_seizures = integer(),
    months_recorded = integer()
  )
} else {
  sz_avgs <- monthly_counts_complete %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::summarise(
      avg_seizures_per_month = mean(.data$seizure_count),
      total_seizures = sum(.data$seizure_count),
      months_with_seizures = sum(.data$seizure_count > 0),
      months_recorded = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(.data$avg_seizures_per_month))
}

if (!dir.exists("output/tabs")) {
  dir.create("output/tabs", recursive = TRUE, showWarnings = FALSE)
}

write.csv(sz_avgs, "output/tabs/sz_avgs.csv", row.names = FALSE)
