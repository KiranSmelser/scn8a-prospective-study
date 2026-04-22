# Average seizures per month per patient

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(jsonlite)
})

current_date <- Sys.Date()
analysis_end <- current_date
analysis_month_cap <- if (is.na(analysis_end)) NA else floor_date(analysis_end, "month")

events <- readr::read_csv("data/events.csv", show_col_types = FALSE) %>%
  dplyr::filter(tolower(.data$type) == "seizure") %>%
  dplyr::mutate(
    date_time = lubridate::ymd_hms(.data$date, quiet = TRUE, tz = "UTC"),
    event_date = as.Date(.data$date_time),
    month = as.Date(lubridate::floor_date(.data$date_time, "month"))
  ) %>%
  dplyr::filter(
    !is.na(.data$patient_id),
    !is.na(.data$month),
    .data$event_date <= analysis_end
  )

usage_records <- if (file.exists("output/tabs/app_usage.json")) {
  jsonlite::fromJSON("output/tabs/app_usage.json", simplifyVector = FALSE)
} else {
  list()
}

usage_months <- purrr::map_dfr(
  usage_records,
  function(rec) {
    patient <- rec$patient_id
    usage <- rec$app_usage_per_month
    if (is.null(patient) || is.null(usage) || length(usage) == 0) {
      return(tibble::tibble())
    }
    usage_vec <- unlist(usage, use.names = TRUE)
    if (is.null(usage_vec) || length(usage_vec) == 0) {
      return(tibble::tibble())
    }
    tibble::tibble(
      patient_id = patient,
      month = names(usage_vec),
      app_usage = as.numeric(usage_vec)
    )
  }
)

usage_months <- usage_months %>%
  dplyr::mutate(
    month = as.Date(paste0(.data$month, "-01")),
    month = lubridate::floor_date(.data$month, "month")
  )

if (!is.na(analysis_month_cap)) {
  usage_months <- usage_months %>%
    dplyr::filter(
      !is.na(.data$patient_id),
      !is.na(.data$month),
      .data$month <= analysis_month_cap,
      .data$app_usage >= 1
    ) %>%
    dplyr::distinct(.data$patient_id, .data$month)
} else {
  usage_months <- tibble::tibble()
}

if (nrow(usage_months) == 0) {
  sz_avgs <- tibble::tibble(
    patient_id = character(),
    avg_seizures_per_month = numeric(),
    total_seizures = integer(),
    months_with_seizures = integer(),
    months_recorded = integer()
  )
} else {
  monthly_counts <- events %>%
    dplyr::semi_join(usage_months, by = c("patient_id", "month")) %>%
    dplyr::group_by(.data$patient_id, .data$month) %>%
    dplyr::summarise(seizure_count = dplyr::n(), .groups = "drop")

  monthly_counts_complete <- usage_months %>%
    dplyr::left_join(monthly_counts, by = c("patient_id", "month")) %>%
    tidyr::replace_na(list(seizure_count = 0))

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

OUTPUT_TAB_DIR <- "output/tabs/seizure_patterns"
if (!dir.exists(OUTPUT_TAB_DIR)) {
  dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
}

write.csv(sz_avgs, file.path(OUTPUT_TAB_DIR, "sz_avgs.csv"), row.names = FALSE)
