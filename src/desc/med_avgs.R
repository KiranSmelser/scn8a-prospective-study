# Average medications per month per patient (usage-filtered)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(jsonlite)
})

current_date <- Sys.Date()
analysis_end <- current_date

collapse_medication_names <- function(values) {
  cleaned <- trimws(values)
  cleaned <- cleaned[!is.na(cleaned)]
  cleaned <- cleaned[nzchar(cleaned)]
  cleaned <- sort(unique(cleaned))
  if (length(cleaned) == 0) {
    NA_character_
  } else {
    paste(cleaned, collapse = "; ")
  }
}

# Load medication data (restrict to epilepsy treatments)
medications <- readr::read_csv("data/medications.csv", show_col_types = FALSE) %>%
  dplyr::mutate(
    reason_clean = tolower(trimws(.data$reason))
  ) %>%
  dplyr::filter(!is.na(.data$reason_clean), .data$reason_clean == "epilepsy") %>%
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
  dplyr::filter(!is.na(.data$med_patient_id), .data$med_patient_id == .data$patient_id) %>%
  dplyr::mutate(patient_id = dplyr::coalesce(.data$patient_id, .data$med_patient_id)) %>%
  dplyr::filter(!is.na(.data$patient_id)) %>%
  dplyr::select(.data$patient_id, .data$medication_id, .data$start_date, .data$end_date) %>%
  dplyr::filter(.data$start_date <= analysis_end) %>%
  dplyr::mutate(
    end_date = dplyr::if_else(.data$end_date > analysis_end, analysis_end, .data$end_date)
  ) %>%
  dplyr::filter(.data$start_date <= .data$end_date) %>%
  dplyr::distinct()

# Identify months with app usage >= 1
analysis_month_cap <- if (is.na(analysis_end)) NA else lubridate::floor_date(analysis_end, "month")

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

usage_months <- if (!is.na(analysis_month_cap)) {
  usage_months %>%
    dplyr::filter(
      !is.na(.data$patient_id),
      !is.na(.data$month),
      .data$month <= analysis_month_cap,
      .data$app_usage >= 1
    ) %>%
    dplyr::mutate(
      month_end = pmin((.data$month %m+% months(1)) - days(1), analysis_end)
    ) %>%
    dplyr::select(.data$patient_id, .data$month, .data$month_end) %>%
    dplyr::distinct()
} else {
  tibble::tibble(
    patient_id = character(),
    month = as.Date(character()),
    month_end = as.Date(character())
  )
}

patient_med_months <- dosage_intervals %>%
  dplyr::mutate(
    month_start = lubridate::floor_date(.data$start_date, "month"),
    month_finish = lubridate::floor_date(.data$end_date, "month"),
    month = purrr::map2(.data$month_start, .data$month_finish, ~ seq(.x, .y, by = "month"))
  ) %>%
  tidyr::unnest(.data$month) %>%
  dplyr::mutate(
    month = as.Date(.data$month),
    month_end = pmin((.data$month %m+% months(1)) - days(1), .data$end_date, analysis_end)
  ) %>%
  dplyr::filter(.data$start_date <= .data$month_end, .data$end_date >= .data$month) %>%
  dplyr::select(.data$patient_id, .data$medication_id, .data$month, .data$month_end) %>%
  dplyr::distinct()

analysis_months <- dplyr::bind_rows(
  usage_months,
  patient_med_months %>% dplyr::select(.data$patient_id, .data$month, .data$month_end)
) %>%
  dplyr::group_by(.data$patient_id, .data$month) %>%
  dplyr::summarise(
    month_end = max(.data$month_end),
    .groups = "drop"
  )

if (nrow(analysis_months) == 0) {
  med_avgs <- tibble::tibble(
    patient_id = character(),
    avg_medications_per_month = numeric(),
    max_concurrent_medications = numeric(),
    months_with_medications = integer(),
    months_recorded = integer(),
    unique_medications = integer(),
    current_medications = integer(),
    current_medication_names = character(),
    weaned_medication_names = character()
  )
} else {
  patient_month_med <- patient_med_months %>%
    dplyr::select(.data$patient_id, .data$medication_id, .data$month) %>%
    dplyr::distinct()

  monthly_med_counts <- patient_month_med %>%
    dplyr::count(.data$patient_id, .data$month, name = "active_medications")

  monthly_med_counts_complete <- analysis_months %>%
    dplyr::left_join(monthly_med_counts, by = c("patient_id", "month")) %>%
    tidyr::replace_na(list(active_medications = 0))

  medication_lookup <- medications %>%
    dplyr::select(.data$medication_id, .data$medication_name) %>%
    dplyr::distinct()

  patient_med_history <- dosage_intervals %>%
    dplyr::group_by(.data$patient_id, .data$medication_id) %>%
    dplyr::summarise(
      is_current = any(.data$end_date == analysis_end),
      was_weaned = any(.data$end_date < analysis_end),
      .groups = "drop"
    ) %>%
    dplyr::left_join(medication_lookup, by = "medication_id") %>%
    dplyr::mutate(
      medication_display_name = dplyr::case_when(
        !is.na(.data$medication_name) & nzchar(trimws(.data$medication_name)) ~ trimws(.data$medication_name),
        TRUE ~ paste0("Unknown medication (", .data$medication_id, ")")
      )
    )

  current_medication_details <- patient_med_history %>%
    dplyr::filter(.data$is_current) %>%
    dplyr::distinct(.data$patient_id, .data$medication_id, .data$medication_display_name) %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::summarise(
      current_medication_list = list(sort(unique(.data$medication_display_name))),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      current_medications = lengths(.data$current_medication_list),
      current_medication_names = purrr::map_chr(
        .data$current_medication_list,
        collapse_medication_names
      )
    ) %>%
    dplyr::select(-.data$current_medication_list)

  weaned_medications <- patient_med_history %>%
    dplyr::filter(!.data$is_current & .data$was_weaned) %>%
    dplyr::distinct(.data$patient_id, .data$medication_id, .data$medication_display_name) %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::summarise(
      weaned_medication_names = collapse_medication_names(.data$medication_display_name),
      .groups = "drop"
    )

  overall_unique_medications <- patient_med_history %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::summarise(
      unique_medications = dplyr::n_distinct(.data$medication_id),
      .groups = "drop"
    )

  max_concurrent_overall <- dosage_intervals %>%
    dplyr::mutate(
      active_day = purrr::map2(.data$start_date, .data$end_date, seq, by = "day")
    ) %>%
    tidyr::unnest(.data$active_day) %>%
    dplyr::group_by(.data$patient_id, .data$active_day) %>%
    dplyr::summarise(
      active_medications = dplyr::n_distinct(.data$medication_id),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::summarise(
      max_concurrent_overall = max(.data$active_medications),
      .groups = "drop"
    )

  med_avgs <- monthly_med_counts_complete %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::summarise(
      avg_medications_per_month = mean(.data$active_medications),
      max_distinct_meds_in_month = max(.data$active_medications),
      months_with_medications = sum(.data$active_medications > 0),
      months_recorded = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::left_join(overall_unique_medications, by = "patient_id") %>%
    dplyr::left_join(current_medication_details, by = "patient_id") %>%
    dplyr::left_join(weaned_medications, by = "patient_id") %>%
    dplyr::left_join(max_concurrent_overall, by = "patient_id") %>%
    dplyr::mutate(
      max_concurrent_medications = dplyr::coalesce(.data$max_concurrent_overall, 0L),
      unique_medications = dplyr::coalesce(.data$unique_medications, 0L),
      current_medications = dplyr::coalesce(.data$current_medications, 0L)
    ) %>%
    dplyr::select(dplyr::all_of(c(
      "patient_id",
      "avg_medications_per_month",
      "max_concurrent_medications",
      "months_with_medications",
      "months_recorded",
      "unique_medications",
      "current_medications",
      "current_medication_names",
      "weaned_medication_names"
    ))) %>%
    dplyr::arrange(dplyr::desc(.data$avg_medications_per_month))
}

if (!dir.exists("output/tabs")) {
  dir.create("output/tabs", recursive = TRUE, showWarnings = FALSE)
}

write.csv(med_avgs, "output/tabs/med_avgs.csv", row.names = FALSE)
