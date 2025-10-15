suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(jsonlite)
})

analysis_end <- Sys.Date()

parse_bool <- function(x) {
  stringr::str_to_lower(trimws(as.character(x))) %in% c("true", "1", "yes", "y")
}

normalise_date <- function(x) {
  dt <- suppressWarnings(lubridate::ymd_hms(x, quiet = TRUE, tz = "UTC"))
  as.Date(dt)
}

usage_parts <- list()

if (file.exists("data/events.csv")) {
  events <- readr::read_csv("data/events.csv", show_col_types = FALSE) %>%
    dplyr::mutate(date_only = normalise_date(date)) %>%
    dplyr::filter(!is.na(patient_id), !is.na(date_only)) %>%
    dplyr::filter(date_only <= analysis_end) %>%
    dplyr::select(patient_id, date_only) %>%
    dplyr::distinct()

  if (nrow(events) > 0) {
    usage_parts <- append(usage_parts, list(events))
  }
}

if (file.exists("data/forms.csv")) {
  forms <- readr::read_csv("data/forms.csv", show_col_types = FALSE) %>%
    dplyr::mutate(date_only = normalise_date(date)) %>%
    dplyr::filter(!is.na(patient_id), !is.na(date_only)) %>%
    dplyr::filter(date_only <= analysis_end) %>%
    dplyr::select(patient_id, date_only) %>%
    dplyr::distinct()

  if (nrow(forms) > 0) {
    usage_parts <- append(usage_parts, list(forms))
  }
}

if (file.exists("data/med_intakes.csv")) {
  med_intakes <- readr::read_csv("data/med_intakes.csv", show_col_types = FALSE) %>%
    dplyr::mutate(
      taken_flag = parse_bool(taken),
      date_only = normalise_date(taken_date)
    ) %>%
    dplyr::filter(taken_flag, !is.na(patient_id), !is.na(date_only)) %>%
    dplyr::filter(date_only <= analysis_end) %>%
    dplyr::select(patient_id, date_only) %>%
    dplyr::distinct()

  if (nrow(med_intakes) > 0) {
    usage_parts <- append(usage_parts, list(med_intakes))
  }
}

if (length(usage_parts) == 0) {
  app_usage_summary <- tibble::tibble(
    patient_id = character(),
    app_usage_days_total = integer(),
    app_usage_days_latest_month = integer(),
    app_usage_per_month = vector("list", 0)
  )
} else {
  usage_all <- dplyr::bind_rows(usage_parts)

  usage_days <- usage_all %>%
    dplyr::distinct(patient_id, date_only)

  usage_totals <- usage_days %>%
    dplyr::count(patient_id, name = "app_usage_days_total")

  usage_monthly <- usage_days %>%
    dplyr::mutate(month_date = lubridate::floor_date(date_only, "month")) %>%
    dplyr::count(patient_id, month_date, name = "usage_days") %>%
    dplyr::arrange(patient_id, month_date)

  latest_month <- usage_monthly %>%
    dplyr::group_by(patient_id) %>%
    dplyr::slice_max(order_by = month_date, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      patient_id,
      app_usage_days_latest_month = as.integer(usage_days)
    )

  usage_per_month <- usage_monthly %>%
    dplyr::mutate(
      month_label = format(month_date, "%Y-%m"),
      usage_days = as.integer(usage_days)
    ) %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(
      app_usage_per_month = list(as.list(stats::setNames(usage_days, month_label))),
      .groups = "drop"
    )

  unique_patients <- usage_days %>%
    dplyr::distinct(patient_id)

  app_usage_summary <- unique_patients %>%
    dplyr::left_join(usage_totals, by = "patient_id") %>%
    dplyr::left_join(latest_month, by = "patient_id") %>%
    dplyr::left_join(usage_per_month, by = "patient_id") %>%
    dplyr::mutate(
      app_usage_days_total = dplyr::coalesce(app_usage_days_total, 0L),
      app_usage_days_latest_month = dplyr::coalesce(app_usage_days_latest_month, 0L),
      app_usage_per_month = purrr::map(
        app_usage_per_month,
        function(x) {
          if (is.null(x) || (length(x) == 1 && all(is.na(x)))) {
            list()
          } else {
            x
          }
        }
      )
    ) %>%
    dplyr::arrange(patient_id)
}

if (!dir.exists("output/tabs")) {
  dir.create("output/tabs", recursive = TRUE, showWarnings = FALSE)
}

records <- purrr::map(seq_len(nrow(app_usage_summary)), function(i) {
  list(
    patient_id = app_usage_summary$patient_id[[i]],
    app_usage_days_total = as.integer(app_usage_summary$app_usage_days_total[[i]]),
    app_usage_days_latest_month = as.integer(app_usage_summary$app_usage_days_latest_month[[i]]),
    app_usage_per_month = app_usage_summary$app_usage_per_month[[i]]
  )
})

jsonlite::write_json(
  records,
  "output/tabs/app_usage.json",
  pretty = TRUE,
  auto_unbox = TRUE,
  null = "null"
)