# Daily seizure counts in wide format.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(writexl)
})

source("src/desc/seizure_type_standardization.R")

START_DATE <- as.Date("2025-01-01")
END_DATE <- analysis_end_date()
EXCLUDED_PATIENT_IDS <- c(
  "5f9943574483f20039e2ec55",
  "5c09cabbe09148182fa713d3",
  "67eafe28773a9ee7d1911e17",
  "59ae7709696ec5111368197e"
)
STANDARDIZED_GROUPS <- c("Tonic-clonic", "Focal", "Tonic", "Myoclonic", "Absence", "Spasms")

OUTPUT_TAB_DIR <- "output/tabs/seizure_patterns"
OUTPUT_XLSX_PATH <- file.path(OUTPUT_TAB_DIR, "seizures_by_day.xlsx")

dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

parse_event_date <- function(x) {
  parsed <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC"))
  parsed <- dplyr::if_else(is.na(parsed), suppressWarnings(lubridate::ymd(x)), parsed)
  as.Date(parsed)
}

patient_names <- readr::read_csv("data/whatsapp_status.csv", show_col_types = FALSE) %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    first_name = as.character(.data$first_name),
    last_name = as.character(.data$last_name)
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c("first_name", "last_name"),
      ~ dplyr::if_else(
        is.na(.x) | stringr::str_to_lower(stringr::str_squish(.x)) %in% c("", "na", "n/a", "nan", "null"),
        "",
        stringr::str_squish(.x)
      )
    ),
    patient_name = stringr::str_squish(stringr::str_trim(paste(.data$first_name, .data$last_name))),
    patient_name = dplyr::na_if(.data$patient_name, "")
  ) %>%
  dplyr::filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    !is.na(.data$patient_id),
    .data$patient_id != "",
    !.data$patient_id %in% EXCLUDED_PATIENT_IDS
  ) %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(patient_name = dplyr::first(.data$patient_name), .groups = "drop")

seizure_events <- readr::read_csv("data/events.csv", show_col_types = FALSE) %>%
  standardize_seizure_events(filter_to_seizure = TRUE) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    Date = as.Date(.data$event_date)
  ) %>%
  dplyr::filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    !is.na(.data$patient_id),
    .data$patient_id != "",
    !.data$patient_id %in% EXCLUDED_PATIENT_IDS,
    !is.na(.data$Date),
    .data$Date >= START_DATE,
    .data$Date <= END_DATE
  )

app_activity_dates <- jsonlite::fromJSON("data/patient_summary_metrics.json") %>%
  tibble::as_tibble() %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    start_date = parse_event_date(.data$first_app_activity_createdAt),
    start_date = dplyr::if_else(!is.na(.data$start_date) & .data$start_date < START_DATE, START_DATE, .data$start_date)
  ) %>%
  dplyr::filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    is.na(.data$start_date) | .data$start_date <= END_DATE
  )

all_daily_counts <- seizure_events %>%
  dplyr::count(.data$Date, .data$patient_id, name = "seizure_count")

date_grid <- tibble::tibble(Date = seq(START_DATE, END_DATE, by = "day"))

get_patient_order <- function(counts_df) {
  if (nrow(counts_df) == 0) {
    return(character())
  }

  counts_df %>%
    dplyr::group_by(.data$patient_id) %>%
    dplyr::summarise(total_seizures = sum(.data$seizure_count), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$total_seizures), .data$patient_id) %>%
    dplyr::pull(.data$patient_id)
}

all_patient_order <- get_patient_order(all_daily_counts)

if (length(all_patient_order) == 0) {
  all_patient_lookup <- tibble::tibble(patient_id = character(), patient_name = character())
  all_patient_starts <- tibble::tibble(patient_id = character(), start_date = as.Date(character()))
} else {
  all_patient_lookup <- tibble::tibble(patient_id = all_patient_order) %>%
    dplyr::left_join(patient_names, by = "patient_id") %>%
    dplyr::mutate(
      patient_name = dplyr::coalesce(.data$patient_name, .data$patient_id),
      patient_name = make.unique(.data$patient_name, sep = "__")
    )

  all_patient_starts <- tibble::tibble(patient_id = all_patient_order) %>%
    dplyr::left_join(app_activity_dates, by = "patient_id")
}

build_wide_daily_table <- function(counts_df, date_grid_df, patient_order_vec, patient_lookup_df, patient_starts_df) {
  if (length(patient_order_vec) == 0) {
    return(date_grid_df)
  }

  tidyr::expand_grid(
    Date = date_grid_df$Date,
    patient_id = patient_order_vec
  ) %>%
    dplyr::left_join(counts_df, by = c("Date", "patient_id")) %>%
    dplyr::left_join(patient_starts_df, by = "patient_id") %>%
    dplyr::mutate(
      seizure_count = dplyr::case_when(
        is.na(.data$start_date) ~ NA_integer_,
        .data$Date < .data$start_date ~ NA_integer_,
        TRUE ~ tidyr::replace_na(as.integer(.data$seizure_count), 0L)
      )
    ) %>%
    dplyr::left_join(patient_lookup_df, by = "patient_id") %>%
    dplyr::select("Date", "patient_name", "seizure_count") %>%
    tidyr::pivot_wider(
      names_from = "patient_name",
      values_from = "seizure_count",
      values_fill = 0L
    ) %>%
    dplyr::mutate(dplyr::across(-"Date", as.integer)) %>%
    dplyr::arrange(.data$Date)
}

sheet_tables <- list(
  All = build_wide_daily_table(
    all_daily_counts,
    date_grid,
    all_patient_order,
    all_patient_lookup,
    all_patient_starts
  )
)

for (group_name in STANDARDIZED_GROUPS) {
  group_daily_counts <- seizure_events %>%
    dplyr::filter(.data$seizure_group == group_name) %>%
    dplyr::count(.data$Date, .data$patient_id, name = "seizure_count")

  group_patient_order <- get_patient_order(group_daily_counts)
  group_patient_lookup <- all_patient_lookup %>%
    dplyr::filter(.data$patient_id %in% group_patient_order)
  group_patient_starts <- all_patient_starts %>%
    dplyr::filter(.data$patient_id %in% group_patient_order)

  sheet_tables[[group_name]] <- build_wide_daily_table(
    group_daily_counts,
    date_grid,
    group_patient_order,
    group_patient_lookup,
    group_patient_starts
  )
}

writexl::write_xlsx(sheet_tables, OUTPUT_XLSX_PATH)
