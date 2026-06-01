# Module overlap tables

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

source("src/analysis_config.R")

EVENTS_PATH <- "data/events.csv"
MED_INTAKES_PATH <- "data/med_intakes.csv"
MOOD_SLEEP_PATH <- "data/mood_sleep.csv"
PROSPECTIVE_SURVEYS_PATH <- "data/prospective_surveys.csv"
FORMS_PATH <- "data/forms.csv"

OUTPUT_TAB_DIR <- "output/tabs/module_overlap"
OUTPUT_INTERSECTIONS_PATH <- file.path(OUTPUT_TAB_DIR, "module_overlap_intersections.csv")
OUTPUT_COUNTS_PATH <- file.path(OUTPUT_TAB_DIR, "module_overlap_counts.csv")
OUTPUT_PAIRWISE_PATH <- file.path(OUTPUT_TAB_DIR, "module_overlap_pairwise.csv")
OUTPUT_MEMBERSHIP_PATH <- file.path(OUTPUT_TAB_DIR, "module_overlap_membership.csv")

normalize_patient_ids <- function(x) {
  x %>%
    as.character() %>%
    trimws() %>%
    .[!is.na(.) & nzchar(.)] %>%
    unique() %>%
    sort()
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    return(tibble::tibble())
  }
  readr::read_csv(path, show_col_types = FALSE)
}

filter_records_on_or_before_cutoff <- function(data, date_column) {
  if (!date_column %in% names(data)) {
    return(data[0, , drop = FALSE])
  }

  data %>%
    dplyr::mutate(
      analysis_record_date = as.Date(suppressWarnings(lubridate::parse_date_time(
        as.character(.data[[date_column]]),
        orders = c("ymd HMS", "ymd HM", "ymd", "mdy HMS", "mdy HM", "mdy"),
        tz = "UTC",
        quiet = TRUE
      )))
    ) %>%
    dplyr::filter(!is.na(.data$analysis_record_date), .data$analysis_record_date <= ANALYSIS_CUTOFF_DATE) %>%
    dplyr::select(-"analysis_record_date")
}

get_seizure_event_patients <- function() {
  events <- safe_read_csv(EVENTS_PATH)
  if (!all(c("patient_id", "type") %in% names(events))) {
    return(character())
  }
  events %>%
    filter_records_on_or_before_cutoff("date") %>%
    dplyr::mutate(type_clean = stringr::str_to_lower(stringr::str_squish(as.character(.data$type)))) %>%
    dplyr::filter(.data$type_clean == "seizure") %>%
    dplyr::pull(.data$patient_id) %>%
    normalize_patient_ids()
}

get_med_intake_patients <- function() {
  med_intakes <- safe_read_csv(MED_INTAKES_PATH)
  if (!"patient_id" %in% names(med_intakes)) {
    return(character())
  }
  med_intakes %>%
    filter_records_on_or_before_cutoff("date") %>%
    dplyr::pull(.data$patient_id) %>%
    normalize_patient_ids()
}

get_mood_sleep_patients <- function() {
  mood_sleep <- safe_read_csv(MOOD_SLEEP_PATH)
  if (!"patient_id" %in% names(mood_sleep)) {
    return(character())
  }
  mood_sleep %>%
    filter_records_on_or_before_cutoff("date") %>%
    dplyr::pull(.data$patient_id) %>%
    normalize_patient_ids()
}

get_prospective_matched_patients <- function() {
  surveys <- safe_read_csv(PROSPECTIVE_SURVEYS_PATH)
  if (!"patient_id" %in% names(surveys)) {
    return(character())
  }
  surveys %>%
    filter_records_on_or_before_cutoff("prospective_study_timestamp_utc") %>%
    dplyr::pull(.data$patient_id) %>%
    normalize_patient_ids()
}

get_weekly_diary_patients <- function() {
  forms <- safe_read_csv(FORMS_PATH)
  if (!all(c("patient_id", "form_name") %in% names(forms))) {
    return(character())
  }
  forms %>%
    filter_records_on_or_before_cutoff("date") %>%
    dplyr::mutate(form_name_clean = stringr::str_to_lower(stringr::str_squish(as.character(.data$form_name)))) %>%
    dplyr::filter(.data$form_name_clean == "scn8a diary completion") %>%
    dplyr::pull(.data$patient_id) %>%
    normalize_patient_ids()
}

module_labels <- c(
  events = "Seizure Events",
  med_intakes = "Medication Intakes",
  mood_sleep = "Mood/Sleep",
  prospective_matched = "Prospective Survey (Matched)",
  weekly_diary = "Weekly SCN8A Diary"
)

module_members <- list(
  events = get_seizure_event_patients(),
  med_intakes = get_med_intake_patients(),
  mood_sleep = get_mood_sleep_patients(),
  prospective_matched = get_prospective_matched_patients(),
  weekly_diary = get_weekly_diary_patients()
)

module_keys <- names(module_members)
all_patients <- sort(unique(unlist(module_members, use.names = FALSE)))

if (length(all_patients) == 0) {
  warning("No patient IDs found across modules; outputs were not generated.")
} else {
  presence <- tibble::tibble(patient_id = all_patients)
  for (key in module_keys) {
    presence[[key]] <- as.integer(presence$patient_id %in% module_members[[key]])
  }

  module_counts <- tibble::tibble(
    module_key = module_keys,
    module_label = unname(module_labels[module_keys]),
    n_patients = purrr::map_int(module_keys, ~ sum(presence[[.x]] == 1L))
  ) %>%
    dplyr::arrange(dplyr::desc(.data$n_patients), .data$module_label)

  intersections <- presence %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(module_keys))) %>%
    dplyr::summarise(
      n_patients = dplyr::n(),
      patient_ids = paste(sort(patient_id), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_modules = rowSums(dplyr::across(dplyr::all_of(module_keys))),
      modules_present = purrr::pmap_chr(
        dplyr::across(dplyr::all_of(module_keys)),
        function(...) {
          flags <- as.logical(c(...))
          labels <- unname(module_labels[module_keys])[flags]
          if (length(labels) == 0) "None" else paste(labels, collapse = " | ")
        }
      )
    ) %>%
    dplyr::arrange(dplyr::desc(.data$n_patients), dplyr::desc(.data$n_modules), .data$modules_present) %>%
    dplyr::mutate(intersection_id = paste0("I", dplyr::row_number())) %>%
    dplyr::select(
      intersection_id,
      n_patients,
      n_modules,
      modules_present,
      dplyr::all_of(module_keys),
      patient_ids
    )

  pairwise <- tidyr::expand_grid(module_a = module_keys, module_b = module_keys) %>%
    dplyr::mutate(
      module_a_label = unname(module_labels[.data$module_a]),
      module_b_label = unname(module_labels[.data$module_b]),
      n_overlap = purrr::map2_int(
        .data$module_a,
        .data$module_b,
        ~ sum(presence[[.x]] == 1L & presence[[.y]] == 1L)
      )
    ) %>%
    dplyr::select(module_a, module_a_label, module_b, module_b_label, n_overlap)

  if (!dir.exists(OUTPUT_TAB_DIR)) {
    dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
  }

  readr::write_csv(intersections, OUTPUT_INTERSECTIONS_PATH)
  readr::write_csv(module_counts, OUTPUT_COUNTS_PATH)
  readr::write_csv(pairwise, OUTPUT_PAIRWISE_PATH)
  readr::write_csv(presence, OUTPUT_MEMBERSHIP_PATH)
}
