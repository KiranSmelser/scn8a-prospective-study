# Shared data correction helpers

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

source("src/analysis_config.R")

if (!exists("standardize_non_rescue_epilepsy_medications")) {
  source("src/desc/medication_standardization.R")
}

MANUAL_CORRECTION_DIR <- "data/manual_corrections"
MEDICATION_ADDITIONS_PATH <- file.path(MANUAL_CORRECTION_DIR, "medication_additions.csv")
MEDICATION_INTERVAL_CORRECTIONS_PATH <- file.path(MANUAL_CORRECTION_DIR, "medication_interval_corrections.csv")

bind_to_template <- function(base_df, additions_df) {
  missing_in_additions <- setdiff(names(base_df), names(additions_df))
  for (column_name in missing_in_additions) {
    additions_df[[column_name]] <- NA
  }

  missing_in_base <- setdiff(names(additions_df), names(base_df))
  for (column_name in missing_in_base) {
    base_df[[column_name]] <- NA
  }

  dplyr::bind_rows(
    base_df,
    additions_df %>% dplyr::select(dplyr::all_of(names(base_df)))
  )
}

read_optional_csv <- function(path) {
  if (!file.exists(path)) {
    return(tibble::tibble())
  }

  readr::read_csv(path, show_col_types = FALSE)
}

read_medications_corrected <- function(
  medications_path = "data/medications.csv",
  additions_path = MEDICATION_ADDITIONS_PATH
) {
  medications <- readr::read_csv(medications_path, show_col_types = FALSE)
  additions <- read_optional_csv(additions_path)

  if (nrow(additions) == 0) {
    return(medications)
  }

  bind_to_template(medications, additions) %>%
    dplyr::distinct(.data$medication_id, .data$patient_id, .keep_all = TRUE)
}

read_medication_interval_corrections <- function(path = MEDICATION_INTERVAL_CORRECTIONS_PATH) {
  corrections <- read_optional_csv(path)

  if (nrow(corrections) == 0) {
    return(tibble::tibble(
      patient_id = character(),
      medication_id = character(),
      name_standardized = character(),
      operation = character(),
      start_date = as.Date(character()),
      end_date = as.Date(character())
    ))
  }

  corrections %>%
    dplyr::transmute(
      patient_id = as.character(.data$patient_id),
      medication_id = as.character(.data$medication_id),
      name_standardized = stringr::str_squish(as.character(.data$name_standardized)),
      operation = stringr::str_to_lower(stringr::str_squish(as.character(.data$operation))),
      start_date = suppressWarnings(lubridate::ymd(.data$start_date)),
      end_date = suppressWarnings(lubridate::ymd(.data$end_date))
    )
}

standardized_medication_lookup <- function(medications_df = read_medications_corrected()) {
  standardize_non_rescue_epilepsy_medications(
    medications_df,
    include_non_drug = FALSE
  ) %>%
    dplyr::mutate(
      patient_id = as.character(.data$patient_id),
      medication_id = as.character(.data$medication_id),
      name_standardized = stringr::str_squish(as.character(.data$name_standardized))
    )
}

read_raw_medication_intervals <- function(
  med_dosages_path = "data/med_dosages.csv",
  med_intakes_path = "data/med_intakes.csv"
) {
  med_dosages <- readr::read_csv(med_dosages_path, show_col_types = FALSE) %>%
    dplyr::transmute(
      patient_id = as.character(.data$patient_id),
      medication_id = as.character(.data$medication_id),
      start_date = suppressWarnings(lubridate::ymd(.data$from)),
      end_date = suppressWarnings(lubridate::ymd(.data$to)),
      interval_source = "dosage"
    )

  med_intakes <- readr::read_csv(med_intakes_path, show_col_types = FALSE) %>%
    dplyr::transmute(
      patient_id = as.character(.data$patient_id),
      medication_id = as.character(.data$medication_id),
      start_date = suppressWarnings(lubridate::ymd(.data$intake_from)),
      end_date = suppressWarnings(lubridate::ymd(.data$intake_to)),
      interval_source = "intake"
    )

  dplyr::bind_rows(med_dosages, med_intakes) %>%
    dplyr::filter(
      !is.na(.data$start_date),
      is.na(.data$end_date) | .data$end_date >= .data$start_date
    ) %>%
    dplyr::distinct()
}

apply_interval_truncations <- function(intervals, corrections) {
  truncations <- corrections %>%
    dplyr::filter(.data$operation == "truncate_at", !is.na(.data$end_date))

  if (nrow(truncations) == 0) {
    return(intervals)
  }

  result <- intervals

  for (i in seq_len(nrow(truncations))) {
    correction <- truncations[i, ]
    target_patient_id <- correction$patient_id
    target_medication_id <- correction$medication_id
    target_name <- correction$name_standardized
    stop_date <- correction$end_date

    result <- result %>%
      dplyr::mutate(
        correction_target = .data$patient_id == target_patient_id &
          (
            (!is.na(target_medication_id) & target_medication_id != "" & .data$medication_id == target_medication_id) |
              (!is.na(target_name) & target_name != "" & .data$name_standardized == target_name)
          ),
        end_date = dplyr::if_else(
          .data$correction_target & (is.na(.data$end_date) | .data$end_date > stop_date),
          stop_date,
          .data$end_date
        )
      ) %>%
      dplyr::filter(!.data$correction_target | .data$start_date <= .data$end_date) %>%
      dplyr::select(-"correction_target")
  }

  result
}

build_added_intervals <- function(corrections, medication_lookup) {
  additions <- corrections %>%
    dplyr::filter(.data$operation == "add_interval", !is.na(.data$start_date)) %>%
    dplyr::select("patient_id", "medication_id", "name_standardized", "start_date", "end_date")

  if (nrow(additions) == 0) {
    return(tibble::tibble(
      patient_id = character(),
      medication_id = character(),
      start_date = as.Date(character()),
      end_date = as.Date(character()),
      interval_source = character(),
      name_standardized = character()
    ))
  }

  additions %>%
    dplyr::left_join(
      medication_lookup %>%
        dplyr::select("patient_id", "medication_id", medication_lookup_name = "name_standardized"),
      by = c("patient_id", "medication_id")
    ) %>%
    dplyr::mutate(
      name_standardized = dplyr::coalesce(.data$name_standardized, .data$medication_lookup_name),
      interval_source = "manual_correction"
    ) %>%
    dplyr::select("patient_id", "medication_id", "start_date", "end_date", "interval_source", "name_standardized")
}

read_medication_intervals_corrected <- function(
  medications_for_analysis = NULL,
  interval_corrections_path = MEDICATION_INTERVAL_CORRECTIONS_PATH,
  med_dosages_path = "data/med_dosages.csv",
  med_intakes_path = "data/med_intakes.csv"
) {
  medication_lookup <- if (is.null(medications_for_analysis)) {
    standardized_medication_lookup()
  } else {
    if (!"name_standardized" %in% names(medications_for_analysis)) {
      medications_for_analysis$name_standardized <- NA_character_
    }
    if (!"name" %in% names(medications_for_analysis)) {
      medications_for_analysis$name <- NA_character_
    }

    medications_for_analysis %>%
      dplyr::mutate(
        patient_id = as.character(.data$patient_id),
        medication_id = as.character(.data$medication_id),
        name_standardized = dplyr::coalesce(
          dplyr::na_if(stringr::str_squish(as.character(.data$name_standardized)), ""),
          dplyr::na_if(stringr::str_squish(as.character(.data$name)), "")
        )
      )
  }

  corrections <- read_medication_interval_corrections(interval_corrections_path)

  intervals <- read_raw_medication_intervals(
    med_dosages_path = med_dosages_path,
    med_intakes_path = med_intakes_path
  ) %>%
    dplyr::semi_join(
      medication_lookup %>% dplyr::select("patient_id", "medication_id"),
      by = c("patient_id", "medication_id")
    ) %>%
    dplyr::left_join(
      medication_lookup %>%
        dplyr::select("patient_id", "medication_id", "name_standardized"),
      by = c("patient_id", "medication_id")
    )

  added_intervals <- build_added_intervals(corrections, medication_lookup)

  dplyr::bind_rows(intervals, added_intervals) %>%
    apply_interval_truncations(corrections) %>%
    dplyr::filter(.data$start_date <= ANALYSIS_CUTOFF_DATE) %>%
    dplyr::mutate(
      end_date = dplyr::if_else(
        !is.na(.data$end_date) & .data$end_date > ANALYSIS_CUTOFF_DATE,
        ANALYSIS_CUTOFF_DATE,
        .data$end_date
      )
    ) %>%
    dplyr::filter(
      !is.na(.data$start_date),
      is.na(.data$end_date) | .data$end_date >= .data$start_date
    ) %>%
    dplyr::distinct()
}
