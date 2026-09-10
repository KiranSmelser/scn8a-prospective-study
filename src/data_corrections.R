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
MEDICATION_EXCLUSIONS_PATH <- file.path(MANUAL_CORRECTION_DIR, "medication_exclusions.csv")
MEDICATION_INTERVAL_CORRECTIONS_PATH <- file.path(MANUAL_CORRECTION_DIR, "medication_interval_corrections.csv")
PROSPECTIVE_SURVEY_SELECTIONS_PATH <- file.path(MANUAL_CORRECTION_DIR, "prospective_survey_selections.csv")
DEVELOPMENT_MILESTONE_ADDITIONS_PATH <- file.path(MANUAL_CORRECTION_DIR, "development_milestone_additions.csv")
SEIZURE_DAILY_COUNT_REPLACEMENTS_PATH <- file.path(MANUAL_CORRECTION_DIR, "seizure_daily_count_replacements.csv")

MICHAEL_SPIEGEL_PATIENT_ID <- "6723f6b40f7cf300409a6197"

# Michael's supplemented history predates his first Helpilepsy activity. Use
# the first corrected seizure as his analysis origin; other patients retain the
# established first-app-activity rule.
EARLIEST_SEIZURE_START_PATIENT_IDS <- c(MICHAEL_SPIEGEL_PATIENT_ID)

select_analysis_start_date <- function(
  patient_id,
  app_activity_date,
  seizure_first_date
) {
  dplyr::if_else(
    as.character(patient_id) %in% EARLIEST_SEIZURE_START_PATIENT_IDS,
    as.Date(seizure_first_date),
    as.Date(app_activity_date)
  )
}

# Michael's manually supplemented medication metadata did not carry the
# Helpilepsy "reason" value used by the analysis medication filter. Keep this
# correction limited to medications that are epilepsy treatments; rescue and
# non-epilepsy/support medications retain their source classifications.
MICHAEL_SPIEGEL_EPILEPSY_MEDICATION_IDS <- c(
  "0f8110ddb75d0b297e4016fc", # phenobarbital
  "11872d30a33b5e22b0d0baff", # brivaracetam
  "18ae973f0bb087d7a0b0cd35", # cenobamate
  "191206958b85f538da0666d4", # lacosamide
  "2414642c156e48d5525134f9", # levetiracetam
  "5282b81a1f6548b6377f16b4", # sodium valproate
  "9d9e60b3161d8460b898a39c", # carbamazepine
  "a4b4b562097134a32dc9805b", # clobazam
  "d78c333e768ea6367a33e01c", # ACTH
  "f101e7079e764e6b24134ab2", # topiramate
  "fcf3c000eb1aa307f7470534"  # clonazepam
)

# The manual history retained the original Helpilepsy medication IDs in
# med_dosages.csv, while medications.csv uses deterministic import IDs. Map the
# former to the latter so the authoritative dosage intervals join to metadata.
MICHAEL_SPIEGEL_MEDICATION_ID_REPLACEMENTS <- tibble::tribble(
  ~source_medication_id, ~medication_id,
  "67249da30f7cf300409a658e", "5282b81a1f6548b6377f16b4", # Depalept
  "67249de00f7cf300409a6594", "f101e7079e764e6b24134ab2", # Topamax
  "67249e300f7cf300409a659e", "fcf3c000eb1aa307f7470534", # Rivotril
  "67249e6e0f7cf300409a65a9", "a4b4b562097134a32dc9805b", # Frisium
  "675c234500855b004cfebe62", "18ae973f0bb087d7a0b0cd35"  # Xcopri
)

parse_corrected_datetime <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt")) {
    return(as.POSIXct(x, tz = tz))
  }
  if (inherits(x, "Date")) {
    return(as.POSIXct(x, tz = tz))
  }

  value <- as.character(x)
  parsed <- suppressWarnings(lubridate::ymd_hms(value, quiet = TRUE, tz = tz))
  parsed <- dplyr::coalesce(
    parsed,
    suppressWarnings(lubridate::ymd(value, quiet = TRUE, tz = tz)),
    suppressWarnings(lubridate::mdy_hms(value, quiet = TRUE, tz = tz)),
    suppressWarnings(lubridate::mdy(value, quiet = TRUE, tz = tz))
  )
  parsed
}

parse_corrected_date <- function(x) {
  as.Date(parse_corrected_datetime(x))
}

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

read_seizure_daily_count_replacements <- function(
  path = SEIZURE_DAILY_COUNT_REPLACEMENTS_PATH
) {
  replacements <- read_optional_csv(path)

  if (nrow(replacements) == 0) {
    return(tibble::tibble(
      patient_id = character(),
      event_date = as.Date(character()),
      seizure_type = character(),
      count = integer(),
      source_seizure_type = character()
    ))
  }

  required_columns <- c(
    "patient_id", "event_date", "seizure_type", "count", "source_seizure_type"
  )
  missing_columns <- setdiff(required_columns, names(replacements))
  if (length(missing_columns) > 0) {
    stop(
      "Seizure daily count replacements are missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  replacements <- replacements %>%
    dplyr::transmute(
      patient_id = stringr::str_squish(as.character(.data$patient_id)),
      event_date = parse_corrected_date(.data$event_date),
      seizure_type = stringr::str_squish(as.character(.data$seizure_type)),
      count = suppressWarnings(as.integer(.data$count)),
      source_seizure_type = stringr::str_squish(as.character(.data$source_seizure_type))
    )

  invalid <- replacements %>%
    dplyr::filter(
      is.na(.data$patient_id) | .data$patient_id == "" |
        is.na(.data$event_date) |
        is.na(.data$seizure_type) | .data$seizure_type == "" |
        is.na(.data$count) | .data$count < 0 |
        is.na(.data$source_seizure_type) | .data$source_seizure_type == ""
    )
  if (nrow(invalid) > 0) {
    stop("Seizure daily count replacements contain invalid values.", call. = FALSE)
  }

  duplicate_keys <- replacements %>%
    dplyr::count(.data$patient_id, .data$event_date, .data$seizure_type) %>%
    dplyr::filter(.data$n > 1)
  if (nrow(duplicate_keys) > 0) {
    stop("Seizure daily count replacements contain duplicate patient/date/type rows.", call. = FALSE)
  }

  replacements
}

read_events_corrected <- function(
  events_path = "data/events.csv",
  replacements_path = SEIZURE_DAILY_COUNT_REPLACEMENTS_PATH
) {
  events <- readr::read_csv(events_path, show_col_types = FALSE) %>%
    dplyr::mutate(
      event_id = as.character(.data$event_id),
      patient_id = as.character(.data$patient_id)
    )
  replacements <- read_seizure_daily_count_replacements(replacements_path)

  if (nrow(replacements) == 0) {
    return(events)
  }

  replacement_keys <- replacements %>%
    dplyr::transmute(
      patient_id = .data$patient_id,
      event_date = .data$event_date,
      source_seizure_type_key = stringr::str_to_lower(.data$source_seizure_type)
    ) %>%
    dplyr::distinct()

  retained_events <- events %>%
    dplyr::mutate(
      event_date = parse_corrected_date(.data$date),
      source_seizure_type_key = stringr::str_to_lower(
        stringr::str_squish(as.character(.data$seizure_type))
      ),
      event_type_key = stringr::str_to_lower(stringr::str_squish(as.character(.data$type)))
    ) %>%
    dplyr::left_join(
      replacement_keys %>% dplyr::mutate(replace_source_event = TRUE),
      by = c("patient_id", "event_date", "source_seizure_type_key")
    ) %>%
    dplyr::filter(!(.data$event_type_key == "seizure" & dplyr::coalesce(.data$replace_source_event, FALSE))) %>%
    dplyr::select(-"event_date", -"source_seizure_type_key", -"event_type_key", -"replace_source_event")

  added_events <- replacements %>%
    dplyr::filter(.data$count > 0) %>%
    tidyr::uncount(.data$count, .id = "event_number") %>%
    dplyr::transmute(
      event_id = paste0(
        "manual_", .data$patient_id, "_",
        format(.data$event_date, "%Y%m%d"), "_",
        stringr::str_pad(.data$event_number, width = 3, pad = "0")
      ),
      patient_id = .data$patient_id,
      type = "seizure",
      date = parse_corrected_datetime(
        paste0(format(.data$event_date, "%Y-%m-%d"), "T12:00:00.000Z")
      ),
      seizure_type = .data$seizure_type,
      createdAt = parse_corrected_datetime(
        paste0(format(.data$event_date, "%Y-%m-%d"), "T12:00:00.000Z")
      ),
      updatedAt = parse_corrected_datetime(
        paste0(format(.data$event_date, "%Y-%m-%d"), "T12:00:00.000Z")
      )
    )

  bind_to_template(retained_events, added_events) %>%
    dplyr::distinct(.data$event_id, .keep_all = TRUE)
}

read_prospective_survey_selections <- function(path = PROSPECTIVE_SURVEY_SELECTIONS_PATH) {
  selections <- read_optional_csv(path)

  if (nrow(selections) == 0) {
    return(tibble::tibble(
      patient_id = character(),
      participant_id = character()
    ))
  }

  required_columns <- c("patient_id", "participant_id")
  missing_columns <- setdiff(required_columns, names(selections))
  if (length(missing_columns) > 0) {
    stop(
      "Prospective survey selections are missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  selections <- selections %>%
    dplyr::transmute(
      patient_id = stringr::str_squish(as.character(.data$patient_id)),
      participant_id = stringr::str_squish(as.character(.data$participant_id))
    ) %>%
    dplyr::filter(
      !is.na(.data$patient_id),
      !is.na(.data$participant_id),
      .data$patient_id != "",
      .data$participant_id != ""
    )

  duplicate_patients <- selections %>%
    dplyr::count(.data$patient_id) %>%
    dplyr::filter(.data$n > 1)
  if (nrow(duplicate_patients) > 0) {
    stop(
      "Prospective survey selections contain multiple entries for patient(s): ",
      paste(duplicate_patients$patient_id, collapse = ", "),
      call. = FALSE
    )
  }

  selections
}

apply_prospective_survey_selections <- function(
  surveys,
  selections = read_prospective_survey_selections()
) {
  required_columns <- c("patient_id", "participant_id")
  missing_columns <- setdiff(required_columns, names(surveys))
  if (length(missing_columns) > 0) {
    stop(
      "Prospective surveys are missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  surveys <- surveys %>%
    dplyr::mutate(
      patient_id = as.character(.data$patient_id),
      participant_id = as.character(.data$participant_id)
    )

  if (nrow(selections) == 0) {
    return(surveys)
  }

  missing_selected_surveys <- selections %>%
    dplyr::anti_join(
      surveys %>% dplyr::distinct(.data$patient_id, .data$participant_id),
      by = c("patient_id", "participant_id")
    )
  if (nrow(missing_selected_surveys) > 0) {
    missing_labels <- paste0(
      missing_selected_surveys$patient_id,
      " -> ",
      missing_selected_surveys$participant_id
    )
    stop(
      "Selected prospective survey(s) were not found: ",
      paste(missing_labels, collapse = ", "),
      call. = FALSE
    )
  }

  surveys %>%
    dplyr::left_join(
      selections %>% dplyr::rename(selected_participant_id = "participant_id"),
      by = "patient_id"
    ) %>%
    dplyr::filter(
      is.na(.data$selected_participant_id) |
        .data$participant_id == .data$selected_participant_id
    ) %>%
    dplyr::select(-"selected_participant_id")
}

read_prospective_surveys_corrected <- function(
  surveys_path = "data/prospective_surveys.csv",
  selections_path = PROSPECTIVE_SURVEY_SELECTIONS_PATH
) {
  surveys <- readr::read_csv(surveys_path, show_col_types = FALSE)
  selections <- read_prospective_survey_selections(selections_path)
  apply_prospective_survey_selections(surveys, selections)
}

read_prospective_child_records_corrected <- function(path, surveys) {
  if (!"survey_instance_id" %in% names(surveys)) {
    stop("Corrected prospective surveys are missing survey_instance_id.", call. = FALSE)
  }

  records <- readr::read_csv(path, show_col_types = FALSE)
  if (!"survey_instance_id" %in% names(records)) {
    stop(path, " is missing survey_instance_id.", call. = FALSE)
  }

  records %>%
    dplyr::mutate(survey_instance_id = as.character(.data$survey_instance_id)) %>%
    dplyr::semi_join(
      surveys %>%
        dplyr::transmute(survey_instance_id = as.character(.data$survey_instance_id)) %>%
        dplyr::distinct(),
      by = "survey_instance_id"
    )
}

read_development_milestones_corrected <- function(
  surveys = read_prospective_surveys_corrected(),
  milestones_path = "data/prospective_development_milestones.csv",
  additions_path = DEVELOPMENT_MILESTONE_ADDITIONS_PATH
) {
  survey_lookup <- surveys %>%
    dplyr::transmute(
      survey_instance_id = as.character(.data$survey_instance_id),
      patient_id = as.character(.data$patient_id),
      event_date = as.Date(suppressWarnings(lubridate::ymd_hms(
        .data$prospective_study_timestamp_utc,
        quiet = TRUE,
        tz = "UTC"
      ))),
      event_date = dplyr::coalesce(
        .data$event_date,
        as.Date(suppressWarnings(lubridate::mdy_hm(
          .data$prospective_study_timestamp,
          quiet = TRUE,
          tz = "UTC"
        )))
      )
    ) %>%
    dplyr::distinct(.data$survey_instance_id, .keep_all = TRUE)

  prospective_records <- read_prospective_child_records_corrected(
    milestones_path,
    surveys
  ) %>%
    dplyr::mutate(survey_instance_id = as.character(.data$survey_instance_id)) %>%
    dplyr::inner_join(survey_lookup, by = "survey_instance_id") %>%
    dplyr::transmute(
      survey_instance_id = .data$survey_instance_id,
      patient_id = .data$patient_id,
      milestone = as.character(.data$milestone),
      status = suppressWarnings(as.numeric(.data$status)),
      event_date = .data$event_date,
      data_source = "Prospective survey"
    )

  additions <- read_optional_csv(additions_path)
  if (nrow(additions) == 0) {
    return(prospective_records)
  }

  required_columns <- c("patient_id", "milestone", "status", "event_date")
  missing_columns <- setdiff(required_columns, names(additions))
  if (length(missing_columns) > 0) {
    stop(
      "Development milestone additions are missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  manual_records <- additions %>%
    dplyr::transmute(
      survey_instance_id = NA_character_,
      patient_id = stringr::str_squish(as.character(.data$patient_id)),
      milestone = stringr::str_squish(as.character(.data$milestone)),
      status = suppressWarnings(as.numeric(.data$status)),
      event_date = suppressWarnings(lubridate::ymd(.data$event_date)),
      data_source = "Manual correction"
    ) %>%
    dplyr::filter(
      !is.na(.data$patient_id),
      .data$patient_id != "",
      !is.na(.data$milestone),
      .data$milestone != "",
      !is.na(.data$status),
      !is.na(.data$event_date)
    )

  dplyr::bind_rows(prospective_records, manual_records) %>%
    dplyr::distinct(
      .data$patient_id,
      .data$milestone,
      .data$status,
      .data$event_date,
      .data$data_source,
      .keep_all = TRUE
    )
}

read_medication_exclusions <- function(path = MEDICATION_EXCLUSIONS_PATH) {
  exclusions <- read_optional_csv(path)

  if (nrow(exclusions) == 0) {
    return(tibble::tibble(
      patient_id = character(),
      medication_id = character()
    ))
  }

  exclusions %>%
    dplyr::transmute(
      patient_id = as.character(.data$patient_id),
      medication_id = as.character(.data$medication_id)
    ) %>%
    dplyr::filter(
      !is.na(.data$patient_id),
      !is.na(.data$medication_id),
      .data$patient_id != "",
      .data$medication_id != ""
    ) %>%
    dplyr::distinct()
}

read_medications_corrected <- function(
  medications_path = "data/medications.csv",
  additions_path = MEDICATION_ADDITIONS_PATH,
  exclusions_path = MEDICATION_EXCLUSIONS_PATH
) {
  medications <- readr::read_csv(medications_path, show_col_types = FALSE) %>%
    dplyr::mutate(
      patient_id = as.character(.data$patient_id),
      medication_id = as.character(.data$medication_id),
      reason = dplyr::if_else(
        .data$patient_id == MICHAEL_SPIEGEL_PATIENT_ID &
          .data$medication_id %in% MICHAEL_SPIEGEL_EPILEPSY_MEDICATION_IDS,
        "epilepsy",
        as.character(.data$reason)
      )
    )
  exclusions <- read_medication_exclusions(exclusions_path)

  if (nrow(exclusions) > 0) {
    medications <- medications %>%
      dplyr::anti_join(exclusions, by = c("patient_id", "medication_id"))
  }

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
      source_medication_id = as.character(.data$medication_id),
      start_date = parse_corrected_date(.data$from),
      end_date = parse_corrected_date(.data$to),
      interval_source = "dosage"
    ) %>%
    dplyr::left_join(
      MICHAEL_SPIEGEL_MEDICATION_ID_REPLACEMENTS,
      by = "source_medication_id"
    ) %>%
    dplyr::mutate(
      medication_id = dplyr::if_else(
        .data$patient_id == MICHAEL_SPIEGEL_PATIENT_ID & !is.na(.data$medication_id),
        .data$medication_id,
        .data$source_medication_id
      )
    ) %>%
    dplyr::select("patient_id", "medication_id", "start_date", "end_date", "interval_source")

  medication_ids_with_authoritative_dosages <- med_dosages %>%
    dplyr::filter(.data$patient_id == MICHAEL_SPIEGEL_PATIENT_ID) %>%
    dplyr::distinct(.data$medication_id)

  med_intakes_raw <- readr::read_csv(med_intakes_path, show_col_types = FALSE) %>%
    dplyr::mutate(
      patient_id = as.character(.data$patient_id),
      medication_id = as.character(.data$medication_id)
    )

  med_intakes <- med_intakes_raw %>%
    dplyr::filter(.data$patient_id != MICHAEL_SPIEGEL_PATIENT_ID) %>%
    dplyr::transmute(
      patient_id = .data$patient_id,
      medication_id = .data$medication_id,
      start_date = parse_corrected_date(.data$intake_from),
      end_date = parse_corrected_date(.data$intake_to),
      interval_source = "intake"
    )

  # Michael's manually supplemented med_intakes rows are daily history records,
  # not open-ended schedule definitions. Collapse them to bounded spans and use
  # dosage intervals instead whenever an authoritative dosage series is present.
  michael_manual_history <- med_intakes_raw %>%
    dplyr::filter(
      .data$patient_id == MICHAEL_SPIEGEL_PATIENT_ID,
      !.data$medication_id %in% medication_ids_with_authoritative_dosages$medication_id
    ) %>%
    dplyr::mutate(
      history_date = parse_corrected_date(.data$taken_date),
      taken_flag = stringr::str_to_lower(stringr::str_squish(as.character(.data$taken))) %in%
        c("true", "1", "yes", "y")
    ) %>%
    dplyr::filter(.data$taken_flag, !is.na(.data$history_date)) %>%
    dplyr::group_by(.data$patient_id, .data$medication_id) %>%
    dplyr::summarise(
      start_date = min(.data$history_date),
      end_date = max(.data$history_date),
      interval_source = "manual_history",
      .groups = "drop"
    )

  dplyr::bind_rows(med_dosages, med_intakes, michael_manual_history) %>%
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

apply_interval_replacements <- function(intervals, corrections) {
  replacement_targets <- corrections %>%
    dplyr::filter(.data$operation == "replace_interval") %>%
    dplyr::distinct(.data$patient_id, .data$medication_id, .data$name_standardized)

  if (nrow(replacement_targets) == 0) {
    return(intervals)
  }

  result <- intervals

  for (i in seq_len(nrow(replacement_targets))) {
    target <- replacement_targets[i, ]
    target_patient_id <- target$patient_id
    target_medication_id <- target$medication_id
    target_name <- target$name_standardized

    result <- result %>%
      dplyr::filter(
        !(
          .data$patient_id == target_patient_id &
            (
              (!is.na(target_medication_id) & target_medication_id != "" & .data$medication_id == target_medication_id) |
                (!is.na(target_name) & target_name != "" & .data$name_standardized == target_name)
            )
        )
      )
  }

  result
}

build_added_intervals <- function(corrections, medication_lookup) {
  additions <- corrections %>%
    dplyr::filter(
      .data$operation %in% c("add_interval", "replace_interval"),
      !is.na(.data$start_date)
    ) %>%
    dplyr::select("patient_id", "medication_id", "name_standardized", "operation", "start_date", "end_date")

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
      interval_source = dplyr::if_else(
        .data$operation == "replace_interval",
        "manual_replacement",
        "manual_correction"
      )
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
    ) %>%
    apply_interval_replacements(corrections)

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
