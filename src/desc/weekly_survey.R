# Weekly diary completion summary table

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

FORMS_PATH <- "data/forms.csv"
ANSWERS_PATH <- "data/form_answers.csv"
OUTPUT_DIR <- "output/tabs"
OUTPUT_PATH <- file.path(OUTPUT_DIR, "weekly_survey.csv")

QUESTION_LABELS <- c(
  "q1_1_did_you_give_your_child_all_medications_every_day_this_week" =
    "medications_administered_all_week",
  "q2_2_does_the_daily_seizure_record_in_the_app_this_week_accurately_reflect_your_child_s_seizure_count_including_if_there_were_none" =
    "seizure_record_matches",
  "q3_3_on_a_scale_of_1_10_how_was_your_child_s_mood_this_week" =
    "child_mood_rating",
  "q4_4_on_a_scale_of_1_10_how_was_your_child_s_sleep_quality_this_week" =
    "sleep_quality_rating"
)

NUMERIC_COLUMNS <- c("child_mood_rating", "sleep_quality_rating")
EXPECTED_ANSWER_COLUMNS <- c(
  "medications_administered_all_week",
  "seizure_record_matches",
  "child_mood_rating",
  "sleep_quality_rating"
)

scn8a_forms <- readr::read_csv(FORMS_PATH, show_col_types = FALSE) %>%
  dplyr::mutate(
    form_name_clean = stringr::str_to_lower(stringr::str_squish(.data$form_name)),
    completion_timestamp_utc = lubridate::ymd_hms(.data$date, quiet = TRUE, tz = "UTC"),
    completion_date = as.Date(.data$completion_timestamp_utc)
  ) %>%
  dplyr::filter(!is.na(.data$form_name_clean), .data$form_name_clean == "scn8a diary completion") %>%
  dplyr::mutate(
    week_start = dplyr::if_else(
      !is.na(.data$completion_date),
      lubridate::floor_date(.data$completion_date, unit = "week", week_start = 1),
      as.Date(NA)
    ),
    week_end = dplyr::if_else(
      !is.na(.data$week_start),
      .data$week_start + lubridate::days(6),
      as.Date(NA)
    )
  ) %>%
  dplyr::transmute(
    form_id,
    patient_id,
    form_label = label,
    completion_timestamp_utc,
    completion_date,
    week_start,
    week_end
  )

if (nrow(scn8a_forms) == 0) {
  weekly_summary <- tibble::tibble(
    form_id = character(),
    patient_id = character(),
    form_label = character(),
    completion_timestamp_utc = as.POSIXct(character()),
    completion_date = as.Date(character()),
    week_start = as.Date(character()),
    week_end = as.Date(character()),
    medications_administered_all_week = character(),
    seizure_record_matches = character(),
    child_mood_rating = numeric(),
    sleep_quality_rating = numeric()
  )
} else {
  answers_filtered <- readr::read_csv(ANSWERS_PATH, show_col_types = FALSE) %>%
    dplyr::semi_join(scn8a_forms %>% dplyr::select(form_id), by = "form_id") %>%
    dplyr::mutate(
      question_key = dplyr::recode(
        .data$question_code,
        !!!QUESTION_LABELS,
        .default = NA_character_,
        .missing = NA_character_
      ),
      answer_value = stringr::str_squish(as.character(.data$answer)),
      answer_value = dplyr::na_if(.data$answer_value, "")
    ) %>%
    dplyr::filter(!is.na(.data$question_key))

  answers_wide <- if (nrow(answers_filtered) == 0) {
    tibble::tibble(form_id = character())
  } else {
    answers_filtered %>%
      dplyr::arrange(.data$form_id, .data$question_index, .data$question_key) %>%
      dplyr::group_by(.data$form_id, .data$question_key) %>%
      dplyr::summarise(
        answer_value = dplyr::coalesce(dplyr::first(.data$answer_value), NA_character_),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(
        id_cols = form_id,
        names_from = question_key,
        values_from = answer_value
      )
  }

  missing_cols <- setdiff(EXPECTED_ANSWER_COLUMNS, names(answers_wide))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      answers_wide[[col]] <- NA_character_
    }
  }

  answers_wide <- answers_wide %>%
    dplyr::mutate(
      dplyr::across(
        tidyselect::any_of(NUMERIC_COLUMNS),
        ~ suppressWarnings(as.numeric(.x))
      )
    )

  weekly_summary <- scn8a_forms %>%
    dplyr::left_join(answers_wide, by = "form_id") %>%
    dplyr::arrange(
      dplyr::desc(.data$completion_timestamp_utc),
      .data$patient_id
    )

  weekly_summary <- weekly_summary %>%
    dplyr::select(
      form_id,
      patient_id,
      form_label,
      completion_timestamp_utc,
      completion_date,
      week_start,
      week_end,
      tidyselect::all_of(EXPECTED_ANSWER_COLUMNS)
    )
}

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

write.csv(weekly_summary, OUTPUT_PATH, row.names = FALSE)
