# Build patient-month seizure panel for the compliant cohort of patients.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

source("src/desc/medication_standardization.R")
source("src/data_corrections.R")

OUTPUT_TAB_DIR <- "output/tabs/modeling"
OUTPUT_PANEL_PATH <- file.path(OUTPUT_TAB_DIR, "patient_month_panel.csv")
MED_EXPOSURE_MIN_PATIENTS <- 4L
MED_EXPOSURE_MIN_PATIENT_MONTHS <- 24L
MED_PAIR_EXPOSURE_MIN_PATIENTS <- 4L
MED_PAIR_EXPOSURE_MIN_PATIENT_MONTHS <- 20L
MED_PAIR_EXPOSURE_MIN_SINGLE_ONLY_MONTHS <- 15L

dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

current_date <- analysis_end_date()

parse_event_date <- function(x) {
  parsed <- suppressWarnings(ymd_hms(x, tz = "UTC"))
  parsed <- if_else(is.na(parsed), suppressWarnings(ymd(x)), parsed)
  as.Date(parsed)
}

parse_bool <- function(x) {
  stringr::str_to_lower(stringr::str_trim(as.character(x))) %in% c("true", "1", "yes", "y")
}

merge_intervals <- function(intervals_df) {
  if (nrow(intervals_df) == 0) {
    return(tibble(start_date = as.Date(character()), end_date = as.Date(character())))
  }

  df <- intervals_df %>%
    mutate(end_date = if_else(is.na(end_date), as.Date(NA), end_date)) %>%
    arrange(start_date, coalesce(end_date, start_date))

  merged_list <- list()
  current_start <- df$start_date[1]
  current_end <- df$end_date[1]
  current_open <- is.na(current_end)

  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      s <- df$start_date[i]
      e <- df$end_date[i]
      open_e <- is.na(e)
      overlaps <- current_open || (!is.na(current_end) && s <= current_end + 1)

      if (overlaps) {
        if (!current_open) {
          if (open_e) {
            current_open <- TRUE
            current_end <- NA
          } else {
            current_end <- max(current_end, e, na.rm = TRUE)
          }
        }
      } else {
        merged_list[[length(merged_list) + 1]] <- tibble(
          start_date = current_start,
          end_date = current_end
        )
        current_start <- s
        current_end <- e
        current_open <- open_e
      }
    }
  }

  merged_list[[length(merged_list) + 1]] <- tibble(
    start_date = current_start,
    end_date = current_end
  )

  bind_rows(merged_list)
}

events <- readr::read_csv("data/events.csv", show_col_types = FALSE)
cleaned_seizures <- readr::read_csv("output/tabs/seizures/seizures.csv", show_col_types = FALSE)
forms <- if (file.exists("data/forms.csv")) {
  readr::read_csv("data/forms.csv", show_col_types = FALSE)
} else {
  tibble()
}
med_intakes <- readr::read_csv("data/med_intakes.csv", show_col_types = FALSE)
medications_raw <- read_medications_corrected()
surveys <- readr::read_csv("data/prospective_surveys.csv", show_col_types = FALSE)
milestones_raw <- readr::read_csv("data/prospective_development_milestones.csv", show_col_types = FALSE)
whatsapp_status <- readr::read_csv("data/whatsapp_status.csv", show_col_types = FALSE)

patient_variant_lookup <- whatsapp_status %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    variant_p = stringr::str_squish(as.character(.data$variant_p)),
    variant_p = na_if(.data$variant_p, ""),
    run_timestamp = as.character(.data$run_timestamp),
    run_timestamp_parsed = suppressWarnings(ymd_hms(stringr::str_replace(.data$run_timestamp, "_", " "), tz = "UTC")),
    run_timestamp_sort = coalesce(.data$run_timestamp_parsed, as.POSIXct("1970-01-01 00:00:00", tz = "UTC"))
  ) %>%
  filter(.data$patient_id %in% TARGET_PATIENT_IDS) %>%
  arrange(.data$patient_id, desc(!is.na(.data$variant_p)), desc(.data$run_timestamp_sort)) %>%
  distinct(.data$patient_id, .keep_all = TRUE) %>%
  select(patient_id, variant_p)

seizure_events <- cleaned_seizures %>%
  mutate(
    patient_id = as.character(.data$patient_id),
    event_date = parse_event_date(.data$date)
  ) %>%
  filter(.data$patient_id %in% TARGET_PATIENT_IDS) %>%
  left_join(PATIENT_START_DATES, by = "patient_id") %>%
  filter(is.na(.data$min_event_date) | (!is.na(.data$event_date) & .data$event_date >= .data$min_event_date)) %>%
  filter(!is.na(.data$event_date), .data$event_date <= current_date) %>%
  select(patient_id, event_date)

app_activity_dates <- jsonlite::fromJSON("data/patient_summary_metrics.json") %>%
  as_tibble() %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    app_activity_date = parse_event_date(.data$first_app_activity_createdAt)
  ) %>%
  filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    is.na(.data$app_activity_date) | .data$app_activity_date <= current_date
  )

app_usage_parts <- list()

if ("date" %in% names(events)) {
  app_usage_parts <- append(
    app_usage_parts,
    list(
      events %>%
        transmute(patient_id = as.character(.data$patient_id), usage_date = parse_event_date(.data$date)) %>%
        filter(.data$patient_id %in% TARGET_PATIENT_IDS, !is.na(.data$usage_date), .data$usage_date <= current_date)
    )
  )
}

if (nrow(forms) > 0 && "date" %in% names(forms)) {
  app_usage_parts <- append(
    app_usage_parts,
    list(
      forms %>%
        transmute(patient_id = as.character(.data$patient_id), usage_date = parse_event_date(.data$date)) %>%
        filter(.data$patient_id %in% TARGET_PATIENT_IDS, !is.na(.data$usage_date), .data$usage_date <= current_date)
    )
  )
}

if (nrow(med_intakes) > 0 && all(c("taken", "taken_date") %in% names(med_intakes))) {
  app_usage_parts <- append(
    app_usage_parts,
    list(
      med_intakes %>%
        transmute(
          patient_id = as.character(.data$patient_id),
          taken_flag = parse_bool(.data$taken),
          usage_date = parse_event_date(.data$taken_date)
        ) %>%
        filter(.data$patient_id %in% TARGET_PATIENT_IDS, .data$taken_flag, !is.na(.data$usage_date), .data$usage_date <= current_date) %>%
        select(patient_id, usage_date)
    )
  )
}

app_usage_daily <- if (length(app_usage_parts) == 0) {
  tibble(patient_id = character(), usage_date = as.Date(character()))
} else {
  bind_rows(app_usage_parts) %>%
    distinct(patient_id, usage_date)
}

medications_for_analysis <- standardize_non_rescue_epilepsy_medications(
  medications_raw,
  include_non_drug = FALSE
) %>%
  mutate(name = .data$name_standardized)

med_intervals_raw <- read_medication_intervals_corrected(
  medications_for_analysis = medications_for_analysis
) %>%
  filter(
    .data$patient_id %in% TARGET_PATIENT_IDS
  ) %>%
  select(patient_id, medication_id, start_date, end_date) %>%
  distinct()

medication_lookup <- medications_for_analysis %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    medication_id = .data$medication_id,
    name = .data$name
  )

med_schedule_summary_raw <- med_intervals_raw %>%
  group_by(patient_id) %>%
  summarise(
    med_schedule_start_max = max(.data$start_date, na.rm = TRUE),
    .groups = "drop"
  )

med_intervals <- med_intervals_raw %>%
  left_join(medication_lookup, by = c("patient_id", "medication_id")) %>%
  mutate(
    med_name = if_else(is.na(.data$name) | .data$name == "", "Unknown medication", .data$name),
    med_name = stringr::str_squish(.data$med_name)
  ) %>%
  select(patient_id, med_name, start_date, end_date) %>%
  group_by(patient_id, med_name) %>%
  group_modify(~ merge_intervals(.x)) %>%
  ungroup()

milestone_labels <- c(
  eye_ps = "Eye Control",
  grasp_ps = "Grasp",
  reach_ps = "Reach",
  pincer_ps = "Pincer Grasp",
  blocks_ps = "Build Block Tower",
  circle_ps = "Draw Circle",
  hc_ps = "Head Control",
  roll_ps = "Roll Over",
  sit_ps = "Sit",
  stand_ps = "Stand Supported",
  walk_ps = "Walk",
  run_ps = "Run",
  smile_ps = "Smile",
  wave_ps = "Wave",
  cup_ps = "Drink Cup",
  fork_ps = "Use Utensils",
  wh_ps = "Wash Hands",
  bt_ps = "Brush Teeth",
  vocalize_ps = "Vocalizing",
  laugh_ps = "Laughing",
  babble_ps = "Babbling",
  namecolors_ps = "Name Colors",
  words_ps = "Use 2 Words",
  phrase_ps = "Phrases",
  reade_ps = "Read"
)

milestone_status_history <- milestones_raw %>%
  mutate(status_numeric = suppressWarnings(as.numeric(.data$status))) %>%
  left_join(
    surveys %>%
      transmute(
        survey_instance_id = .data$survey_instance_id,
        patient_id = as.character(.data$patient_id),
        prospective_study_timestamp_utc = .data$prospective_study_timestamp_utc,
        prospective_study_timestamp = .data$prospective_study_timestamp
      ),
    by = "survey_instance_id"
  ) %>%
  mutate(
    event_date = parse_event_date(.data$prospective_study_timestamp_utc),
    event_date = if_else(
      is.na(.data$event_date),
      as.Date(suppressWarnings(mdy_hm(.data$prospective_study_timestamp, tz = "UTC"))),
      .data$event_date
    ),
    milestone_label = recode(.data$milestone, !!!milestone_labels, .default = stringr::str_to_title(stringr::str_replace_all(.data$milestone, "_", " ")))
  ) %>%
  filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    !is.na(.data$event_date),
    .data$event_date <= current_date
  )

milestone_intervals <- tibble(
  patient_id = character(),
  milestone_label = character(),
  start_date = as.Date(character()),
  end_date = as.Date(character())
)

if (nrow(milestone_status_history) > 0) {
  milestone_runs <- milestone_status_history %>%
    arrange(patient_id, milestone_label, event_date) %>%
    group_by(patient_id, milestone_label) %>%
    mutate(
      status_flag = .data$status_numeric %in% c(2, 3, 5),
      change_flag = .data$status_flag != lag(.data$status_flag),
      run_id = cumsum(if_else(is.na(.data$change_flag) | .data$change_flag, 1L, 0L))
    ) %>%
    group_by(patient_id, milestone_label, run_id) %>%
    summarise(
      status_flag = first(.data$status_flag),
      run_start = min(.data$event_date),
      .groups = "drop"
    ) %>%
    arrange(patient_id, milestone_label, run_start) %>%
    group_by(patient_id, milestone_label) %>%
    mutate(next_run_start = lead(.data$run_start)) %>%
    ungroup()

  milestone_intervals <- milestone_runs %>%
    filter(.data$status_flag) %>%
    transmute(
      patient_id = .data$patient_id,
      milestone_label = .data$milestone_label,
      start_date = .data$run_start,
      end_date = coalesce(.data$next_run_start, current_date)
    )
}

seizure_summary <- seizure_events %>%
  group_by(patient_id) %>%
  summarise(
    seizure_first_date = min(.data$event_date, na.rm = TRUE),
    seizure_last_date = max(.data$event_date, na.rm = TRUE),
    .groups = "drop"
  )

med_summary <- med_intervals %>%
  group_by(patient_id) %>%
  summarise(
    med_end_max = if (all(is.na(.data$end_date))) as.Date(NA) else max(.data$end_date, na.rm = TRUE),
    .groups = "drop"
  )

milestone_summary <- milestone_intervals %>%
  group_by(patient_id) %>%
  summarise(
    milestone_last_date = max(.data$end_date, na.rm = TRUE),
    .groups = "drop"
  )

app_usage_summary <- app_usage_daily %>%
  group_by(patient_id) %>%
  summarise(
    app_usage_last_date = max(.data$usage_date, na.rm = TRUE),
    .groups = "drop"
  )

timeline_windows <- tibble(patient_id = TARGET_PATIENT_IDS) %>%
  left_join(patient_variant_lookup, by = "patient_id") %>%
  left_join(PATIENT_START_DATES, by = "patient_id") %>%
  left_join(app_activity_dates, by = "patient_id") %>%
  left_join(seizure_summary, by = "patient_id") %>%
  left_join(med_summary, by = "patient_id") %>%
  left_join(med_schedule_summary_raw, by = "patient_id") %>%
  left_join(milestone_summary, by = "patient_id") %>%
  left_join(app_usage_summary, by = "patient_id") %>%
  rowwise() %>%
  mutate(
    timeline_start_date = .data$app_activity_date,
    timeline_end_date = {
      dates <- c(
        .data$med_end_max,
        .data$med_schedule_start_max,
        .data$milestone_last_date,
        .data$app_usage_last_date
      )
      dates <- dates[!is.na(dates)]
      if (length(dates) == 0) as.Date(NA) else max(dates)
    }
  ) %>%
  ungroup() %>%
  mutate(
    timeline_end_date = if_else(is.na(.data$timeline_end_date), .data$timeline_end_date, pmin(.data$timeline_end_date, current_date)),
    study_start_date = case_when(
      !is.na(.data$timeline_start_date) & !is.na(.data$min_event_date) ~ pmax(.data$timeline_start_date, .data$min_event_date),
      !is.na(.data$timeline_start_date) ~ .data$timeline_start_date,
      !is.na(.data$min_event_date) ~ .data$min_event_date,
      TRUE ~ as.Date(NA)
    ),
    study_end_date = .data$timeline_end_date
  ) %>%
  mutate(
    study_start_date = if_else(!is.na(.data$study_end_date) & !is.na(.data$study_start_date) & .data$study_start_date > .data$study_end_date, as.Date(NA), .data$study_start_date)
  ) %>%
  select(
    patient_id,
    variant_p,
    seizure_first_date,
    seizure_last_date,
    study_start_date,
    study_end_date
  )

seizure_events_in_window <- seizure_events %>%
  left_join(
    timeline_windows %>% select(patient_id, study_start_date, study_end_date),
    by = "patient_id"
  ) %>%
  filter(
    !is.na(.data$study_start_date),
    !is.na(.data$study_end_date),
    .data$event_date >= .data$study_start_date,
    .data$event_date <= .data$study_end_date
  )

monthly_seizure_counts <- seizure_events_in_window %>%
  mutate(month = as.Date(floor_date(.data$event_date, "month"))) %>%
  count(patient_id, month, name = "seizure_count")

med_intervals_with_windows <- med_intervals %>%
  left_join(
    timeline_windows %>% select(patient_id, study_start_date, study_end_date),
    by = "patient_id"
  ) %>%
  filter(
    !is.na(.data$study_start_date),
    !is.na(.data$study_end_date),
    !is.na(.data$start_date)
  )

med_active_months <- med_intervals_with_windows %>%
  mutate(
    active_start_date = pmax(.data$start_date, .data$study_start_date),
    active_end_date = case_when(
      is.na(.data$end_date) ~ .data$study_end_date,
      TRUE ~ pmin(.data$end_date, .data$study_end_date)
    )
  ) %>%
  filter(!is.na(.data$active_end_date), .data$active_start_date <= .data$active_end_date) %>%
  mutate(
    month_start = as.Date(floor_date(.data$active_start_date, "month")),
    month_end = as.Date(floor_date(.data$active_end_date, "month")),
    month = purrr::map2(.data$month_start, .data$month_end, ~ seq(.x, .y, by = "month"))
  ) %>%
  select(patient_id, med_name, month) %>%
  tidyr::unnest(month) %>%
  distinct(patient_id, month, med_name)

monthly_active_med_counts <- med_active_months %>%
  count(patient_id, month, name = "active_med_count")

med_exposure_summary <- med_active_months %>%
  group_by(med_name) %>%
  summarise(
    n_patients_exposed = n_distinct(patient_id),
    n_patient_months_exposed = n(),
    .groups = "drop"
  ) %>%
  filter(
    n_patients_exposed >= MED_EXPOSURE_MIN_PATIENTS,
    n_patient_months_exposed >= MED_EXPOSURE_MIN_PATIENT_MONTHS
  ) %>%
  arrange(desc(n_patients_exposed), desc(n_patient_months_exposed), med_name) %>%
  mutate(
    med_col_base = stringr::str_to_lower(med_name),
    med_col_base = stringr::str_replace_all(med_col_base, "[^a-z0-9]+", "_"),
    med_col_base = stringr::str_replace_all(med_col_base, "^_+|_+$", ""),
    med_col_base = if_else(med_col_base == "", "unknown", med_col_base),
    med_col = paste0("med_exposed_", med_col_base)
  ) %>%
  group_by(med_col) %>%
  mutate(
    med_col = if_else(dplyr::n() > 1, paste0(med_col, "_", row_number()), med_col)
  ) %>%
  ungroup() %>%
  select(med_name, med_col)

med_specific_cols <- med_exposure_summary$med_col

med_month_exposure <- if (length(med_specific_cols) == 0) {
  tibble(patient_id = character(), month = as.Date(character()), med_col = character())
} else {
  med_active_months %>%
    inner_join(med_exposure_summary, by = "med_name") %>%
    distinct(patient_id, month, med_col)
}

monthly_med_specific <- if (length(med_specific_cols) == 0) {
  tibble(patient_id = character(), month = as.Date(character()))
} else {
  med_month_exposure %>%
    mutate(exposed = 1L) %>%
    tidyr::pivot_wider(
      id_cols = c(patient_id, month),
      names_from = med_col,
      values_from = exposed,
      values_fill = 0L
    )
}

single_med_month_counts <- med_month_exposure %>%
  count(med_col, name = "n_months_exposed")

med_pair_month_exposure <- if (length(med_specific_cols) < 2) {
  tibble(
    patient_id = character(),
    month = as.Date(character()),
    med_col_a = character(),
    med_col_b = character()
  )
} else {
  med_month_exposure %>%
    inner_join(
      med_month_exposure,
      by = c("patient_id", "month"),
      suffix = c("_a", "_b"),
      relationship = "many-to-many"
    ) %>%
    filter(.data$med_col_a < .data$med_col_b) %>%
    distinct(patient_id, month, med_col_a, med_col_b)
}

med_pair_exposure_summary <- if (nrow(med_pair_month_exposure) == 0) {
  tibble(
    med_col_a = character(),
    med_col_b = character(),
    pair_col = character()
  )
} else {
  med_pair_month_exposure %>%
    group_by(med_col_a, med_col_b) %>%
    summarise(
      n_pair_months_exposed = n(),
      n_patients_pair_exposed = n_distinct(patient_id),
      .groups = "drop"
    ) %>%
    left_join(
      single_med_month_counts %>% rename(med_col_a = med_col, n_months_a = n_months_exposed),
      by = "med_col_a"
    ) %>%
    left_join(
      single_med_month_counts %>% rename(med_col_b = med_col, n_months_b = n_months_exposed),
      by = "med_col_b"
    ) %>%
    mutate(
      n_months_a_only = .data$n_months_a - .data$n_pair_months_exposed,
      n_months_b_only = .data$n_months_b - .data$n_pair_months_exposed
    ) %>%
    filter(
      .data$n_pair_months_exposed >= MED_PAIR_EXPOSURE_MIN_PATIENT_MONTHS,
      .data$n_patients_pair_exposed >= MED_PAIR_EXPOSURE_MIN_PATIENTS,
      .data$n_months_a_only >= MED_PAIR_EXPOSURE_MIN_SINGLE_ONLY_MONTHS,
      .data$n_months_b_only >= MED_PAIR_EXPOSURE_MIN_SINGLE_ONLY_MONTHS
    ) %>%
    arrange(desc(.data$n_pair_months_exposed), desc(.data$n_patients_pair_exposed), .data$med_col_a, .data$med_col_b) %>%
    mutate(
      pair_col = paste0(
        "med_pair_exposed_",
        stringr::str_remove(.data$med_col_a, "^med_exposed_"),
        "__",
        stringr::str_remove(.data$med_col_b, "^med_exposed_")
      )
    ) %>%
    select(med_col_a, med_col_b, pair_col)
}

med_pair_specific_cols <- med_pair_exposure_summary$pair_col

monthly_med_pair_specific <- if (length(med_pair_specific_cols) == 0) {
  tibble(patient_id = character(), month = as.Date(character()))
} else {
  med_pair_month_exposure %>%
    inner_join(med_pair_exposure_summary, by = c("med_col_a", "med_col_b")) %>%
    select(patient_id, month, pair_col) %>%
    distinct(patient_id, month, pair_col) %>%
    mutate(exposed = 1L) %>%
    tidyr::pivot_wider(
      id_cols = c(patient_id, month),
      names_from = pair_col,
      values_from = exposed,
      values_fill = 0L
    )
}

monthly_med_start_flags <- med_intervals_with_windows %>%
  filter(
    .data$start_date >= .data$study_start_date,
    .data$start_date <= .data$study_end_date
  ) %>%
  transmute(
    patient_id,
    month = as.Date(floor_date(.data$start_date, "month")),
    med_name
  ) %>%
  distinct(patient_id, month, med_name) %>%
  count(patient_id, month, name = "n_meds_started") %>%
  mutate(med_started_flag = as.integer(.data$n_meds_started > 0L)) %>%
  select(patient_id, month, med_started_flag)

monthly_med_stop_flags <- med_intervals_with_windows %>%
  filter(
    !is.na(.data$end_date),
    .data$end_date >= .data$study_start_date,
    .data$end_date <= .data$study_end_date
  ) %>%
  transmute(
    patient_id,
    month = as.Date(floor_date(.data$end_date, "month")),
    med_name
  ) %>%
  distinct(patient_id, month, med_name) %>%
  count(patient_id, month, name = "n_meds_stopped") %>%
  mutate(med_stopped_flag = as.integer(.data$n_meds_stopped > 0L)) %>%
  select(patient_id, month, med_stopped_flag)

monthly_med_flags <- full_join(monthly_med_start_flags, monthly_med_stop_flags, by = c("patient_id", "month")) %>%
  mutate(
    med_started_flag = replace_na(.data$med_started_flag, 0L),
    med_stopped_flag = replace_na(.data$med_stopped_flag, 0L)
  )

patient_month_panel <- timeline_windows %>%
  filter(!is.na(.data$study_start_date), !is.na(.data$study_end_date), .data$study_start_date <= .data$study_end_date) %>%
  mutate(
    start_month = as.Date(floor_date(.data$study_start_date, "month")),
    end_month = as.Date(floor_date(.data$study_end_date, "month")),
    month = purrr::map2(.data$start_month, .data$end_month, ~ seq(.x, .y, by = "month"))
  ) %>%
  select(-start_month, -end_month) %>%
  tidyr::unnest(month) %>%
  left_join(monthly_seizure_counts, by = c("patient_id", "month")) %>%
  left_join(monthly_active_med_counts, by = c("patient_id", "month")) %>%
  left_join(monthly_med_flags, by = c("patient_id", "month")) %>%
  left_join(monthly_med_specific, by = c("patient_id", "month")) %>%
  left_join(monthly_med_pair_specific, by = c("patient_id", "month")) %>%
  mutate(
    seizure_count = replace_na(.data$seizure_count, 0L),
    active_med_count = replace_na(.data$active_med_count, 0L),
    med_started_flag = replace_na(.data$med_started_flag, 0L),
    med_stopped_flag = replace_na(.data$med_stopped_flag, 0L)
  ) %>%
  mutate(
    across(all_of(med_specific_cols), ~ replace_na(.x, 0L)),
    across(all_of(med_pair_specific_cols), ~ replace_na(.x, 0L))
  ) %>%
  group_by(.data$patient_id) %>%
  arrange(.data$month, .by_group = TRUE) %>%
  mutate(
    month_index = row_number() - 1L,
    lag_active_med_count = lag(.data$active_med_count, n = 1L),
    lag_med_started_flag = lag(.data$med_started_flag, n = 1L),
    lag_med_stopped_flag = lag(.data$med_stopped_flag, n = 1L)
  ) %>%
  ungroup() %>%
  arrange(.data$patient_id, .data$month_index)

readr::write_csv(patient_month_panel, OUTPUT_PANEL_PATH)
