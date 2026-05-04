# Create timeline plots for Helpilepsy patients

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

source("src/desc/medication_standardization.R")
source("src/desc/seizure_type_standardization.R")

dir.create("output/figs", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figs/timelines", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tabs", recursive = TRUE, showWarnings = FALSE)
existing_timeline_pdfs <- list.files("output/figs/timelines", pattern = "\\.pdf$", full.names = TRUE)
if (length(existing_timeline_pdfs) > 0) {
  invisible(file.remove(existing_timeline_pdfs))
}
if (file.exists("output/figs/helpilepsy_timelines.pdf")) {
  invisible(file.remove("output/figs/helpilepsy_timelines.pdf"))
}
current_date <- Sys.Date()
CHANGE_POINT_INPUT_PATH <- "output/tabs/changepoints/patient_change_points.csv"
CHANGE_POINT_REQUIRED_COLUMNS <- c(
  "patient_id",
  "candidate_week",
  "significant"
)

# Milestone labels (from Citizen timelines)
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

# Helper to parse date formats
parse_event_date <- function(x) {
  parsed <- suppressWarnings(ymd_hms(x, tz = "UTC"))
  parsed <- if_else(is.na(parsed), suppressWarnings(ymd(x)), parsed)
  as.Date(parsed)
}

parse_bool <- function(x) {
  str_to_lower(str_trim(as.character(x))) %in% c("true", "1", "yes", "y")
}

parse_yes_no <- function(x) {
  normalized <- str_to_lower(str_squish(as.character(x)))
  case_when(
    normalized %in% c("yes", "y", "true", "1") ~ "yes",
    normalized %in% c("no", "n", "false", "0") ~ "no",
    TRUE ~ NA_character_
  )
}

parse_change_point_flag <- function(x) {
  str_to_lower(str_squish(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

change_point_date_column <- function(data, column_name) {
  if (column_name %in% names(data)) {
    parse_event_date(data[[column_name]])
  } else {
    as.Date(rep(NA, nrow(data)))
  }
}

change_point_character_column <- function(data, column_name) {
  if (column_name %in% names(data)) {
    as.character(data[[column_name]])
  } else {
    rep(NA_character_, nrow(data))
  }
}

min_date_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) as.Date(NA) else min(x)
}

max_date_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) as.Date(NA) else max(x)
}

# Merge overlapping medication intervals
merge_intervals <- function(intervals_df) {
  if (nrow(intervals_df) == 0) {
    return(tibble(start_date = as.Date(character()), end_date = as.Date(character()), is_ongoing = logical()))
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
          end_date = current_end,
          is_ongoing = current_open
        )
        current_start <- s
        current_end <- e
        current_open <- open_e
      }
    }
  }

  merged_list[[length(merged_list) + 1]] <- tibble(
    start_date = current_start,
    end_date = current_end,
    is_ongoing = current_open
  )

  bind_rows(merged_list)
}

# Import data

patients <- readr::read_csv("data/patients.csv", show_col_types = FALSE) %>%
  select(patient_id, first_name, last_name)

events <- readr::read_csv("data/events.csv", show_col_types = FALSE)
medications_raw <- readr::read_csv("data/medications.csv", show_col_types = FALSE)
medications_for_analysis <- standardize_non_rescue_epilepsy_medications(
  medications_raw,
  include_non_drug = FALSE
) %>%
  mutate(name = name_standardized)
forms <- if (file.exists("data/forms.csv")) {
  readr::read_csv("data/forms.csv", show_col_types = FALSE)
} else {
  tibble()
}
form_answers <- if (file.exists("data/form_answers.csv")) {
  readr::read_csv("data/form_answers.csv", show_col_types = FALSE)
} else {
  tibble()
}
med_dosages <- readr::read_csv("data/med_dosages.csv", show_col_types = FALSE)
med_intakes <- readr::read_csv("data/med_intakes.csv", show_col_types = FALSE)
surveys <- readr::read_csv("data/prospective_surveys.csv", show_col_types = FALSE)
milestones_raw <- readr::read_csv("data/prospective_development_milestones.csv", show_col_types = FALSE)
app_activity_dates <- jsonlite::fromJSON("data/patient_summary_metrics.json") %>%
  as_tibble() %>%
  transmute(
    patient_id,
    app_activity_date = parse_event_date(first_app_activity_createdAt)
  )

change_point_markers <- if (file.exists(CHANGE_POINT_INPUT_PATH)) {
  change_point_raw <- readr::read_csv(
    CHANGE_POINT_INPUT_PATH,
    col_types = readr::cols(.default = readr::col_guess(), patient_id = readr::col_character())
  )
  missing_change_point_columns <- setdiff(CHANGE_POINT_REQUIRED_COLUMNS, names(change_point_raw))
  if (length(missing_change_point_columns) > 0) {
    stop(
      "Change-point output is missing required columns: ",
      paste(missing_change_point_columns, collapse = ", ")
    )
  }

  change_point_raw %>%
    transmute(
      patient_id = as.character(patient_id),
      change_point_date = parse_event_date(candidate_week),
      pre_segment_start_date = change_point_date_column(change_point_raw, "pre_segment_start_date"),
      post_segment_end_date = change_point_date_column(change_point_raw, "post_segment_end_date"),
      change_direction = str_to_lower(str_squish(change_point_character_column(change_point_raw, "direction"))),
      significant = parse_change_point_flag(significant)
    ) %>%
    filter(
      significant,
      !is.na(patient_id),
      !is.na(change_point_date)
    ) %>%
    distinct(patient_id, change_point_date, pre_segment_start_date, post_segment_end_date, change_direction)
} else {
  tibble(
    patient_id = character(),
    pre_segment_start_date = as.Date(character()),
    change_point_date = as.Date(character()),
    post_segment_end_date = as.Date(character()),
    change_direction = character()
  )
}

prospective_survey_completion <- surveys %>%
  mutate(
    survey_complete_flag = suppressWarnings(as.numeric(prospective_study_complete)) == 2,
    survey_date = parse_event_date(prospective_study_timestamp_utc),
    survey_date = if_else(
      is.na(survey_date),
      as.Date(suppressWarnings(mdy_hm(prospective_study_timestamp, tz = "UTC"))),
      survey_date
    )
  ) %>%
  filter(!is.na(patient_id), survey_complete_flag, !is.na(survey_date)) %>%
  distinct(patient_id)

# Daily app usage

app_usage_parts <- list()

if ("date" %in% names(events)) {
  app_usage_parts <- append(
    app_usage_parts,
    list(
      events %>%
        transmute(patient_id, usage_date = parse_event_date(date)) %>%
        filter(!is.na(patient_id), !is.na(usage_date), usage_date <= current_date)
    )
  )
}

if (nrow(forms) > 0 && "date" %in% names(forms)) {
  app_usage_parts <- append(
    app_usage_parts,
    list(
      forms %>%
        transmute(patient_id, usage_date = parse_event_date(date)) %>%
        filter(!is.na(patient_id), !is.na(usage_date), usage_date <= current_date)
    )
  )
}

if (nrow(med_intakes) > 0 && all(c("taken", "taken_date") %in% names(med_intakes))) {
  app_usage_parts <- append(
    app_usage_parts,
    list(
      med_intakes %>%
        transmute(
          patient_id,
          taken_flag = parse_bool(taken),
          usage_date = parse_event_date(taken_date)
        ) %>%
        filter(taken_flag, !is.na(patient_id), !is.na(usage_date), usage_date <= current_date) %>%
        select(patient_id, usage_date)
    )
  )
}

app_usage_daily <- if (length(app_usage_parts) == 0) {
  tibble(patient_id = character(), usage_date = as.Date(character()), app_actions = integer())
} else {
  bind_rows(app_usage_parts) %>%
    distinct(patient_id, usage_date) %>%
    mutate(app_actions = 1L)
}

# Weekly survey completion status

scn8a_diary_q1_code <- "q1_1_did_you_give_your_child_all_medications_every_day_this_week"
scn8a_diary_q2_code <- "q2_2_does_the_daily_seizure_record_in_the_app_this_week_accurately_reflect_your_child_s_seizure_count_including_if_there_were_none"

scn8a_diary_forms <- if (
  nrow(forms) > 0 &&
    all(c("form_id", "patient_id", "form_name", "date") %in% names(forms))
) {
  forms %>%
    mutate(
      form_name_clean = str_to_lower(str_squish(form_name)),
      diary_date = parse_event_date(date)
    ) %>%
    filter(
      !is.na(form_id),
      !is.na(patient_id),
      form_name_clean == "scn8a diary completion",
      !is.na(diary_date),
      diary_date <= current_date
    ) %>%
    select(form_id, patient_id, diary_date)
} else {
  tibble(
    form_id = character(),
    patient_id = character(),
    diary_date = as.Date(character())
  )
}

diary_answers_wide <- if (
  nrow(scn8a_diary_forms) > 0 &&
    nrow(form_answers) > 0 &&
    all(c("form_id", "question_code", "answer") %in% names(form_answers))
) {
  form_answers %>%
    semi_join(scn8a_diary_forms %>% select(form_id), by = "form_id") %>%
    filter(question_code %in% c(scn8a_diary_q1_code, scn8a_diary_q2_code)) %>%
    mutate(answer_yes_no = parse_yes_no(answer)) %>%
    arrange(form_id, question_code) %>%
    group_by(form_id, question_code) %>%
    summarise(answer_yes_no = first(answer_yes_no[!is.na(answer_yes_no)]), .groups = "drop") %>%
    pivot_wider(
      id_cols = form_id,
      names_from = question_code,
      values_from = answer_yes_no
    )
} else {
  tibble(
    form_id = character(),
    !!scn8a_diary_q1_code := character(),
    !!scn8a_diary_q2_code := character()
  )
}

scn8a_diary_events <- scn8a_diary_forms %>%
  left_join(diary_answers_wide, by = "form_id") %>%
  mutate(
    diary_status = case_when(
      .data[[scn8a_diary_q1_code]] == "yes" & .data[[scn8a_diary_q2_code]] == "yes" ~ "All yes",
      .data[[scn8a_diary_q1_code]] == "no" | .data[[scn8a_diary_q2_code]] == "no" ~ "Any no",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(diary_status))

# Seizure events

seizure_events <- events %>%
  standardize_seizure_events(filter_to_seizure = TRUE) %>%
  mutate(
    event_date = coalesce(event_date, parse_event_date(date)),
    seizure_label = if_else(
      is.na(seizure_type_standardized) | seizure_type_standardized == "",
      "Unspecified",
      seizure_type_standardized
    )
  ) %>%
  filter(!is.na(event_date)) %>%
  select(
    patient_id, event_date, seizure_label, seizure_group,
    seizure_type_raw, seizure_type_standardized, seizure_type_primary
  )

# Medication intervals

med_intervals_raw <- bind_rows(
  med_dosages %>%
    transmute(patient_id, medication_id, start_date = ymd(from), end_date = ymd(to)),
  med_intakes %>%
    transmute(patient_id, medication_id, start_date = ymd(intake_from), end_date = ymd(intake_to))
) %>%
  mutate(end_date = if_else(!is.na(end_date) & end_date < start_date, as.Date(NA), end_date)) %>%
  filter(!is.na(start_date)) %>%
  semi_join(
    medications_for_analysis %>% select(patient_id, medication_id),
    by = c("patient_id", "medication_id")
  ) %>%
  distinct()

medication_lookup <- medications_for_analysis %>%
  transmute(
    patient_id,
    medication_id,
    name
  )

medication_refs_missing_metadata <- med_intervals_raw %>%
  anti_join(medication_lookup %>% select(patient_id, medication_id), by = c("patient_id", "medication_id")) %>%
  distinct(patient_id, medication_id) %>%
  arrange(patient_id, medication_id)

med_schedule_summary_raw <- med_intervals_raw %>%
  group_by(patient_id) %>%
  summarise(
    med_schedule_start_min = min(start_date, na.rm = TRUE),
    med_schedule_start_max = max(start_date, na.rm = TRUE),
    med_schedule_end_max = if (all(is.na(end_date))) as.Date(NA) else max(end_date, na.rm = TRUE),
    .groups = "drop"
  )

med_intervals <- med_intervals_raw %>%
  left_join(medication_lookup, by = c("patient_id", "medication_id")) %>%
  mutate(
    med_name = if_else(is.na(name) | name == "", "Unknown medication", name),
    med_name = str_squish(med_name)
  ) %>%
  select(patient_id, med_name, start_date, end_date) %>%
  group_by(patient_id, med_name) %>%
  group_modify(~ merge_intervals(.x)) %>%
  ungroup()

# Developmental milestones

milestone_status_history <- milestones_raw %>%
  mutate(status_numeric = suppressWarnings(as.numeric(status))) %>%
  left_join(
    surveys %>% select(survey_instance_id, patient_id, prospective_study_timestamp_utc, prospective_study_timestamp),
    by = "survey_instance_id"
  ) %>%
  mutate(
    event_date = parse_event_date(prospective_study_timestamp_utc),
    event_date = if_else(
      is.na(event_date),
      as.Date(suppressWarnings(mdy_hm(prospective_study_timestamp, tz = "UTC"))),
      event_date
    ),
    milestone_label = recode(milestone, !!!milestone_labels, .default = str_to_title(str_replace_all(milestone, "_", " ")))
  ) %>%
  filter(!is.na(patient_id), !is.na(event_date))

development_module_answered <- milestone_status_history %>%
  distinct(patient_id)

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
      status_flag = status_numeric %in% c(2, 3, 5),
      change_flag = status_flag != lag(status_flag),
      run_id = cumsum(if_else(is.na(change_flag) | change_flag, 1L, 0L))
    ) %>%
    group_by(patient_id, milestone_label, run_id) %>%
    summarise(
      status_flag = first(status_flag),
      run_start = min(event_date),
      .groups = "drop"
    ) %>%
    arrange(patient_id, milestone_label, run_start) %>%
    group_by(patient_id, milestone_label) %>%
    mutate(next_run_start = lead(run_start)) %>%
    ungroup()

  milestone_intervals <- milestone_runs %>%
    filter(status_flag) %>%
    transmute(
      patient_id,
      milestone_label,
      start_date = run_start,
      end_date = coalesce(next_run_start, current_date)
    )
}

# Patient names + variants

whatsapp_patient_metadata <- readr::read_csv("data/whatsapp_status.csv", show_col_types = FALSE) %>%
  mutate(
    patient_id = as.character(patient_id),
    whatsapp_name = str_squish(str_trim(paste(first_name, last_name))),
    whatsapp_name = na_if(whatsapp_name, ""),
    variant_p = str_squish(as.character(variant_p)),
    variant_p = na_if(variant_p, ""),
    run_timestamp = as.character(run_timestamp)
  ) %>%
  filter(!is.na(patient_id), patient_id != "") %>%
  arrange(patient_id, desc(!is.na(variant_p)), desc(run_timestamp)) %>%
  group_by(patient_id) %>%
  summarise(
    whatsapp_name = first(whatsapp_name[!is.na(whatsapp_name)]),
    variant_p = first(variant_p),
    .groups = "drop"
  )

whatsapp_names <- whatsapp_patient_metadata %>%
  select(patient_id, whatsapp_name) %>%
  filter(!is.na(whatsapp_name), whatsapp_name != "")

patients_with_data <- union(
  union(unique(seizure_events$patient_id), unique(med_intervals$patient_id)),
  unique(milestone_intervals$patient_id)
)
patients_with_data <- union(patients_with_data, unique(prospective_survey_completion$patient_id))
patients_with_data <- union(patients_with_data, unique(app_usage_daily$patient_id))
patients_with_data <- union(patients_with_data, unique(scn8a_diary_events$patient_id))
patients_with_data <- patients_with_data[!is.na(patients_with_data)]
patients_with_data <- intersect(patients_with_data, whatsapp_names$patient_id)

milestone_row_counts <- milestone_intervals %>%
  group_by(patient_id) %>%
  summarise(n_milestone_rows = n_distinct(milestone_label), .groups = "drop")
seizure_row_counts <- seizure_events %>%
  group_by(patient_id) %>%
  summarise(n_seizure_rows = n_distinct(seizure_label), .groups = "drop")
med_row_counts <- med_intervals %>%
  group_by(patient_id) %>%
  summarise(n_med_rows = n_distinct(med_name), .groups = "drop")
diary_row_counts <- scn8a_diary_events %>%
  group_by(patient_id) %>%
  summarise(n_diary_rows = 1L, .groups = "drop")
patient_row_slots <- tibble(patient_id = patients_with_data) %>%
  left_join(milestone_row_counts, by = "patient_id") %>%
  left_join(seizure_row_counts, by = "patient_id") %>%
  left_join(med_row_counts, by = "patient_id") %>%
  left_join(diary_row_counts, by = "patient_id") %>%
  mutate(
    n_milestone_rows = replace_na(n_milestone_rows, 0L),
    n_seizure_rows = replace_na(n_seizure_rows, 0L),
    n_med_rows = replace_na(n_med_rows, 0L),
    n_diary_rows = replace_na(n_diary_rows, 0L),
    n_app_usage_rows = 1L,
    n_no_skills_rows = if_else(
      patient_id %in% prospective_survey_completion$patient_id & n_milestone_rows == 0L,
      1L,
      0L
    ),
    n_total_y_rows = n_milestone_rows + n_no_skills_rows + n_diary_rows + n_app_usage_rows + n_seizure_rows + n_med_rows
  )

max_y_rows <- max(patient_row_slots$n_total_y_rows, na.rm = TRUE)
if (!is.finite(max_y_rows) || max_y_rows < 1) {
  max_y_rows <- 1L
}

# Outputs for quality checking

patient_order <- tibble(
  page = seq_along(patients_with_data),
  patient_id = patients_with_data
)

patient_names <- patients %>%
  mutate(patient_name = str_squish(str_trim(paste(first_name, last_name)))) %>%
  group_by(patient_id) %>%
  summarise(patient_name = first(patient_name[patient_name != ""]), .groups = "drop")
patient_names <- patient_names %>%
  full_join(whatsapp_patient_metadata %>% select(patient_id, whatsapp_name, variant_p), by = "patient_id") %>%
  mutate(
    display_name = coalesce(whatsapp_name, patient_name),
    display_name = if_else(is.na(display_name) | display_name == "", patient_id, display_name),
    variant_p = na_if(str_squish(as.character(variant_p)), "")
  ) %>%
  select(patient_id, display_name, variant_p)

patient_timeline_filenames <- patient_names %>%
  mutate(
    timeline_file_base = str_to_lower(str_replace_all(coalesce(display_name, ""), "[^A-Za-z0-9]+", "_")),
    timeline_file_base = str_replace_all(timeline_file_base, "^_+|_+$", ""),
    timeline_file_base = if_else(timeline_file_base == "", patient_id, timeline_file_base)
  ) %>%
  arrange(patient_id) %>%
  group_by(timeline_file_base) %>%
  mutate(
    timeline_file_base = if_else(n() > 1, paste0(timeline_file_base, "_", row_number()), timeline_file_base)
  ) %>%
  ungroup() %>%
  mutate(timeline_filename = paste0(timeline_file_base, "_timeline.pdf")) %>%
  select(patient_id, timeline_filename)

seizure_summary <- seizure_events %>%
  group_by(patient_id) %>%
  summarise(
    n_seizure_events = n(),
    n_seizure_types = n_distinct(seizure_label),
    seizure_first_date = min(event_date, na.rm = TRUE),
    seizure_last_date = max(event_date, na.rm = TRUE),
    .groups = "drop"
  )

med_summary <- med_intervals %>%
  group_by(patient_id) %>%
  summarise(
    n_med_intervals = n(),
    n_medications = n_distinct(med_name),
    n_medications_ongoing = sum(is.na(end_date)),
    med_start_min = min(start_date, na.rm = TRUE),
    med_start_max = max(start_date, na.rm = TRUE),
    med_end_max = if (all(is.na(end_date))) as.Date(NA) else max(end_date, na.rm = TRUE),
    .groups = "drop"
  )

milestone_summary <- milestone_intervals %>%
  group_by(patient_id) %>%
  summarise(
    n_milestones_achieved = n_distinct(milestone_label),
    milestone_first_date = min(start_date, na.rm = TRUE),
    milestone_last_date = max(end_date, na.rm = TRUE),
    .groups = "drop"
  )
app_usage_summary <- app_usage_daily %>%
  group_by(patient_id) %>%
  summarise(
    n_app_usage_days = n(),
    app_usage_first_date = min(usage_date, na.rm = TRUE),
    app_usage_last_date = max(usage_date, na.rm = TRUE),
    .groups = "drop"
  )
diary_summary <- scn8a_diary_events %>%
  group_by(patient_id) %>%
  summarise(
    n_diary_submissions = n(),
    n_diary_all_yes = sum(diary_status == "All yes"),
    n_diary_any_no = sum(diary_status == "Any no"),
    diary_first_date = min(diary_date, na.rm = TRUE),
    diary_last_date = max(diary_date, na.rm = TRUE),
    .groups = "drop"
  )
change_point_summary <- change_point_markers %>%
  group_by(patient_id) %>%
  summarise(
    n_significant_change_points = n(),
    change_point_first_date = min_date_or_na(change_point_date),
    pre_segment_start_min = min_date_or_na(pre_segment_start_date),
    post_segment_end_max = max_date_or_na(post_segment_end_date),
    .groups = "drop"
  )

medications_missing_schedule <- medications_for_analysis %>%
  anti_join(med_intervals_raw %>% distinct(patient_id, medication_id), by = c("patient_id", "medication_id")) %>%
  select(
    patient_id, medication_id, name, standardized_components, standardization_status,
    requires_review, review_reason, reason, treatment_type, intake_type,
    is_deleted, createdAt, updatedAt
  ) %>%
  arrange(patient_id, name)

medications_missing_schedule_summary <- medications_missing_schedule %>%
  group_by(patient_id) %>%
  summarise(n_med_records_missing_schedule = n_distinct(medication_id), .groups = "drop")

medication_refs_missing_metadata_summary <- medication_refs_missing_metadata %>%
  group_by(patient_id) %>%
  summarise(n_med_refs_missing_metadata = n_distinct(medication_id), .groups = "drop")

timeline_qc <- patient_order %>%
  left_join(patient_names, by = "patient_id") %>%
  left_join(app_activity_dates, by = "patient_id") %>%
  left_join(seizure_summary, by = "patient_id") %>%
  left_join(med_summary, by = "patient_id") %>%
  left_join(med_schedule_summary_raw, by = "patient_id") %>%
  left_join(milestone_summary, by = "patient_id") %>%
  left_join(app_usage_summary, by = "patient_id") %>%
  left_join(diary_summary, by = "patient_id") %>%
  left_join(change_point_summary, by = "patient_id") %>%
  left_join(medications_missing_schedule_summary, by = "patient_id") %>%
  left_join(medication_refs_missing_metadata_summary, by = "patient_id") %>%
  mutate(
    n_seizure_events = replace_na(n_seizure_events, 0L),
    n_seizure_types = replace_na(n_seizure_types, 0L),
    n_med_intervals = replace_na(n_med_intervals, 0L),
    n_medications = replace_na(n_medications, 0L),
    n_medications_ongoing = replace_na(n_medications_ongoing, 0L),
    n_milestones_achieved = replace_na(n_milestones_achieved, 0L),
    n_app_usage_days = replace_na(n_app_usage_days, 0L),
    n_diary_submissions = replace_na(n_diary_submissions, 0L),
    n_diary_all_yes = replace_na(n_diary_all_yes, 0L),
    n_diary_any_no = replace_na(n_diary_any_no, 0L),
    n_significant_change_points = replace_na(n_significant_change_points, 0L),
    n_med_records_missing_schedule = replace_na(n_med_records_missing_schedule, 0L),
    n_med_refs_missing_metadata = replace_na(n_med_refs_missing_metadata, 0L)
  ) %>%
  rowwise() %>%
  mutate(
    timeline_start_date = {
      dates <- c(
        seizure_first_date, med_start_min, med_schedule_start_min, milestone_first_date,
        app_usage_first_date, diary_first_date, pre_segment_start_min, change_point_first_date
      )
      dates <- dates[!is.na(dates)]
      start_date <- if (length(dates) == 0) as.Date(NA) else min(dates)
      if (!is.na(app_activity_date) && n_significant_change_points == 0L) {
        max(start_date, app_activity_date)
      } else {
        start_date
      }
    },
    timeline_end_date = {
      dates <- c(
        seizure_last_date, med_end_max, med_start_max, med_schedule_start_max,
        med_schedule_end_max, milestone_last_date, app_usage_last_date,
        diary_last_date, change_point_first_date, post_segment_end_max
      )
      dates <- dates[!is.na(dates)]
      if (length(dates) == 0) as.Date(NA) else max(dates)
    }
  ) %>%
  ungroup() %>%
  mutate(
    timeline_end_date = if_else(is.na(timeline_end_date), timeline_end_date, pmin(timeline_end_date, current_date))
  )

readr::write_csv(timeline_qc, "output/tabs/helpilepsy_timeline_qc.csv")
readr::write_csv(medications_missing_schedule, "output/tabs/helpilepsy_medications_missing_schedule.csv")
readr::write_csv(medication_refs_missing_metadata, "output/tabs/helpilepsy_medication_refs_missing_metadata.csv")

# Plotting

plot_patient_timeline <- function(pt_id) {
  pt_info <- patients %>% filter(patient_id == pt_id) %>% slice_head(n = 1)
  pt_name <- patient_names %>% filter(patient_id == pt_id) %>% pull(display_name) %>% first()
  pt_variant <- patient_names %>% filter(patient_id == pt_id) %>% pull(variant_p) %>% first()
  if (is.na(pt_name) || pt_name == "") {
    pt_name <- str_trim(paste(pt_info$first_name, pt_info$last_name))
  }
  variant_label <- ifelse(is.na(pt_variant) || pt_variant == "", "Unknown", pt_variant)
  base_title <- ifelse(is.na(pt_name) || pt_name == "", pt_id, paste0(pt_name, " (", pt_id, ")"))
  title_text <- paste0(base_title, " ", variant_label)

  pt_meds <- med_intervals %>% filter(patient_id == pt_id)
  pt_med_schedule <- med_schedule_summary_raw %>% filter(patient_id == pt_id) %>% slice_head(n = 1)
  pt_seizures <- seizure_events %>% filter(patient_id == pt_id)
  pt_milestones <- milestone_intervals %>% filter(patient_id == pt_id)
  pt_app_usage <- app_usage_daily %>% filter(patient_id == pt_id)
  pt_diary <- scn8a_diary_events %>% filter(patient_id == pt_id)
  pt_change_points <- change_point_markers %>% filter(patient_id == pt_id)
  pt_completed_prospective <- prospective_survey_completion %>% filter(patient_id == pt_id) %>%
    pull(patient_id) %>% length() > 0
  pt_entered_development_module <- development_module_answered %>% filter(patient_id == pt_id) %>%
    pull(patient_id) %>% length() > 0
  show_no_skills_label <- pt_completed_prospective && pt_entered_development_module && nrow(pt_milestones) == 0
  show_no_entered_skills_label <- pt_completed_prospective && !pt_entered_development_module && nrow(pt_milestones) == 0
  pt_app_activity_date <- app_activity_dates %>% filter(patient_id == pt_id) %>%
    pull(app_activity_date) %>% first()

  if (
    nrow(pt_meds) == 0 &&
      nrow(pt_seizures) == 0 &&
      nrow(pt_milestones) == 0 &&
      nrow(pt_app_usage) == 0 &&
      nrow(pt_diary) == 0 &&
      nrow(pt_change_points) == 0 &&
      !show_no_skills_label &&
      !show_no_entered_skills_label
  ) {
    return(NULL)
  }

  pt_med_schedule_start_min <- pt_med_schedule$med_schedule_start_min %>% first()
  pt_med_schedule_start_max <- pt_med_schedule$med_schedule_start_max %>% first()
  pt_med_schedule_end_max <- pt_med_schedule$med_schedule_end_max %>% first()

  latest_date <- max(
    c(
      pt_meds$end_date, pt_meds$start_date,
      pt_med_schedule_start_max, pt_med_schedule_end_max,
      pt_seizures$event_date,
      pt_milestones$start_date, pt_milestones$end_date,
      pt_app_usage$usage_date,
      pt_diary$diary_date,
      pt_change_points$pre_segment_start_date,
      pt_change_points$change_point_date,
      pt_change_points$post_segment_end_date
    ),
    na.rm = TRUE
  )
  earliest_date <- min(
    c(
      pt_meds$start_date, pt_med_schedule_start_min, pt_seizures$event_date,
      pt_milestones$start_date, pt_app_usage$usage_date, pt_diary$diary_date,
      pt_change_points$pre_segment_start_date, pt_change_points$change_point_date
    ),
    na.rm = TRUE
  )
  if (!is.na(pt_app_activity_date) && nrow(pt_change_points) == 0) {
    earliest_date <- max(earliest_date, pt_app_activity_date, na.rm = TRUE)
  }

  if (!is.finite(latest_date)) {
    latest_date <- current_date
  }
  if (!is.finite(earliest_date)) {
    earliest_date <- if (!is.na(pt_app_activity_date)) pt_app_activity_date else latest_date - 30
  }

  plot_end_date <- min(latest_date, current_date)
  if (earliest_date > plot_end_date) {
    earliest_date <- plot_end_date - 30
  }
  pt_change_point_markers <- pt_change_points %>%
    transmute(
      marker_date = change_point_date,
      marker_direction = change_direction
    ) %>%
    filter(!is.na(marker_date), marker_date >= earliest_date, marker_date <= plot_end_date)

  pt_meds <- pt_meds %>%
    mutate(
      med_label = med_name,
      med_status = if_else(is.na(end_date), "Ongoing", "Ended"),
      end_plot_date = case_when(
        is.na(end_date) ~ plot_end_date,
        end_date > plot_end_date ~ plot_end_date,
        TRUE ~ end_date
      )
    ) %>%
    mutate(
      start_plot_date = pmax(start_date, earliest_date),
      end_plot_date = if_else(end_plot_date == start_plot_date,
                              pmin(end_plot_date + 1, plot_end_date),
                              end_plot_date),
      start_plot_date = if_else(end_plot_date == start_plot_date,
                                pmax(start_plot_date - 1, earliest_date),
                                start_plot_date)
    ) %>%
    filter(start_plot_date <= plot_end_date, end_plot_date >= earliest_date)

  pt_seizures <- pt_seizures %>%
    mutate(
      seizure_label = if_else(is.na(seizure_label) | seizure_label == "", "Unspecified", seizure_label)
    ) %>%
    filter(event_date >= earliest_date, event_date <= plot_end_date)
  pt_app_usage_daily <- tibble(usage_date = seq(earliest_date, plot_end_date, by = "day")) %>%
    left_join(
      pt_app_usage %>%
        group_by(usage_date) %>%
        summarise(
          app_actions = as.integer(any(!is.na(app_actions) & app_actions > 0L)),
          .groups = "drop"
        ),
      by = "usage_date"
    ) %>%
    mutate(
      app_actions = coalesce(as.integer(app_actions), 0L),
      app_usage_bin = as.character(app_actions),
      app_usage_label = "App usage"
    )
  pt_diary <- pt_diary %>%
    filter(diary_date >= earliest_date, diary_date <= plot_end_date) %>%
    mutate(
      diary_label = "Weekly Survey",
      diary_dot_color = case_when(
        diary_status == "All yes" ~ "#2E8B57",
        diary_status == "Any no" ~ "#C80813",
        TRUE ~ "#8A9197"
      )
    )
  n_months_plotted <- length(seq(
    floor_date(earliest_date, unit = "month"),
    floor_date(plot_end_date, unit = "month"),
    by = "1 month"
  ))
  if (n_months_plotted < 1) {
    n_months_plotted <- 1L
  }
  avg_seizures_per_month <- nrow(pt_seizures) / n_months_plotted
  avg_seizure_text <- sprintf(
    "Average seizures/month: %.2f",
    avg_seizures_per_month
  )

  y_levels <- character(0)
  if (nrow(pt_milestones) > 0) {
    milestone_levels <- pt_milestones %>%
      arrange(start_date) %>%
      pull(milestone_label) %>%
      unique()
    y_levels <- c(y_levels, milestone_levels)
  } else if (show_no_skills_label) {
    y_levels <- c(y_levels, "No skills")
  } else if (show_no_entered_skills_label) {
    y_levels <- c(y_levels, "No entered skills")
  }
  if (nrow(pt_app_usage_daily) > 0) {
    y_levels <- c(y_levels, "App usage")
  }
  if (nrow(pt_diary) > 0) {
    y_levels <- c(y_levels, "Weekly Survey")
  }
  if (nrow(pt_seizures) > 0) {
    seizure_levels <- pt_seizures %>% count(seizure_label, sort = TRUE) %>% pull(seizure_label)
    y_levels <- c(y_levels, seizure_levels)
  }
  if (nrow(pt_meds) > 0) {
    med_levels <- pt_meds %>% arrange(start_date) %>% pull(med_label) %>% unique()
    y_levels <- c(y_levels, rev(med_levels))
  }

  p <- ggplot()

  if (nrow(pt_meds) > 0) {
    p <- p +
      geom_segment(
        data = pt_meds,
        aes(x = start_plot_date, xend = end_plot_date, y = med_label, yend = med_label, color = med_status),
        linewidth = 2
      )
  }

  if (nrow(pt_app_usage_daily) > 0) {
    p <- p +
      geom_tile(
        data = pt_app_usage_daily,
        aes(x = usage_date, y = app_usage_label, fill = app_usage_bin),
        width = 0.95,
        height = 0.8
      )
  }

  if (nrow(pt_diary) > 0) {
    p <- p +
      geom_point(
        data = pt_diary,
        aes(x = diary_date, y = diary_label),
        color = pt_diary$diary_dot_color,
        size = 3,
        alpha = 0.9
      )
  }

  if (nrow(pt_seizures) > 0) {
    p <- p +
      geom_point(
        data = pt_seizures,
        aes(x = event_date, y = seizure_label),
        color = "#C80813",
        size = 2.5,
        alpha = 0.7
      )
  }

  if (nrow(pt_milestones) > 0) {
    p <- p +
      geom_segment(
        data = pt_milestones,
        aes(x = start_date, xend = end_date, y = milestone_label, yend = milestone_label),
        color = "#1A9993",
        linewidth = 2
      )
  }

  pt_increase_change_point_markers <- pt_change_point_markers %>%
    filter(marker_direction == "increase")
  pt_decrease_change_point_markers <- pt_change_point_markers %>%
    filter(marker_direction == "decrease")

  if (nrow(pt_increase_change_point_markers) > 0) {
    p <- p +
      geom_vline(
        data = pt_increase_change_point_markers,
        aes(xintercept = marker_date),
        color = "#C80813",
        linewidth = 0.8,
        alpha = 0.7
      )
  }

  if (nrow(pt_decrease_change_point_markers) > 0) {
    p <- p +
      geom_vline(
        data = pt_decrease_change_point_markers,
        aes(xintercept = marker_date),
        color = "#2B6CB0",
        linewidth = 0.8,
        alpha = 0.7
      )
  }

  p +
    scale_y_discrete(limits = y_levels) +
    scale_color_manual(values = c("Ongoing" = "#709AE1", "Ended" = "#8A9197"), drop = FALSE) +
    scale_fill_manual(
      values = c("0" = "#F2F4F7", "1" = "#3182BD"),
      drop = FALSE
    ) +
    scale_x_date(
      date_labels = "%b %Y",
      date_breaks = "1 month",
      minor_breaks = NULL,
      limits = c(earliest_date, plot_end_date),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = title_text,
      subtitle = avg_seizure_text,
      x = "Date",
      y = "",
      color = "Medication Status",
      fill = "App used/day"
    ) +
    theme_linedraw() +
    theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor.x = element_blank()
    )
}

for (pt in patients_with_data) {
  plt <- plot_patient_timeline(pt)
  if (!is.null(plt)) {
    pt_row_count <- patient_row_slots %>%
      filter(patient_id == pt) %>%
      pull(n_total_y_rows) %>%
      first()
    if (is.na(pt_row_count) || pt_row_count < 1) {
      pt_row_count <- 1
    }
    min_plot_height_frac <- 0.2
    plot_height_frac <- max(pt_row_count / max_y_rows, min_plot_height_frac)

    pt_filename <- patient_timeline_filenames %>%
      filter(patient_id == pt) %>%
      pull(timeline_filename) %>%
      first()
    if (is.na(pt_filename) || pt_filename == "") {
      pt_filename <- paste0(pt, "_timeline.pdf")
    }
    pt_output_pdf <- file.path("output/figs/timelines", pt_filename)
    grDevices::pdf(pt_output_pdf, width = 16, height = 10)
    suppressWarnings(
      print(
        plt,
        newpage = TRUE,
        vp = grid::viewport(x = 0.5, y = 0.5, width = 1, height = plot_height_frac)
      )
    )
    grDevices::dev.off()
  }
}
