# Create timeline plots for Helpilepsy patients

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

dir.create("output/figs", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tabs", recursive = TRUE, showWarnings = FALSE)
current_date <- Sys.Date()

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
medications <- readr::read_csv("data/medications.csv", show_col_types = FALSE)
med_dosages <- readr::read_csv("data/med_dosages.csv", show_col_types = FALSE)
med_intakes <- readr::read_csv("data/med_intakes.csv", show_col_types = FALSE)
surveys <- readr::read_csv("data/prospective_surveys.csv", show_col_types = FALSE)
milestones_raw <- readr::read_csv("data/prospective_development_milestones.csv", show_col_types = FALSE)

# Seizure events

seizure_events <- events %>%
  filter(tolower(type) == "seizure") %>%
  mutate(
    event_date = parse_event_date(date),
    seizure_type = if_else(is.na(seizure_type) | seizure_type == "", "Unspecified", seizure_type),
    seizure_type = str_replace_all(seizure_type, "_", " "),
    seizure_type = str_squish(seizure_type),
    seizure_label = str_to_title(seizure_type)
  ) %>%
  filter(!is.na(event_date)) %>%
  select(patient_id, seizure_type, seizure_label, event_date)

# Medication intervals

med_intervals_raw <- bind_rows(
  med_dosages %>%
    transmute(patient_id, medication_id, start_date = ymd(from), end_date = ymd(to)),
  med_intakes %>%
    transmute(patient_id, medication_id, start_date = ymd(intake_from), end_date = ymd(intake_to))
) %>%
  mutate(end_date = if_else(!is.na(end_date) & end_date < start_date, as.Date(NA), end_date)) %>%
  filter(!is.na(start_date)) %>%
  distinct()

med_schedule_summary_raw <- med_intervals_raw %>%
  group_by(patient_id) %>%
  summarise(
    med_schedule_start_min = min(start_date, na.rm = TRUE),
    med_schedule_start_max = max(start_date, na.rm = TRUE),
    med_schedule_end_max = if (all(is.na(end_date))) as.Date(NA) else max(end_date, na.rm = TRUE),
    .groups = "drop"
  )

med_intervals <- med_intervals_raw %>%
  left_join(medications %>% select(patient_id, medication_id, name), by = c("patient_id", "medication_id")) %>%
  mutate(
    med_name = if_else(is.na(name) | name == "", medication_id, name),
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
      status_flag = (status_numeric == 1),
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

# Patient names

whatsapp_names <- readr::read_csv("data/whatsapp_status.csv", col_names = c(
  "patient_id", "wh_first_name", "wh_last_name", "status", "source", "run_timestamp"
), show_col_types = FALSE) %>%
  mutate(
    whatsapp_name = str_squish(str_trim(paste(wh_first_name, wh_last_name)))
  ) %>%
  select(patient_id, whatsapp_name) %>%
  filter(whatsapp_name != "")

patients_with_data <- union(
  union(unique(seizure_events$patient_id), unique(med_intervals$patient_id)),
  unique(milestone_intervals$patient_id)
)
patients_with_data <- patients_with_data[!is.na(patients_with_data)]

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
  left_join(whatsapp_names, by = "patient_id") %>%
  mutate(
    display_name = coalesce(whatsapp_name, patient_name),
    display_name = if_else(is.na(display_name) | display_name == "", patient_id, display_name)
  ) %>%
  select(patient_id, display_name)

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

medications_missing_schedule <- medications %>%
  anti_join(med_intervals_raw %>% distinct(patient_id, medication_id), by = c("patient_id", "medication_id")) %>%
  select(patient_id, medication_id, name, reason, treatment_type, intake_type, createdAt, updatedAt) %>%
  arrange(patient_id, name)

medications_missing_schedule_summary <- medications_missing_schedule %>%
  group_by(patient_id) %>%
  summarise(n_med_records_missing_schedule = n_distinct(medication_id), .groups = "drop")

timeline_qc <- patient_order %>%
  left_join(patient_names, by = "patient_id") %>%
  left_join(seizure_summary, by = "patient_id") %>%
  left_join(med_summary, by = "patient_id") %>%
  left_join(med_schedule_summary_raw, by = "patient_id") %>%
  left_join(milestone_summary, by = "patient_id") %>%
  left_join(medications_missing_schedule_summary, by = "patient_id") %>%
  mutate(
    n_seizure_events = replace_na(n_seizure_events, 0L),
    n_seizure_types = replace_na(n_seizure_types, 0L),
    n_med_intervals = replace_na(n_med_intervals, 0L),
    n_medications = replace_na(n_medications, 0L),
    n_medications_ongoing = replace_na(n_medications_ongoing, 0L),
    n_milestones_achieved = replace_na(n_milestones_achieved, 0L),
    n_med_records_missing_schedule = replace_na(n_med_records_missing_schedule, 0L)
  ) %>%
  rowwise() %>%
  mutate(
    timeline_start_date = {
      dates <- c(seizure_first_date, med_start_min, med_schedule_start_min, milestone_first_date)
      dates <- dates[!is.na(dates)]
      if (length(dates) == 0) as.Date(NA) else min(dates)
    },
    timeline_end_date = {
      dates <- c(seizure_last_date, med_end_max, med_start_max, med_schedule_start_max, med_schedule_end_max, milestone_last_date)
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

# Plotting

plot_patient_timeline <- function(pt_id) {
  pt_info <- patients %>% filter(patient_id == pt_id) %>% slice_head(n = 1)
  pt_name <- patient_names %>% filter(patient_id == pt_id) %>% pull(display_name) %>% first()
  if (is.na(pt_name) || pt_name == "") {
    pt_name <- str_trim(paste(pt_info$first_name, pt_info$last_name))
  }
  title_text <- ifelse(is.na(pt_name) || pt_name == "", pt_id, paste0(pt_name, " (", pt_id, ")"))

  pt_meds <- med_intervals %>% filter(patient_id == pt_id)
  pt_med_schedule <- med_schedule_summary_raw %>% filter(patient_id == pt_id) %>% slice_head(n = 1)
  pt_seizures <- seizure_events %>% filter(patient_id == pt_id)
  pt_milestones <- milestone_intervals %>% filter(patient_id == pt_id)

  if (nrow(pt_meds) == 0 && nrow(pt_seizures) == 0 && nrow(pt_milestones) == 0) {
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
      pt_milestones$start_date, pt_milestones$end_date
    ),
    na.rm = TRUE
  )
  earliest_date <- min(
    c(pt_meds$start_date, pt_med_schedule_start_min, pt_seizures$event_date, pt_milestones$start_date),
    na.rm = TRUE
  )

  if (!is.finite(latest_date)) {
    latest_date <- current_date
  }
  if (!is.finite(earliest_date)) {
    earliest_date <- latest_date - 30
  }

  plot_end_date <- min(latest_date, current_date)
  if (earliest_date > plot_end_date) {
    earliest_date <- plot_end_date - 30
  }

  pt_meds <- pt_meds %>%
    mutate(
      med_label = med_name,
      med_status = if_else(is.na(end_date), "Ongoing", "Ended"),
      end_plot_date = case_when(
        is.na(end_date) ~ plot_end_date,
        end_date > plot_end_date ~ plot_end_date,
        TRUE ~ end_date
      )
    )

  pt_seizures <- pt_seizures %>%
    mutate(seizure_label = if_else(is.na(seizure_label) | seizure_label == "", "Unspecified", seizure_label))

  y_levels <- character(0)
  if (nrow(pt_milestones) > 0) {
    milestone_levels <- pt_milestones %>%
      arrange(start_date) %>%
      pull(milestone_label) %>%
      unique()
    y_levels <- c(y_levels, milestone_levels)
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
        aes(x = start_date, xend = end_plot_date, y = med_label, yend = med_label, color = med_status),
        linewidth = 2
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

  p +
    scale_y_discrete(limits = y_levels) +
    scale_color_manual(values = c("Ongoing" = "#709AE1", "Ended" = "#8A9197"), drop = FALSE) +
    scale_x_date(
      date_labels = "%b %Y",
      limits = c(earliest_date, plot_end_date),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = title_text,
      x = "Date",
      y = "",
      color = "Medication Status"
    ) +
    theme_linedraw() +
    theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

pdf("output/figs/helpilepsy_timelines.pdf", width = 16, height = 10)
for (pt in patients_with_data) {
  plt <- plot_patient_timeline(pt)
  if (!is.null(plt)) {
    suppressWarnings(print(plt))
  }
}
dev.off()
