# Developmental milestone attainment

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(stringr)
})

source("src/analysis_config.R")

CLUSTER_ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"

if (!file.exists(CLUSTER_ASSIGNMENTS_INPUT_PATH)) {
  stop(
    "Cluster assignments not found: ",
    CLUSTER_ASSIGNMENTS_INPUT_PATH,
    ". Run src/clustering/clustering.R before src/desc/dev_attainment.R."
  )
}

survey_lookup_on_or_before_cutoff <- readr::read_csv("data/prospective_surveys.csv", show_col_types = FALSE) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    survey_instance_id = as.character(.data$survey_instance_id),
    survey_date = as.Date(suppressWarnings(lubridate::ymd_hms(.data$prospective_study_timestamp_utc, quiet = TRUE, tz = "UTC"))),
    survey_date = dplyr::coalesce(
      .data$survey_date,
      as.Date(suppressWarnings(lubridate::mdy_hm(.data$prospective_study_timestamp, quiet = TRUE, tz = "UTC")))
    )
  ) %>%
  dplyr::filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    !is.na(.data$survey_date),
    .data$survey_date <= ANALYSIS_CUTOFF_DATE
  ) %>%
  dplyr::distinct(.data$survey_instance_id, .data$patient_id)

cluster_assignments <- readr::read_csv(CLUSTER_ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    pam_cluster = as.character(.data$pam_k3)
  )

milestones_raw <- readr::read_csv("data/prospective_development_milestones.csv", show_col_types = FALSE) %>%
  dplyr::mutate(survey_instance_id = as.character(.data$survey_instance_id)) %>%
  dplyr::inner_join(survey_lookup_on_or_before_cutoff, by = "survey_instance_id") %>%
  dplyr::inner_join(cluster_assignments, by = "patient_id") %>%
  dplyr::mutate(
    milestone = as.character(.data$milestone),
    status_numeric = suppressWarnings(as.numeric(.data$status)),
    achieved = .data$status_numeric %in% c(2, 3, 5)
  )

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

patient_milestone_attainment <- milestones_raw %>%
  dplyr::filter(!is.na(.data$status_numeric)) %>%
  dplyr::group_by(.data$patient_id, .data$pam_cluster, .data$milestone) %>%
  dplyr::summarise(
    achieved = any(.data$achieved, na.rm = TRUE),
    .groups = "drop"
  )

milestones_summary <- patient_milestone_attainment %>%
  dplyr::group_by(.data$milestone) %>%
  dplyr::mutate(patients_with_status = dplyr::n_distinct(.data$patient_id)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(.data$milestone, .data$pam_cluster) %>%
  dplyr::summarise(
    patients_with_status = dplyr::first(.data$patients_with_status),
    achieved_patients = dplyr::n_distinct(.data$patient_id[.data$achieved]),
    pct_achieved = dplyr::if_else(.data$patients_with_status > 0, .data$achieved_patients / .data$patients_with_status, 0),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    milestone_label = dplyr::recode(.data$milestone, !!!milestone_labels, .default = str_to_title(str_replace_all(.data$milestone, "_", " "))),
    pam_cluster = factor(.data$pam_cluster, levels = sort(unique(.data$pam_cluster))),
    pct_achieved = tidyr::replace_na(.data$pct_achieved, 0)
  ) %>%
  dplyr::group_by(.data$milestone, .data$milestone_label) %>%
  dplyr::mutate(total_pct_achieved = sum(.data$pct_achieved, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(.data$total_pct_achieved), .data$pam_cluster)

p_dev_attainment <- milestones_summary %>%
  dplyr::mutate(
    milestone_label = stats::reorder(.data$milestone_label, .data$total_pct_achieved),
    pct_label = dplyr::if_else(.data$pct_achieved > 0, scales::percent(.data$pct_achieved, accuracy = 1), NA_character_)
  ) %>%
  ggplot(aes(x = .data$milestone_label, y = .data$pct_achieved, fill = .data$pam_cluster)) +
  geom_col(color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = .data$pct_label),
    position = position_stack(vjust = 0.5),
    size = 2.8,
    color = "#1F2933",
    na.rm = TRUE
  ) +
  coord_flip() +
  scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_brewer(palette = "Set2", name = "Cluster") +
  labs(
    x = "Developmental Milestone",
    y = "Patients",
    title = "Developmental Milestone Attainment",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "bottom"
  )

png_output_path <- "output/figs/dev_attainment.png"
if (file.exists(png_output_path)) {
  invisible(file.remove(png_output_path))
}

ggsave(filename = "output/figs/dev_attainment.pdf", plot = p_dev_attainment, width = 10, height = 7, bg = "white")
