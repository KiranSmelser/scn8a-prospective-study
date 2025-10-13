# Developmental milestone attainment bar chart

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(stringr)
})

milestones_raw <- readr::read_csv("data/prospective_development_milestones.csv", show_col_types = FALSE) %>%
  dplyr::mutate(
    milestone = as.character(.data$milestone),
    survey_instance_id = as.character(.data$survey_instance_id),
    status_numeric = suppressWarnings(as.numeric(.data$status)),
    achieved = .data$status_numeric == 1
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

milestones_summary <- milestones_raw %>%
  dplyr::filter(!is.na(.data$status_numeric)) %>%
  dplyr::group_by(.data$milestone) %>%
  dplyr::summarise(
    respondents = dplyr::n_distinct(.data$survey_instance_id),
    achieved_patients = dplyr::n_distinct(.data$survey_instance_id[.data$achieved]),
    pct_achieved = dplyr::if_else(respondents > 0, achieved_patients / respondents, 0),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    milestone_label = dplyr::recode(.data$milestone, !!!milestone_labels, .default = str_to_title(str_replace_all(.data$milestone, "_", " "))),
    pct_achieved = tidyr::replace_na(.data$pct_achieved, 0)
  ) %>%
  dplyr::arrange(dplyr::desc(.data$pct_achieved))

p_dev_attainment <- milestones_summary %>%
  dplyr::mutate(milestone_label = stats::reorder(.data$milestone_label, .data$pct_achieved)) %>%
  ggplot(aes(x = .data$milestone_label, y = .data$pct_achieved)) +
  geom_col(fill = "#5A9CA8") +
  coord_flip() +
  scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Developmental Milestone",
    y = "Patients",
    title = "Developmental Milestone Attainment",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank()
  )

ggsave(filename = "output/figs/dev_attainment.png", plot = p_dev_attainment, width = 10, height = 7, dpi = 600)