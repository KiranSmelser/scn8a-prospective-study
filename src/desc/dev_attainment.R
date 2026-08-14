# Developmental milestone attainment

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
  library(stringr)
})

source("src/data_corrections.R")

CLUSTER_ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
REGISTRY_INPUT_PATH <- "data/registry.csv"

if (!file.exists(CLUSTER_ASSIGNMENTS_INPUT_PATH)) {
  stop(
    "Cluster assignments not found: ",
    CLUSTER_ASSIGNMENTS_INPUT_PATH,
    ". Run src/clustering/clustering.R before src/desc/dev_attainment.R."
  )
}

corrected_surveys <- read_prospective_surveys_corrected()

survey_lookup_on_or_before_cutoff <- corrected_surveys %>%
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

registry_dev_map <- tibble::tribble(
  ~registry_column, ~milestone,
  "dev_skill_try_to_follow_things_with_eyes_and_recognize_people_at_a_distance_visual_tracking", "eye_ps",
  "dev_skill_grasp_an_object", "grasp_ps",
  "dev_skill_try_to_get_things_that_are_out_of_reach", "reach_ps",
  "dev_skill_grasp_an_object_e_g_picks_up_things_like_cereal_o_s_between_thumb_and_index_finger", "pincer_ps",
  "dev_skill_build_towers_of_more_than_6_blocks", "blocks_ps",
  "dev_skill_draw_a_circle", "circle_ps",
  "dev_skill_hold_their_head_steady_unsupported_up_to_90_degrees", "hc_ps",
  "dev_skill_roll_over_in_both_directions_front_to_back_back_to_front", "roll_ps",
  "dev_skill_sit_unsupported", "sit_ps",
  "dev_skill_stand_with_support", "stand_ps",
  "dev_skill_walk", "walk_ps",
  "dev_skill_run", "run_ps",
  "dev_skill_smile_spontaneously_especially_at_people", "smile_ps",
  "dev_skill_wave_bye_bye", "wave_ps",
  "dev_skill_drink_from_a_cup", "cup_ps",
  "dev_skill_use_a_spoon_fork", "fork_ps",
  "dev_skill_wash_and_dry_hands", "wh_ps",
  "dev_skill_brush_teeth_with_no_help", "bt_ps",
  "dev_skill_vocalized_sounds_like_cooing_and_gurgling", "vocalize_ps",
  "dev_skill_laughed", "laugh_ps",
  "dev_skill_babbled_e_g_bababa_dadada", "babble_ps",
  "dev_skill_used_a_2_word_combination", "words_ps",
  "dev_skill_spoken_in_phrases", "phrase_ps",
  "dev_skill_name_colors", "namecolors_ps",
  "dev_skill_read", "reade_ps"
)

registry <- readr::read_csv(REGISTRY_INPUT_PATH, show_col_types = FALSE)
missing_registry_columns <- setdiff(
  c("patient_id", registry_dev_map$registry_column),
  names(registry)
)
if (length(missing_registry_columns) > 0) {
  stop(
    "Registry is missing required development column(s): ",
    paste(missing_registry_columns, collapse = ", ")
  )
}

prospective_milestones_raw <- read_prospective_child_records_corrected(
  "data/prospective_development_milestones.csv",
  corrected_surveys
) %>%
  dplyr::mutate(survey_instance_id = as.character(.data$survey_instance_id)) %>%
  dplyr::inner_join(survey_lookup_on_or_before_cutoff, by = "survey_instance_id") %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    milestone = as.character(.data$milestone),
    status_numeric = suppressWarnings(as.numeric(.data$status)),
    achieved = .data$status_numeric %in% c(2, 3, 5),
    data_source = "Prospective survey"
  )

registry_milestones_raw <- registry %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    dplyr::across(dplyr::all_of(registry_dev_map$registry_column))
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(registry_dev_map$registry_column),
    names_to = "registry_column",
    values_to = "status_numeric"
  ) %>%
  dplyr::inner_join(registry_dev_map, by = "registry_column") %>%
  dplyr::mutate(
    status_numeric = suppressWarnings(as.numeric(.data$status_numeric)),
    achieved = .data$status_numeric == 1,
    data_source = "Registry"
  ) %>%
  dplyr::filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    !is.na(.data$status_numeric)
  ) %>%
  dplyr::select("patient_id", "milestone", "status_numeric", "achieved", "data_source")

unexpected_registry_statuses <- registry_milestones_raw %>%
  dplyr::filter(!.data$status_numeric %in% c(0, 1)) %>%
  dplyr::distinct(.data$status_numeric) %>%
  dplyr::pull(.data$status_numeric)
if (length(unexpected_registry_statuses) > 0) {
  stop(
    "Registry development fields contain value(s) other than 0, 1, or missing: ",
    paste(unexpected_registry_statuses, collapse = ", ")
  )
}

milestones_raw <- dplyr::bind_rows(
  prospective_milestones_raw,
  registry_milestones_raw
) %>%
  dplyr::inner_join(cluster_assignments, by = "patient_id")

dev_category_order <- c("Fine Motor", "Gross Motor", "Social", "Language")

dev_domain_map <- tibble::tribble(
  ~milestone, ~dev_category, ~dev_skill,
  "eye_ps", "Fine Motor", "Visual tracking",
  "grasp_ps", "Fine Motor", "Grasping",
  "reach_ps", "Fine Motor", "Reaching",
  "pincer_ps", "Fine Motor", "Pincer grasp",
  "blocks_ps", "Fine Motor", "Stacking blocks",
  "circle_ps", "Fine Motor", "Drawing circles",
  "hc_ps", "Gross Motor", "Head control",
  "roll_ps", "Gross Motor", "Rolling over",
  "sit_ps", "Gross Motor", "Sitting independently",
  "stand_ps", "Gross Motor", "Standing independently",
  "walk_ps", "Gross Motor", "Walking independently",
  "run_ps", "Gross Motor", "Running",
  "smile_ps", "Social", "Social smiling",
  "wave_ps", "Social", "Waving",
  "cup_ps", "Social", "Using cup",
  "fork_ps", "Social", "Using spoon",
  "wh_ps", "Social", "Washing hands",
  "bt_ps", "Social", "Brushing teeth",
  "vocalize_ps", "Language", "Vocalizing",
  "laugh_ps", "Language", "Laughing",
  "babble_ps", "Language", "Babbling",
  "words_ps", "Language", "Two-word phrases",
  "phrase_ps", "Language", "Full phrases",
  "namecolors_ps", "Language", "Naming colors",
  "reade_ps", "Language", "Reading"
) %>%
  dplyr::mutate(
    dev_category = factor(.data$dev_category, levels = dev_category_order),
    dev_skill_order = dplyr::row_number()
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
  dplyr::left_join(dev_domain_map, by = "milestone") %>%
  dplyr::mutate(
    dev_category = tidyr::replace_na(.data$dev_category, "Unmapped"),
    dev_category = factor(.data$dev_category, levels = c(dev_category_order, "Unmapped")),
    dev_skill = dplyr::coalesce(.data$dev_skill, str_to_title(str_replace_all(.data$milestone, "_", " "))),
    dev_skill_order = tidyr::replace_na(.data$dev_skill_order, 999L),
    pam_cluster = factor(.data$pam_cluster, levels = sort(unique(.data$pam_cluster))),
    pct_achieved = tidyr::replace_na(.data$pct_achieved, 0)
  ) %>%
  dplyr::group_by(.data$milestone, .data$dev_category, .data$dev_skill, .data$dev_skill_order) %>%
  dplyr::mutate(total_pct_achieved = sum(.data$pct_achieved, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(.data$dev_category, .data$dev_skill_order, .data$pam_cluster)

dev_skill_levels <- milestones_summary %>%
  dplyr::distinct(.data$dev_category, .data$dev_skill_order, .data$dev_skill) %>%
  dplyr::arrange(.data$dev_category, .data$dev_skill_order, .data$dev_skill) %>%
  dplyr::pull(.data$dev_skill)

p_dev_attainment <- milestones_summary %>%
  dplyr::mutate(
    dev_skill = factor(.data$dev_skill, levels = rev(dev_skill_levels)),
    pct_label = dplyr::if_else(.data$pct_achieved > 0, scales::percent(.data$pct_achieved, accuracy = 1), NA_character_)
  ) %>%
  ggplot(aes(x = .data$pct_achieved, y = .data$dev_skill, fill = .data$pam_cluster)) +
  geom_col(color = "white", linewidth = 0.2) +
  geom_text(
    aes(
      label = .data$pct_label,
      color = dplyr::if_else(as.character(.data$pam_cluster) == "3", "white", "#1F2933")
    ),
    position = position_stack(vjust = 0.5),
    size = 2.8,
    na.rm = TRUE
  ) +
  scale_color_identity() +
  scale_x_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(
    values = EPILEPSIA_CLUSTER_COLORS,
    breaks = names(EPILEPSIA_CLUSTER_COLORS),
    name = "Cluster"
  ) +
  facet_grid(
    rows = vars(dev_category),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  labs(
    x = "Patients",
    y = NULL,
    title = "Developmental Milestone Attainment",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold")
  )

png_output_path <- "output/figs/dev_attainment.png"
if (file.exists(png_output_path)) {
  invisible(file.remove(png_output_path))
}

ggsave(filename = "output/figs/dev_attainment.pdf", plot = p_dev_attainment, width = 10, height = 7, bg = "white")
