# Summarize registry phenotype features included in the combined model by cluster.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

REGISTRY_INPUT_PATH <- "data/registry.csv"
CLUSTER_ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
OUTPUT_DIR <- "output/tabs/modeling/prediction/combined"
OUTPUT_PATH <- file.path(OUTPUT_DIR, "combined_feature_table.csv")

TARGET_CLUSTERS <- c("1", "2", "3")

TARGETED_REGISTRY_SOURCE_COLUMNS <- c(
  "age_seizure_onset_months",
  "seizure_type_tonic_clonic_grand_mal",
  "seizure_type_focal_aware_simple_partial_seizure",
  "seizure_type_focal_impaired_awareness_complex_partial_seizure_limbic_psychomotor",
  "dev_skill_sit_unsupported",
  "dev_skill_brush_teeth_with_no_help",
  "dev_skill_name_colors",
  "dev_skill_wash_and_dry_hands",
  "dev_skill_used_a_2_word_combination",
  "dev_skill_spoken_in_phrases",
  "dev_skill_read"
)

required_files <- c(REGISTRY_INPUT_PATH, CLUSTER_ASSIGNMENTS_INPUT_PATH)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Required input file(s) not found: ", paste(missing_files, collapse = ", "), call. = FALSE)
}

check_columns <- function(data, required_columns, label) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

row_max_binary <- function(...) {
  values <- cbind(...)
  values <- apply(values, 2, function(x) suppressWarnings(as.numeric(x)))
  row_value <- apply(
    values,
    1,
    function(x) {
      if (all(is.na(x))) {
        return(0)
      }
      max(x, na.rm = TRUE)
    }
  )
  as.numeric(row_value > 0)
}

row_mean_score <- function(...) {
  values <- cbind(...)
  values <- apply(values, 2, function(x) suppressWarnings(as.numeric(x)))
  row_value <- rowMeans(values, na.rm = TRUE)
  row_value[is.nan(row_value)] <- 0
  row_value
}

registry <- readr::read_csv(REGISTRY_INPUT_PATH, show_col_types = FALSE)
cluster_assignments <- readr::read_csv(CLUSTER_ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE)

check_columns(registry, c("patient_id", TARGETED_REGISTRY_SOURCE_COLUMNS), "Registry input")
check_columns(cluster_assignments, c("patient_id", "pam_k3"), "Cluster assignments")

onset_months <- suppressWarnings(as.numeric(registry$age_seizure_onset_months))
onset_fill_value <- stats::median(onset_months[is.finite(onset_months)], na.rm = TRUE)
if (!is.finite(onset_fill_value)) {
  onset_fill_value <- 0
}
onset_months_imputed <- onset_months
onset_months_imputed[!is.finite(onset_months_imputed)] <- onset_fill_value

targeted_registry_features <- registry %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    registry_age_seizure_onset_source_present = as.numeric(is.finite(onset_months)),
    registry_tonic_clonic_history = suppressWarnings(as.numeric(.data$seizure_type_tonic_clonic_grand_mal)),
    registry_focal_aware_or_impaired = row_max_binary(
      .data$seizure_type_focal_aware_simple_partial_seizure,
      .data$seizure_type_focal_impaired_awareness_complex_partial_seizure_limbic_psychomotor
    ),
    registry_log1p_age_seizure_onset_months = log1p(onset_months_imputed),
    registry_higher_dev_language_adl_score = row_mean_score(
      .data$dev_skill_brush_teeth_with_no_help,
      .data$dev_skill_name_colors,
      .data$dev_skill_wash_and_dry_hands,
      .data$dev_skill_used_a_2_word_combination,
      .data$dev_skill_spoken_in_phrases,
      .data$dev_skill_read
    ),
    registry_dev_sit_unsupported = suppressWarnings(as.numeric(.data$dev_skill_sit_unsupported)),
    registry_dev_name_colors = suppressWarnings(as.numeric(.data$dev_skill_name_colors))
  ) %>%
  filter(!is.na(.data$patient_id)) %>%
  mutate(
    across(starts_with("registry_"), ~ tidyr::replace_na(suppressWarnings(as.numeric(.x)), 0))
  ) %>%
  distinct(.data$patient_id, .keep_all = TRUE)

feature_definitions <- tibble::tribble(
  ~feature, ~feature_type,
  "registry_tonic_clonic_history", "binary",
  "registry_focal_aware_or_impaired", "binary",
  "registry_log1p_age_seizure_onset_months", "continuous",
  "registry_higher_dev_language_adl_score", "continuous",
  "registry_dev_sit_unsupported", "binary",
  "registry_dev_name_colors", "binary"
)

analysis_data <- targeted_registry_features %>%
  inner_join(
    cluster_assignments %>%
      transmute(patient_id = as.character(.data$patient_id), pam_k3 = as.character(.data$pam_k3)),
    by = "patient_id"
  ) %>%
  filter(.data$pam_k3 %in% TARGET_CLUSTERS)

feature_counts <- analysis_data %>%
  pivot_longer(
    cols = all_of(feature_definitions$feature),
    names_to = "feature",
    values_to = "value"
  ) %>%
  left_join(feature_definitions, by = "feature") %>%
  group_by(.data$feature, .data$feature_type, .data$pam_k3) %>%
  summarise(
    summary_value = if_else(
      first(.data$feature_type) == "continuous",
      round(mean(.data$value, na.rm = TRUE), 3),
      sum(.data$value > 0, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    cluster = paste0("cluster_", .data$pam_k3)
  ) %>%
  select(
    "feature",
    "cluster",
    "summary_value"
  ) %>%
  pivot_wider(
    names_from = "cluster",
    values_from = "summary_value"
  ) %>%
  rename(
    "cluster 1" = "cluster_1",
    "cluster 2" = "cluster_2",
    "cluster 3" = "cluster_3"
  ) %>%
  arrange(.data$feature)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(feature_counts, OUTPUT_PATH)
