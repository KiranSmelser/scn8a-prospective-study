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
  "initial_tonic",
  "seizure_type_focal_aware_simple_partial_seizure",
  "seizure_type_typical_absence_petit_mal",
  "dev_skill_brush_teeth_with_no_help"
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

registry <- readr::read_csv(REGISTRY_INPUT_PATH, show_col_types = FALSE)
cluster_assignments <- readr::read_csv(CLUSTER_ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE)

check_columns(registry, c("patient_id", TARGETED_REGISTRY_SOURCE_COLUMNS), "Registry input")
check_columns(cluster_assignments, c("patient_id", "pam_k3"), "Cluster assignments")

targeted_registry_features <- registry %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    registry_brush_teeth_independently = suppressWarnings(as.numeric(.data$dev_skill_brush_teeth_with_no_help)),
    registry_cluster2_signature = row_max_binary(
      .data$initial_tonic,
      .data$seizure_type_focal_aware_simple_partial_seizure,
      .data$seizure_type_typical_absence_petit_mal
    )
  ) %>%
  filter(!is.na(.data$patient_id)) %>%
  mutate(
    across(starts_with("registry_"), ~ tidyr::replace_na(suppressWarnings(as.numeric(.x)), 0))
  ) %>%
  distinct(.data$patient_id, .keep_all = TRUE)

feature_definitions <- tibble::tribble(
  ~feature, ~feature_type,
  "registry_brush_teeth_independently", "binary",
  "registry_cluster2_signature", "binary"
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
