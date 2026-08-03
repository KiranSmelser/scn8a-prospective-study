suppressPackageStartupMessages({
  library(cluster)
  library(dplyr)
  library(lubridate)
  library(readr)
})

source("src/analysis_config.R")

FEATURES_INPUT_PATH <- "output/tabs/clustering/seizure_freq_features.csv"
PANEL_INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
SEIZURES_INPUT_PATH <- "output/tabs/seizures/seizures.csv"
OUTPUT_DIR <- "output/tabs/clustering"
ASSIGNMENTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "cluster_assignments.csv")
SILHOUETTE_OUTPUT_PATH <- file.path(OUTPUT_DIR, "cluster_silhouette_scores.csv")
CLUSTER_COUNT <- 3L
PCA_COMPONENT_COUNT <- 2L
RATE_FEATURES <- c(
  "mean_monthly_seizure_rate",
  "iqr_monthly_seizure_rate"
)
DIRECT_SCALE_FEATURES <- c("proportion_zero_seizure_months")
REQUIRED_COLUMNS <- c("patient_id", "variant_p", RATE_FEATURES, DIRECT_SCALE_FEATURES)
PANEL_REQUIRED_COLUMNS <- c("patient_id", "variant_p", "month", "study_start_date", "study_end_date")
SEIZURE_REQUIRED_COLUMNS <- c("patient_id", "type", "date")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (input_path in c(FEATURES_INPUT_PATH, PANEL_INPUT_PATH, SEIZURES_INPUT_PATH)) {
  if (!file.exists(input_path)) {
    stop("Input file not found: ", input_path)
  }
}

robust_scale <- function(x) {
  center <- stats::median(x, na.rm = TRUE)
  spread <- stats::IQR(x, na.rm = TRUE)

  if (is.na(spread) || spread == 0) {
    spread <- stats::mad(x, center = center, constant = 1, na.rm = TRUE)
  }

  if (is.na(spread) || spread == 0) {
    spread <- 1
  }

  (x - center) / spread
}

parse_event_date <- function(x) {
  parsed <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC"))
  parsed <- dplyr::if_else(
    is.na(parsed),
    suppressWarnings(as.POSIXct(lubridate::ymd(x), tz = "UTC")),
    parsed
  )
  as.Date(parsed)
}

check_columns <- function(data, required_columns, data_label) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(data_label, " is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

build_monthly_rates <- function(patient_month_panel, monthly_counts) {
  patient_month_panel %>%
    left_join(monthly_counts, by = c("patient_id", "month")) %>%
    mutate(
      seizure_count = if_else(is.na(.data$seizure_count), 0, .data$seizure_count),
      next_month_start = as.Date(format(.data$month + 32, "%Y-%m-01")),
      month_end = .data$next_month_start - 1,
      observed_start = as.Date(pmax(.data$study_start_date, .data$month), origin = "1970-01-01"),
      observed_end = as.Date(pmin(.data$study_end_date, .data$month_end), origin = "1970-01-01"),
      observed_days_in_month = pmax(as.integer(.data$observed_end - .data$observed_start + 1), 0L)
    ) %>%
    filter(
      !is.na(.data$patient_id),
      !is.na(.data$month),
      !is.na(.data$study_start_date),
      !is.na(.data$study_end_date),
      !is.na(.data$seizure_count),
      .data$observed_days_in_month > 0
    ) %>%
    mutate(
      monthly_seizure_rate = .data$seizure_count / .data$observed_days_in_month * STANDARD_MONTH_DAYS,
      zero_seizure_month_flag = as.integer(.data$seizure_count == 0)
    )
}

build_feature_table <- function(patient_month_panel, monthly_counts) {
  build_monthly_rates(patient_month_panel, monthly_counts) %>%
    group_by(.data$patient_id) %>%
    summarise(
      variant_p = first_non_missing(.data$variant_p),
      mean_monthly_seizure_rate = mean(.data$monthly_seizure_rate),
      iqr_monthly_seizure_rate = stats::IQR(.data$monthly_seizure_rate),
      proportion_zero_seizure_months = mean(.data$zero_seizure_month_flag),
      .groups = "drop"
    )
}

build_type_feature_table <- function(patient_month_panel, seizure_events, seizure_type) {
  type_events <- seizure_events %>%
    filter(.data$type == seizure_type)

  eligible_patients <- type_events %>%
    distinct(.data$patient_id)

  monthly_counts <- type_events %>%
    mutate(month = as.Date(lubridate::floor_date(.data$event_date, "month"))) %>%
    count(.data$patient_id, .data$month, name = "seizure_count")

  patient_month_panel %>%
    semi_join(eligible_patients, by = "patient_id") %>%
    build_feature_table(monthly_counts)
}

scaled_feature_names <- c(
  paste0(RATE_FEATURES, "_log1p_scaled"),
  paste0(DIRECT_SCALE_FEATURES, "_scaled")
)

relabel_clusters_by_burden <- function(clustering_input, cluster_assignment) {
  cluster_order <- clustering_input %>%
    mutate(original_cluster = as.integer(cluster_assignment)) %>%
    group_by(.data$original_cluster) %>%
    summarise(
      severity_mean_rate = mean(.data$mean_monthly_seizure_rate),
      severity_iqr_rate = mean(.data$iqr_monthly_seizure_rate),
      severity_zero_months = mean(.data$proportion_zero_seizure_months),
      .groups = "drop"
    ) %>%
    arrange(
      desc(.data$severity_mean_rate),
      desc(.data$severity_iqr_rate),
      .data$severity_zero_months,
      .data$original_cluster
    ) %>%
    mutate(relabeled_cluster = row_number())

  cluster_map <- stats::setNames(cluster_order$relabeled_cluster, cluster_order$original_cluster)
  as.integer(cluster_map[as.character(cluster_assignment)])
}

cluster_with_pam <- function(feature_table, assignment_column, data_label) {
  check_columns(feature_table, REQUIRED_COLUMNS, data_label)

  clustering_input <- feature_table %>%
    transmute(
      patient_id = as.character(.data$patient_id),
      variant_p = as.character(.data$variant_p),
      across(all_of(RATE_FEATURES), ~ suppressWarnings(as.numeric(.x))),
      across(all_of(DIRECT_SCALE_FEATURES), ~ suppressWarnings(as.numeric(.x)))
    ) %>%
    filter(
      !is.na(.data$patient_id),
      if_all(all_of(c(RATE_FEATURES, DIRECT_SCALE_FEATURES)), ~ !is.na(.x))
    ) %>%
    mutate(
      across(all_of(RATE_FEATURES), log1p, .names = "{.col}_log1p"),
      across(all_of(c(paste0(RATE_FEATURES, "_log1p"), DIRECT_SCALE_FEATURES)), robust_scale, .names = "{.col}_scaled")
    )

  if (nrow(clustering_input) < CLUSTER_COUNT) {
    stop(
      "Not enough patients for requested cluster count in ",
      data_label,
      ". Patients available: ",
      nrow(clustering_input),
      "; requested k: ",
      CLUSTER_COUNT
    )
  }

  pca_input_matrix <- clustering_input %>%
    select(all_of(scaled_feature_names)) %>%
    as.matrix()

  if (any(!is.finite(pca_input_matrix))) {
    stop("Scaled feature matrix contains non-finite values for ", data_label, ".")
  }

  pca_fit <- stats::prcomp(pca_input_matrix, center = FALSE, scale. = FALSE)

  if (ncol(pca_fit$x) < PCA_COMPONENT_COUNT) {
    stop(
      "PCA produced fewer than ",
      PCA_COMPONENT_COUNT,
      " components for ",
      data_label,
      ". Components available: ",
      ncol(pca_fit$x)
    )
  }

  pc_feature_names <- paste0("pc", seq_len(PCA_COMPONENT_COUNT))
  pc_feature_matrix <- pca_fit$x[, seq_len(PCA_COMPONENT_COUNT), drop = FALSE]
  colnames(pc_feature_matrix) <- pc_feature_names

  distance_matrix <- stats::dist(pc_feature_matrix)
  pam_fit <- cluster::pam(pc_feature_matrix, k = CLUSTER_COUNT, metric = "euclidean", stand = FALSE)
  raw_cluster_assignment <- as.integer(pam_fit$clustering)
  relabeled_cluster_assignment <- relabel_clusters_by_burden(clustering_input, raw_cluster_assignment)
  silhouette_widths <- cluster::silhouette(raw_cluster_assignment, distance_matrix)

  assignments <- clustering_input %>%
    select("patient_id", "variant_p", all_of(RATE_FEATURES), all_of(DIRECT_SCALE_FEATURES)) %>%
    mutate(
      "{assignment_column}" := relabeled_cluster_assignment
    )

  list(
    assignments = assignments,
    silhouette = tibble::tibble(
      clustering_scope = data_label,
      assignment_column = assignment_column,
      n_patients = nrow(clustering_input),
      n_clusters = CLUSTER_COUNT,
      mean_silhouette_score = mean(silhouette_widths[, "sil_width"])
    )
  )
}

feature_table <- readr::read_csv(FEATURES_INPUT_PATH, show_col_types = FALSE)
patient_month_panel <- readr::read_csv(PANEL_INPUT_PATH, show_col_types = FALSE)
seizure_events <- readr::read_csv(SEIZURES_INPUT_PATH, show_col_types = FALSE)

check_columns(feature_table, REQUIRED_COLUMNS, "Feature input")
check_columns(patient_month_panel, PANEL_REQUIRED_COLUMNS, "Patient-month panel")
check_columns(seizure_events, SEIZURE_REQUIRED_COLUMNS, "Seizure input")

patient_month_panel <- patient_month_panel %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    variant_p = as.character(.data$variant_p),
    month = as.Date(.data$month),
    study_start_date = as.Date(.data$study_start_date),
    study_end_date = pmin(as.Date(.data$study_end_date), ANALYSIS_CUTOFF_DATE)
  ) %>%
  filter(
    !is.na(.data$patient_id),
    !is.na(.data$month),
    .data$month <= ANALYSIS_CUTOFF_DATE,
    !is.na(.data$study_start_date),
    !is.na(.data$study_end_date)
  )

seizure_events <- seizure_events %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    type = as.character(.data$type),
    event_date = parse_event_date(.data$date)
  ) %>%
  filter(
    !is.na(.data$patient_id),
    !is.na(.data$type),
    !is.na(.data$event_date)
  ) %>%
  inner_join(
    patient_month_panel %>%
      distinct(.data$patient_id, .data$study_start_date, .data$study_end_date),
    by = "patient_id"
  ) %>%
  filter(
    .data$event_date >= .data$study_start_date,
    .data$event_date <= .data$study_end_date
  )

all_type_clustering <- cluster_with_pam(feature_table, "pam_k3", "all seizure types")
tonic_clonic_clustering <- build_type_feature_table(patient_month_panel, seizure_events, "Tonic-clonic") %>%
  cluster_with_pam("pam_k3_tonic_clonic", "tonic-clonic seizures")
focal_clustering <- build_type_feature_table(patient_month_panel, seizure_events, "Focal") %>%
  cluster_with_pam("pam_k3_focal", "focal seizures")

cluster_assignments <- all_type_clustering$assignments %>%
  left_join(
    tonic_clonic_clustering$assignments %>% select("patient_id", pam_k3_tonic_clonic),
    by = "patient_id"
  ) %>%
  left_join(
    focal_clustering$assignments %>% select("patient_id", pam_k3_focal),
    by = "patient_id"
  ) %>%
  arrange(.data$patient_id)

silhouette_scores <- bind_rows(
  all_type_clustering$silhouette,
  tonic_clonic_clustering$silhouette,
  focal_clustering$silhouette
)

readr::write_csv(cluster_assignments, ASSIGNMENTS_OUTPUT_PATH)
readr::write_csv(silhouette_scores, SILHOUETTE_OUTPUT_PATH)
