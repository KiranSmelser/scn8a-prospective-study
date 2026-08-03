suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

source("src/analysis_config.R")

FEATURES_INPUT_PATH <- "output/tabs/clustering/seizure_freq_features.csv"
ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
SEGMENTS_INPUT_PATH <- "output/tabs/changepoints/change_point_segments.csv"
CHANGEPOINTS_INPUT_PATH <- "output/tabs/changepoints/patient_change_points.csv"
OUTPUT_DIR <- "output/tabs/changepoints"
OUTPUT_PATH <- file.path(OUTPUT_DIR, "cluster_fluctuation.csv")

FEATURE_REQUIRED_COLUMNS <- c(
  "patient_id",
  "mean_monthly_seizure_rate",
  "iqr_monthly_seizure_rate",
  "proportion_zero_seizure_months"
)
ASSIGNMENT_REQUIRED_COLUMNS <- c("patient_id", "pam_k3")
SEGMENT_REQUIRED_COLUMNS <- c("patient_id", "observed_days", "seizure_count")
CHANGEPOINT_REQUIRED_COLUMNS <- c(
  "patient_id",
  "candidate_week",
  "significant",
  "pre_observed_days",
  "post_observed_days",
  "pre_seizure_count",
  "post_seizure_count"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

check_input <- function(path, required_columns, label) {
  if (!file.exists(path)) {
    stop(label, " input file not found: ", path, call. = FALSE)
  }

  data <- readr::read_csv(path, show_col_types = FALSE)
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      label,
      " input is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  data
}

median_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }

  stats::median(x)
}

iqr_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }

  stats::IQR(x)
}

median_available <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }

  stats::median(x)
}

format_interpretation <- function(burden_rank, n_clusters) {
  if (n_clusters == 1L) {
    return("single observed cluster")
  }

  if (burden_rank == 1L) {
    return("high burden / high fluctuation")
  }

  if (burden_rank == n_clusters) {
    return("low burden / intermittent")
  }

  "moderate burden / moderate fluctuation"
}

features <- check_input(FEATURES_INPUT_PATH, FEATURE_REQUIRED_COLUMNS, "Seizure-frequency feature")
assignments <- check_input(ASSIGNMENTS_INPUT_PATH, ASSIGNMENT_REQUIRED_COLUMNS, "Cluster assignment")
segments <- check_input(SEGMENTS_INPUT_PATH, SEGMENT_REQUIRED_COLUMNS, "Changepoint segment")
change_points <- check_input(CHANGEPOINTS_INPUT_PATH, CHANGEPOINT_REQUIRED_COLUMNS, "Patient changepoint")

cluster_lookup <- assignments %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    pam_k3 = suppressWarnings(as.integer(.data$pam_k3))
  ) %>%
  filter(!is.na(.data$patient_id), !is.na(.data$pam_k3)) %>%
  distinct(.data$patient_id, .keep_all = TRUE)

monthly_summary <- features %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    median_monthly_rate_per_28d = suppressWarnings(as.numeric(.data$mean_monthly_seizure_rate)),
    patient_monthly_iqr_per_28d = suppressWarnings(as.numeric(.data$iqr_monthly_seizure_rate)),
    zero_month_proportion = suppressWarnings(as.numeric(.data$proportion_zero_seizure_months))
  ) %>%
  inner_join(cluster_lookup, by = "patient_id") %>%
  filter(
    if_all(
      all_of(c(
        "median_monthly_rate_per_28d",
        "patient_monthly_iqr_per_28d",
        "zero_month_proportion"
      )),
      is.finite
    )
  ) %>%
  group_by(.data$pam_k3) %>%
  summarise(
    n_patients = n_distinct(.data$patient_id),
    median_monthly_rate_per_28d = median_finite(.data$median_monthly_rate_per_28d),
    median_patient_monthly_iqr_per_28d = median_finite(.data$patient_monthly_iqr_per_28d),
    median_zero_month_proportion = median_finite(.data$zero_month_proportion),
    .groups = "drop"
  )

segment_summary <- segments %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    observed_days = suppressWarnings(as.numeric(.data$observed_days)),
    seizure_count = suppressWarnings(as.numeric(.data$seizure_count))
  ) %>%
  inner_join(cluster_lookup, by = "patient_id") %>%
  mutate(
    segment_rate_per_28d = if_else(
      is.finite(.data$observed_days) & .data$observed_days > 0,
      STANDARD_MONTH_DAYS * .data$seizure_count / .data$observed_days,
      NA_real_
    )
  ) %>%
  group_by(.data$pam_k3) %>%
  summarise(
    n_changepoint_analyzable_patients = n_distinct(.data$patient_id),
    n_changepoint_segments = n(),
    median_segment_rate_per_28d = median_finite(.data$segment_rate_per_28d),
    segment_rate_iqr_per_28d = iqr_finite(.data$segment_rate_per_28d),
    .groups = "drop"
  )

changepoint_summary <- change_points %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    candidate_week = as.Date(.data$candidate_week),
    significant = .data$significant %in% TRUE | as.character(.data$significant) == "TRUE",
    pre_observed_days = suppressWarnings(as.numeric(.data$pre_observed_days)),
    post_observed_days = suppressWarnings(as.numeric(.data$post_observed_days)),
    pre_seizure_count = suppressWarnings(as.numeric(.data$pre_seizure_count)),
    post_seizure_count = suppressWarnings(as.numeric(.data$post_seizure_count))
  ) %>%
  inner_join(cluster_lookup, by = "patient_id") %>%
  filter(!is.na(.data$candidate_week), .data$candidate_week <= ANALYSIS_CUTOFF_DATE) %>%
  mutate(
    pre_rate_per_28d = if_else(
      is.finite(.data$pre_observed_days) & .data$pre_observed_days > 0,
      STANDARD_MONTH_DAYS * .data$pre_seizure_count / .data$pre_observed_days,
      NA_real_
    ),
    post_rate_per_28d = if_else(
      is.finite(.data$post_observed_days) & .data$post_observed_days > 0,
      STANDARD_MONTH_DAYS * .data$post_seizure_count / .data$post_observed_days,
      NA_real_
    ),
    rate_ratio = case_when(
      !is.finite(.data$pre_rate_per_28d) | !is.finite(.data$post_rate_per_28d) ~ NA_real_,
      .data$pre_rate_per_28d == 0 & .data$post_rate_per_28d == 0 ~ 1,
      .data$pre_rate_per_28d == 0 & .data$post_rate_per_28d > 0 ~ Inf,
      .data$pre_rate_per_28d > 0 ~ .data$post_rate_per_28d / .data$pre_rate_per_28d,
      TRUE ~ NA_real_
    ),
    fold_change = if_else(
      !is.na(.data$rate_ratio),
      case_when(
        .data$rate_ratio == 0 ~ Inf,
        is.infinite(.data$rate_ratio) ~ Inf,
        .data$rate_ratio > 0 ~ pmax(.data$rate_ratio, 1 / .data$rate_ratio),
        TRUE ~ NA_real_
      ),
      NA_real_
    )
  ) %>%
  group_by(.data$pam_k3) %>%
  summarise(
    patients_with_significant_changepoint = n_distinct(.data$patient_id[.data$significant]),
    n_significant_changepoints = sum(.data$significant, na.rm = TRUE),
    median_fold_change_at_significant_changepoints = median_available(.data$fold_change[.data$significant]),
    .groups = "drop"
  )

cluster_count <- dplyr::n_distinct(cluster_lookup$pam_k3)

presentation_table <- monthly_summary %>%
  full_join(segment_summary, by = "pam_k3") %>%
  full_join(changepoint_summary, by = "pam_k3") %>%
  mutate(
    across(
      all_of(c(
        "n_patients",
        "n_changepoint_analyzable_patients",
        "n_changepoint_segments",
        "patients_with_significant_changepoint",
        "n_significant_changepoints"
      )),
      ~ tidyr::replace_na(.x, 0)
    ),
    burden_rank = rank(-.data$median_monthly_rate_per_28d, ties.method = "first"),
    interpretation = vapply(
      .data$burden_rank,
      format_interpretation,
      character(1),
      n_clusters = cluster_count
    )
  ) %>%
  transmute(
    pam_k3 = .data$pam_k3,
    n_patients = as.integer(.data$n_patients),
    median_monthly_rate_per_28d = round(.data$median_monthly_rate_per_28d, 2),
    median_patient_monthly_iqr_per_28d = round(.data$median_patient_monthly_iqr_per_28d, 2),
    median_zero_month_proportion = round(.data$median_zero_month_proportion, 2),
    median_segment_rate_per_28d = round(.data$median_segment_rate_per_28d, 2),
    segment_rate_iqr_per_28d = round(.data$segment_rate_iqr_per_28d, 2),
    patients_with_significant_changepoint = as.integer(.data$patients_with_significant_changepoint),
    n_significant_changepoints = as.integer(.data$n_significant_changepoints),
    median_fold_change_at_significant_changepoints = round(
      .data$median_fold_change_at_significant_changepoints,
      2
    ),
    interpretation = .data$interpretation
  ) %>%
  arrange(.data$pam_k3)

readr::write_csv(presentation_table, OUTPUT_PATH)
