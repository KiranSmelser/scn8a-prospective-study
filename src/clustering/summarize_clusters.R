suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

FEATURES_INPUT_PATH <- "output/tabs/clustering/seizure_freq_features.csv"
ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
OUTPUT_DIR <- "output/tabs/clustering"
SUMMARY_OUTPUT_PATH <- file.path(OUTPUT_DIR, "cluster_summary.csv")

FEATURE_REQUIRED_COLUMNS <- c(
  "patient_id",
  "mean_monthly_seizure_rate",
  "iqr_monthly_seizure_rate",
  "proportion_zero_seizure_months"
)
ASSIGNMENT_REQUIRED_COLUMNS <- c("patient_id", "pam_k3")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(FEATURES_INPUT_PATH)) {
  stop("Feature input file not found: ", FEATURES_INPUT_PATH)
}

if (!file.exists(ASSIGNMENTS_INPUT_PATH)) {
  stop("Assignment input file not found: ", ASSIGNMENTS_INPUT_PATH)
}

features <- readr::read_csv(FEATURES_INPUT_PATH, show_col_types = FALSE)
assignments <- readr::read_csv(ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE)

missing_feature_columns <- setdiff(FEATURE_REQUIRED_COLUMNS, names(features))
if (length(missing_feature_columns) > 0) {
  stop("Feature input is missing required columns: ", paste(missing_feature_columns, collapse = ", "))
}

missing_assignment_columns <- setdiff(ASSIGNMENT_REQUIRED_COLUMNS, names(assignments))
if (length(missing_assignment_columns) > 0) {
  stop("Assignment input is missing required columns: ", paste(missing_assignment_columns, collapse = ", "))
}

clustered_patients <- features %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    mean_monthly_seizure_rate = suppressWarnings(as.numeric(.data$mean_monthly_seizure_rate)),
    iqr_monthly_seizure_rate = suppressWarnings(as.numeric(.data$iqr_monthly_seizure_rate)),
    proportion_zero_seizure_months = suppressWarnings(as.numeric(.data$proportion_zero_seizure_months))
  ) %>%
  inner_join(
    assignments %>%
      transmute(
        patient_id = as.character(.data$patient_id),
        pam_k3 = suppressWarnings(as.integer(.data$pam_k3))
      ),
    by = "patient_id"
  ) %>%
  filter(
    !is.na(.data$patient_id),
    !is.na(.data$pam_k3),
    if_all(
      all_of(c(
        "mean_monthly_seizure_rate",
        "iqr_monthly_seizure_rate",
        "proportion_zero_seizure_months"
      )),
      ~ !is.na(.x)
    )
  )

cluster_summary <- clustered_patients %>%
  group_by(.data$pam_k3) %>%
  summarise(
    n_patients = n(),
    mean_of_mean_monthly_seizure_rate = mean(.data$mean_monthly_seizure_rate),
    median_of_mean_monthly_seizure_rate = stats::median(.data$mean_monthly_seizure_rate),
    mean_iqr_monthly_seizure_rate = mean(.data$iqr_monthly_seizure_rate),
    iqr_monthly_seizure_rate_median = stats::median(.data$iqr_monthly_seizure_rate),
    mean_proportion_zero_seizure_months = mean(.data$proportion_zero_seizure_months),
    .groups = "drop"
  ) %>%
  arrange(.data$pam_k3)

readr::write_csv(cluster_summary, SUMMARY_OUTPUT_PATH)
