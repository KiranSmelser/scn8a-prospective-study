suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(readr)
})

source("src/analysis_config.R")

FEATURES_INPUT_PATH <- "output/tabs/clustering/seizure_freq_features.csv"
ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
PANEL_INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
SEIZURES_INPUT_PATH <- "output/tabs/seizures/seizures.csv"
OUTPUT_DIR <- "output/tabs/clustering"
SUMMARY_OUTPUT_PATH <- file.path(OUTPUT_DIR, "cluster_summary.csv")

FEATURE_REQUIRED_COLUMNS <- c(
  "patient_id",
  "mean_monthly_seizure_rate",
  "iqr_monthly_seizure_rate",
  "proportion_zero_seizure_months"
)
ASSIGNMENT_COLUMNS <- c("pam_k3", "pam_k3_tonic_clonic", "pam_k3_focal")
ASSIGNMENT_SEIZURE_TYPES <- c(
  pam_k3 = NA_character_,
  pam_k3_tonic_clonic = "Tonic-clonic",
  pam_k3_focal = "Focal"
)
ASSIGNMENT_LABELS <- c(
  pam_k3 = "all seizure types",
  pam_k3_tonic_clonic = "tonic-clonic seizures",
  pam_k3_focal = "focal seizures"
)
ASSIGNMENT_REQUIRED_COLUMNS <- c("patient_id", ASSIGNMENT_COLUMNS)
PANEL_REQUIRED_COLUMNS <- c("patient_id", "variant_p", "month", "study_start_date", "study_end_date")
SEIZURE_REQUIRED_COLUMNS <- c("patient_id", "type", "date")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(FEATURES_INPUT_PATH)) {
  stop("Feature input file not found: ", FEATURES_INPUT_PATH)
}

if (!file.exists(ASSIGNMENTS_INPUT_PATH)) {
  stop("Assignment input file not found: ", ASSIGNMENTS_INPUT_PATH)
}

if (!file.exists(PANEL_INPUT_PATH)) {
  stop("Patient-month panel file not found: ", PANEL_INPUT_PATH)
}

if (!file.exists(SEIZURES_INPUT_PATH)) {
  stop("Seizure input file not found: ", SEIZURES_INPUT_PATH)
}

features <- readr::read_csv(FEATURES_INPUT_PATH, show_col_types = FALSE)
assignments <- readr::read_csv(ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE)
patient_month_panel <- readr::read_csv(PANEL_INPUT_PATH, show_col_types = FALSE)
seizure_events <- readr::read_csv(SEIZURES_INPUT_PATH, show_col_types = FALSE)

missing_feature_columns <- setdiff(FEATURE_REQUIRED_COLUMNS, names(features))
if (length(missing_feature_columns) > 0) {
  stop("Feature input is missing required columns: ", paste(missing_feature_columns, collapse = ", "))
}

missing_assignment_columns <- setdiff(ASSIGNMENT_REQUIRED_COLUMNS, names(assignments))
if (length(missing_assignment_columns) > 0) {
  stop("Assignment input is missing required columns: ", paste(missing_assignment_columns, collapse = ", "))
}

missing_panel_columns <- setdiff(PANEL_REQUIRED_COLUMNS, names(patient_month_panel))
if (length(missing_panel_columns) > 0) {
  stop("Patient-month panel is missing required columns: ", paste(missing_panel_columns, collapse = ", "))
}

missing_seizure_columns <- setdiff(SEIZURE_REQUIRED_COLUMNS, names(seizure_events))
if (length(missing_seizure_columns) > 0) {
  stop("Seizure input is missing required columns: ", paste(missing_seizure_columns, collapse = ", "))
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

build_feature_data <- function(patient_month_panel, monthly_counts) {
  build_monthly_rates(patient_month_panel, monthly_counts) %>%
    group_by(.data$patient_id) %>%
    summarise(
      mean_monthly_seizure_rate = mean(.data$monthly_seizure_rate),
      iqr_monthly_seizure_rate = stats::IQR(.data$monthly_seizure_rate),
      proportion_zero_seizure_months = mean(.data$zero_seizure_month_flag),
      .groups = "drop"
    )
}

build_type_feature_data <- function(patient_month_panel, seizure_events, seizure_type) {
  type_events <- seizure_events %>%
    filter(.data$type == seizure_type)

  eligible_patients <- type_events %>%
    distinct(.data$patient_id)

  monthly_counts <- type_events %>%
    mutate(month = as.Date(lubridate::floor_date(.data$event_date, "month"))) %>%
    count(.data$patient_id, .data$month, name = "seizure_count")

  patient_month_panel %>%
    semi_join(eligible_patients, by = "patient_id") %>%
    build_feature_data(monthly_counts)
}

all_feature_data <- features %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    mean_monthly_seizure_rate = suppressWarnings(as.numeric(.data$mean_monthly_seizure_rate)),
    iqr_monthly_seizure_rate = suppressWarnings(as.numeric(.data$iqr_monthly_seizure_rate)),
    proportion_zero_seizure_months = suppressWarnings(as.numeric(.data$proportion_zero_seizure_months))
  ) %>%
  filter(
    !is.na(.data$patient_id),
    if_all(
      all_of(c(
        "mean_monthly_seizure_rate",
        "iqr_monthly_seizure_rate",
        "proportion_zero_seizure_months"
      )),
      ~ !is.na(.x)
    )
  )

patient_month_panel <- patient_month_panel %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    variant_p = as.character(.data$variant_p),
    month = as.Date(.data$month),
    study_start_date = as.Date(.data$study_start_date),
    study_end_date = as.Date(.data$study_end_date)
  ) %>%
  filter(
    !is.na(.data$patient_id),
    !is.na(.data$month),
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

summarise_cluster_column <- function(cluster_column) {
  seizure_type <- unname(ASSIGNMENT_SEIZURE_TYPES[[cluster_column]])
  feature_data <- if (is.na(seizure_type)) {
    all_feature_data
  } else {
    build_type_feature_data(patient_month_panel, seizure_events, seizure_type)
  }

  cluster_lookup <- assignments %>%
    transmute(
      patient_id = as.character(.data$patient_id),
      cluster = suppressWarnings(as.integer(.data[[cluster_column]]))
    ) %>%
    filter(!is.na(.data$patient_id), !is.na(.data$cluster))

  feature_data %>%
    inner_join(cluster_lookup, by = "patient_id") %>%
    group_by(.data$cluster) %>%
    summarise(
      n_patients = n(),
      mean_of_mean_monthly_seizure_rate = mean(.data$mean_monthly_seizure_rate),
      median_of_mean_monthly_seizure_rate = stats::median(.data$mean_monthly_seizure_rate),
      mean_iqr_monthly_seizure_rate = mean(.data$iqr_monthly_seizure_rate),
      iqr_monthly_seizure_rate_median = stats::median(.data$iqr_monthly_seizure_rate),
      mean_proportion_zero_seizure_months = mean(.data$proportion_zero_seizure_months),
      .groups = "drop"
    ) %>%
    mutate(
      cluster_type = unname(ASSIGNMENT_LABELS[[cluster_column]]),
      cluster_order = match(cluster_column, ASSIGNMENT_COLUMNS),
      .before = "cluster"
    )
}

cluster_summary <- dplyr::bind_rows(
  lapply(ASSIGNMENT_COLUMNS, summarise_cluster_column)
) %>%
  arrange(.data$cluster_order, .data$cluster) %>%
  select(-"cluster_order")

readr::write_csv(cluster_summary, SUMMARY_OUTPUT_PATH)
