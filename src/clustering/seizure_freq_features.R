suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(readr)
  library(tibble)
})

source("src/analysis_config.R")

INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
OUTPUT_DIR <- "output/tabs/clustering"
OUTPUT_PATH <- file.path(OUTPUT_DIR, "seizure_freq_features.csv")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_PATH)) {
  stop("Input file not found: ", INPUT_PATH)
}

required_columns <- c(
  "patient_id",
  "variant_p",
  "month",
  "study_start_date",
  "study_end_date",
  "seizure_count"
)

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

patient_month_panel <- readr::read_csv(INPUT_PATH, show_col_types = FALSE)

missing_columns <- setdiff(required_columns, names(patient_month_panel))
if (length(missing_columns) > 0) {
  stop("Input data is missing required columns: ", paste(missing_columns, collapse = ", "))
}

patient_month_rates <- patient_month_panel %>%
  mutate(
    patient_id = as.character(.data$patient_id),
    month = as.Date(.data$month),
    study_start_date = as.Date(.data$study_start_date),
    study_end_date = pmin(as.Date(.data$study_end_date), ANALYSIS_CUTOFF_DATE),
    seizure_count = suppressWarnings(as.numeric(.data$seizure_count)),
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
    .data$month <= ANALYSIS_CUTOFF_DATE,
    .data$observed_days_in_month > 0
  ) %>%
  mutate(
    monthly_seizure_rate = .data$seizure_count / .data$observed_days_in_month * STANDARD_MONTH_DAYS,
    zero_seizure_month_flag = as.integer(.data$seizure_count == 0)
  )

seizure_freq_features <- patient_month_rates %>%
  group_by(.data$patient_id) %>%
  summarise(
    variant_p = first_non_missing(.data$variant_p),
    mean_monthly_seizure_rate = mean(.data$monthly_seizure_rate),
    median_monthly_seizure_rate = stats::median(.data$monthly_seizure_rate),
    max_monthly_seizure_rate = max(.data$monthly_seizure_rate),
    iqr_monthly_seizure_rate = stats::IQR(.data$monthly_seizure_rate),
    proportion_zero_seizure_months = mean(.data$zero_seizure_month_flag),
    .groups = "drop"
  ) %>%
  arrange(desc(.data$mean_monthly_seizure_rate), .data$patient_id)

readr::write_csv(seizure_freq_features, OUTPUT_PATH)
