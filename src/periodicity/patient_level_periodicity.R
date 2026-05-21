# Patient-level changepoint-stratified permutation tests for seizure periodicity.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(lubridate)
})

PANEL_INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
SEIZURE_INPUT_PATH <- "output/tabs/seizures/seizures.csv"
CHANGEPOINT_SEGMENTS_INPUT_PATH <- "output/tabs/changepoints/change_point_segments.csv"
OUTPUT_DIR <- "output/tabs/periodicity"
PATIENT_RESULTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "patient_level_periodicity.csv")
LAG_RESULTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "patient_level_periodicity_lags.csv")

EXPECTED_COHORT_SIZE <- as.integer(Sys.getenv("PERIODICITY_EXPECTED_COHORT_SIZE", "23"))
MAX_LAG_DAYS <- as.integer(Sys.getenv("PERIODICITY_MAX_LAG_DAYS", "90"))
MIN_LAG_DAYS <- as.integer(Sys.getenv("PERIODICITY_MIN_LAG_DAYS", "1"))
PERMUTATION_REPLICATES <- as.integer(Sys.getenv("PERIODICITY_ACF_PERMUTATIONS", "1000"))
RANDOM_SEED <- as.integer(Sys.getenv("PERIODICITY_ACF_RANDOM_SEED", "20260521"))
ALPHA <- 0.05

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PANEL_INPUT_PATH)) {
  stop("Input file not found: ", PANEL_INPUT_PATH, call. = FALSE)
}

if (!file.exists(SEIZURE_INPUT_PATH)) {
  stop("Input file not found: ", SEIZURE_INPUT_PATH, call. = FALSE)
}

if (!file.exists(CHANGEPOINT_SEGMENTS_INPUT_PATH)) {
  stop("Input file not found: ", CHANGEPOINT_SEGMENTS_INPUT_PATH, call. = FALSE)
}

if (!is.finite(MAX_LAG_DAYS) || MAX_LAG_DAYS < 1L) {
  stop("PERIODICITY_MAX_LAG_DAYS must be a positive integer.", call. = FALSE)
}

if (!is.finite(MIN_LAG_DAYS) || MIN_LAG_DAYS < 1L) {
  stop("PERIODICITY_MIN_LAG_DAYS must be a positive integer.", call. = FALSE)
}

if (MIN_LAG_DAYS > MAX_LAG_DAYS) {
  stop("PERIODICITY_MIN_LAG_DAYS cannot exceed PERIODICITY_MAX_LAG_DAYS.", call. = FALSE)
}

if (!is.finite(PERMUTATION_REPLICATES) || PERMUTATION_REPLICATES < 1L) {
  stop("PERIODICITY_ACF_PERMUTATIONS must be a positive integer.", call. = FALSE)
}

parse_event_date <- function(x) {
  parsed <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC"))
  parsed <- dplyr::if_else(is.na(parsed), suppressWarnings(lubridate::ymd(x)), parsed)
  as.Date(parsed)
}

autocorrelation_by_lag <- function(seizure_count, min_lag_days, max_lag_days) {
  seizure_count <- as.numeric(seizure_count)
  n_days <- length(seizure_count)
  max_lag <- min(max_lag_days, n_days - 1L)

  if (max_lag < min_lag_days) {
    return(tibble::tibble(lag_days = integer(), autocorrelation = numeric()))
  }

  lag_sequence <- seq.int(min_lag_days, max_lag)

  centered_counts <- seizure_count - mean(seizure_count, na.rm = TRUE)
  denominator <- sum(centered_counts^2)

  if (!is.finite(denominator) || denominator <= 0) {
    return(tibble::tibble(
      lag_days = lag_sequence,
      autocorrelation = rep(NA_real_, length(lag_sequence))
    ))
  }

  tibble::tibble(
    lag_days = lag_sequence,
    autocorrelation = purrr::map_dbl(.data$lag_days, function(lag_days) {
      numerator <- sum(
        centered_counts[seq_len(n_days - lag_days)] *
          centered_counts[(lag_days + 1L):n_days]
      )
      numerator / denominator
    })
  )
}

strongest_lags <- function(lag_results) {
  finite_lags <- lag_results %>%
    dplyr::filter(is.finite(.data$autocorrelation))

  if (nrow(finite_lags) == 0) {
    return(tibble::tibble(
      strongest_lag_days = NA_integer_,
      strongest_autocorrelation = NA_real_,
      strongest_positive_lag_days = NA_integer_,
      strongest_positive_autocorrelation = NA_real_
    ))
  }

  strongest <- finite_lags %>%
    dplyr::arrange(dplyr::desc(.data$autocorrelation), .data$lag_days) %>%
    dplyr::slice(1L)

  strongest_positive <- finite_lags %>%
    dplyr::filter(.data$autocorrelation > 0) %>%
    dplyr::arrange(dplyr::desc(.data$autocorrelation), .data$lag_days) %>%
    dplyr::slice(1L)

  tibble::tibble(
    strongest_lag_days = strongest$lag_days[[1]],
    strongest_autocorrelation = strongest$autocorrelation[[1]],
    strongest_positive_lag_days = if (nrow(strongest_positive) > 0) strongest_positive$lag_days[[1]] else NA_integer_,
    strongest_positive_autocorrelation = if (nrow(strongest_positive) > 0) strongest_positive$autocorrelation[[1]] else NA_real_
  )
}

max_positive_autocorrelation <- function(seizure_count, min_lag_days, max_lag_days) {
  lag_results <- autocorrelation_by_lag(seizure_count, min_lag_days, max_lag_days)
  positive_values <- lag_results$autocorrelation[
    is.finite(lag_results$autocorrelation) & lag_results$autocorrelation > 0
  ]

  if (length(positive_values) == 0) {
    return(NA_real_)
  }

  max(positive_values)
}

permute_counts_within_segments <- function(seizure_count, segment_id) {
  permuted <- seizure_count
  segment_levels <- unique(segment_id)

  for (segment_value in segment_levels) {
    index <- which(segment_id == segment_value)
    permuted[index] <- sample(seizure_count[index], length(index), replace = FALSE)
  }

  permuted
}

patient_periodicity_test <- function(patient_daily_data) {
  seizure_count <- patient_daily_data$seizure_count
  segment_id <- patient_daily_data$segment_id
  observed_lags <- autocorrelation_by_lag(seizure_count, MIN_LAG_DAYS, MAX_LAG_DAYS)
  observed_strongest <- strongest_lags(observed_lags)

  base_summary <- tibble::tibble(
    observed_days = nrow(patient_daily_data),
    observed_seizure_days = sum(patient_daily_data$seizure_day),
    observed_seizure_events = sum(seizure_count),
    n_segments = dplyr::n_distinct(segment_id),
    segment_sources = paste(sort(unique(patient_daily_data$segment_source)), collapse = ";"),
    min_lag_days = MIN_LAG_DAYS,
    max_lag_days = MAX_LAG_DAYS,
    lags_tested = nrow(observed_lags),
    permutation_replicates = PERMUTATION_REPLICATES
  ) %>%
    dplyr::bind_cols(observed_strongest)

  if (sum(seizure_count) == 0) {
    return(base_summary %>%
      dplyr::mutate(
        valid_permutation_replicates = 0L,
        null_max_positive_autocorrelation_mean = NA_real_,
        null_max_positive_autocorrelation_lower_95 = NA_real_,
        null_max_positive_autocorrelation_upper_95 = NA_real_,
        p_value = NA_real_,
        analysis_status = "not_run_no_seizure_events"
      ))
  }

  if (!is.finite(stats::sd(seizure_count)) || stats::sd(seizure_count) == 0) {
    return(base_summary %>%
      dplyr::mutate(
        valid_permutation_replicates = 0L,
        null_max_positive_autocorrelation_mean = NA_real_,
        null_max_positive_autocorrelation_lower_95 = NA_real_,
        null_max_positive_autocorrelation_upper_95 = NA_real_,
        p_value = NA_real_,
        analysis_status = "not_run_no_count_variation"
      ))
  }

  observed_statistic <- observed_strongest$strongest_positive_autocorrelation[[1]]

  if (!is.finite(observed_statistic)) {
    return(base_summary %>%
      dplyr::mutate(
        valid_permutation_replicates = 0L,
        null_max_positive_autocorrelation_mean = NA_real_,
        null_max_positive_autocorrelation_lower_95 = NA_real_,
        null_max_positive_autocorrelation_upper_95 = NA_real_,
        p_value = NA_real_,
        analysis_status = "not_run_no_positive_autocorrelation"
      ))
  }

  null_statistics <- purrr::map_dbl(seq_len(PERMUTATION_REPLICATES), function(replicate_id) {
    permuted_counts <- permute_counts_within_segments(seizure_count, segment_id)
    max_positive_autocorrelation(permuted_counts, MIN_LAG_DAYS, MAX_LAG_DAYS)
  })
  valid_null_statistics <- null_statistics[is.finite(null_statistics)]
  valid_count <- length(valid_null_statistics)

  base_summary %>%
    dplyr::mutate(
      valid_permutation_replicates = valid_count,
      null_max_positive_autocorrelation_mean = if (valid_count > 0) mean(valid_null_statistics) else NA_real_,
      null_max_positive_autocorrelation_lower_95 = if (valid_count > 0) as.numeric(stats::quantile(valid_null_statistics, probs = 0.025, names = FALSE)) else NA_real_,
      null_max_positive_autocorrelation_upper_95 = if (valid_count > 0) as.numeric(stats::quantile(valid_null_statistics, probs = 0.975, names = FALSE)) else NA_real_,
      p_value = dplyr::if_else(
        valid_count > 0,
        (1 + sum(valid_null_statistics >= observed_statistic)) / (1 + valid_count),
        NA_real_
      ),
      analysis_status = dplyr::if_else(valid_count > 0, "ok", "not_run_no_valid_permutations")
    )
}

patient_windows <- readr::read_csv(PANEL_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    study_start_date = as.Date(.data$study_start_date),
    study_end_date = as.Date(.data$study_end_date)
  ) %>%
  dplyr::filter(!is.na(.data$study_start_date), !is.na(.data$study_end_date), .data$study_start_date <= .data$study_end_date) %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(
    variant_p = dplyr::first(.data$variant_p),
    study_start_date = min(.data$study_start_date),
    study_end_date = max(.data$study_end_date),
    .groups = "drop"
  ) %>%
  dplyr::arrange(.data$patient_id)

if (nrow(patient_windows) != EXPECTED_COHORT_SIZE) {
  stop(
    "Expected ",
    EXPECTED_COHORT_SIZE,
    " analysis patients from ",
    PANEL_INPUT_PATH,
    " but found ",
    nrow(patient_windows),
    ".",
    call. = FALSE
  )
}

seizure_events <- readr::read_csv(SEIZURE_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    event_date = parse_event_date(.data$date)
  ) %>%
  dplyr::inner_join(patient_windows, by = "patient_id") %>%
  dplyr::filter(
    !is.na(.data$event_date),
    .data$event_date >= .data$study_start_date,
    .data$event_date <= .data$study_end_date
  )

daily_seizure_counts <- seizure_events %>%
  dplyr::count(.data$patient_id, date = .data$event_date, name = "seizure_count")

daily_panel <- patient_windows %>%
  dplyr::mutate(date = purrr::map2(.data$study_start_date, .data$study_end_date, ~ seq(.x, .y, by = "day"))) %>%
  tidyr::unnest("date") %>%
  dplyr::left_join(daily_seizure_counts, by = c("patient_id", "date")) %>%
  tidyr::replace_na(list(seizure_count = 0L)) %>%
  dplyr::mutate(
    seizure_count = as.integer(.data$seizure_count),
    seizure_day = .data$seizure_count > 0
  )

changepoint_segments <- readr::read_csv(CHANGEPOINT_SEGMENTS_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    segment_id = as.integer(.data$segment_id),
    segment_start_date = as.Date(.data$segment_start_date),
    segment_end_date = as.Date(.data$segment_end_date),
    segment_source = "changepoint"
  ) %>%
  dplyr::filter(
    .data$patient_id %in% patient_windows$patient_id,
    !is.na(.data$segment_start_date),
    !is.na(.data$segment_end_date),
    .data$segment_start_date <= .data$segment_end_date
  )

patients_with_changepoint_segments <- unique(changepoint_segments$patient_id)

full_window_segments <- patient_windows %>%
  dplyr::filter(!(.data$patient_id %in% patients_with_changepoint_segments)) %>%
  dplyr::transmute(
    patient_id = .data$patient_id,
    segment_id = 1L,
    segment_start_date = .data$study_start_date,
    segment_end_date = .data$study_end_date,
    segment_source = "full_study_window"
  )

analysis_segments <- dplyr::bind_rows(changepoint_segments, full_window_segments) %>%
  dplyr::inner_join(
    patient_windows %>% dplyr::select("patient_id", "study_start_date", "study_end_date"),
    by = "patient_id"
  ) %>%
  dplyr::mutate(
    segment_start_date = pmax(.data$segment_start_date, .data$study_start_date),
    segment_end_date = pmin(.data$segment_end_date, .data$study_end_date)
  ) %>%
  dplyr::filter(.data$segment_start_date <= .data$segment_end_date) %>%
  dplyr::select(-"study_start_date", -"study_end_date") %>%
  dplyr::arrange(.data$patient_id, .data$segment_start_date, .data$segment_id)

daily_panel_segmented <- daily_panel %>%
  dplyr::inner_join(analysis_segments, by = "patient_id", relationship = "many-to-many") %>%
  dplyr::filter(.data$date >= .data$segment_start_date, .data$date <= .data$segment_end_date) %>%
  dplyr::arrange(.data$patient_id, .data$date, .data$segment_id) %>%
  dplyr::group_by(.data$patient_id, .data$date) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup()

if (nrow(daily_panel_segmented) != nrow(daily_panel)) {
  stop(
    "Daily panel segmentation changed the number of patient-days from ",
    nrow(daily_panel),
    " to ",
    nrow(daily_panel_segmented),
    ". Check changepoint segment coverage.",
    call. = FALSE
  )
}

set.seed(RANDOM_SEED)

patient_results <- daily_panel_segmented %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::group_modify(function(patient_data, patient_key) {
    patient_periodicity_test(patient_data)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(patient_windows, by = "patient_id") %>%
  dplyr::mutate(
    q_value = NA_real_,
    significant = FALSE
  ) %>%
  dplyr::relocate(
    "patient_id",
    "variant_p",
    "study_start_date",
    "study_end_date",
    "observed_days",
    "observed_seizure_days",
    "observed_seizure_events",
    "n_segments",
    "segment_sources"
  )

valid_p_values <- is.finite(patient_results$p_value)
patient_results$q_value[valid_p_values] <- stats::p.adjust(
  patient_results$p_value[valid_p_values],
  method = "BH"
)
patient_results$significant <- !is.na(patient_results$q_value) & patient_results$q_value < ALPHA

lag_results <- daily_panel_segmented %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::group_modify(function(patient_data, patient_key) {
    autocorrelation_by_lag(patient_data$seizure_count, MIN_LAG_DAYS, MAX_LAG_DAYS)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(
    patient_results %>%
      dplyr::select(
        "patient_id",
        "variant_p",
        "strongest_positive_lag_days",
        "strongest_positive_autocorrelation",
        "p_value",
        "q_value",
        "significant",
        "analysis_status"
      ),
    by = "patient_id"
  ) %>%
  dplyr::mutate(
    selected_strongest_positive_lag = .data$lag_days == .data$strongest_positive_lag_days
  ) %>%
  dplyr::relocate("patient_id", "variant_p", "lag_days", "autocorrelation")

readr::write_csv(patient_results, PATIENT_RESULTS_OUTPUT_PATH)
readr::write_csv(lag_results, LAG_RESULTS_OUTPUT_PATH)
