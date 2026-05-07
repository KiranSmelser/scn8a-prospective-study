# Patient-level changepoint-stratified permutation tests for post-seizure risk.

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
PATIENT_RESULTS_OUTPUT_PATH <- file.path(
  OUTPUT_DIR,
  "patient_level_changepoint_stratified_permutation.csv"
)
DAILY_PANEL_OUTPUT_PATH <- file.path(OUTPUT_DIR, "patient_day_post_seizure_panel.csv")
SEGMENTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "patient_day_changepoint_segments.csv")

EXPECTED_COHORT_SIZE <- as.integer(Sys.getenv("PERIODICITY_EXPECTED_COHORT_SIZE", "23"))
PERMUTATION_REPLICATES <- as.integer(Sys.getenv("PERIODICITY_PERMUTATIONS", "1000"))
RANDOM_SEED <- as.integer(Sys.getenv("PERIODICITY_RANDOM_SEED", "20260507"))
ALPHA <- 0.05

POST_SEIZURE_WINDOWS <- tibble::tibble(
  window = c("lag_1_day", "lag_2_3_days", "lag_4_7_days", "lag_8_14_days"),
  lag_start_days = c(1L, 2L, 4L, 8L),
  lag_end_days = c(1L, 3L, 7L, 14L)
)

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

if (!is.finite(PERMUTATION_REPLICATES) || PERMUTATION_REPLICATES < 1L) {
  stop("PERIODICITY_PERMUTATIONS must be a positive integer.", call. = FALSE)
}

parse_event_date <- function(x) {
  parsed <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC"))
  parsed <- dplyr::if_else(is.na(parsed), suppressWarnings(lubridate::ymd(x)), parsed)
  as.Date(parsed)
}

make_post_seizure_exposure <- function(seizure_day, lag_start_days, lag_end_days) {
  exposure <- rep(FALSE, length(seizure_day))
  n_days <- length(seizure_day)

  for (lag_days in seq.int(lag_start_days, lag_end_days)) {
    if (lag_days < n_days) {
      exposure[(lag_days + 1L):n_days] <- exposure[(lag_days + 1L):n_days] | seizure_day[1L:(n_days - lag_days)]
    }
  }

  exposure
}

safe_risk_ratio <- function(exposed_risk, unexposed_risk) {
  dplyr::case_when(
    is.na(exposed_risk) | is.na(unexposed_risk) ~ NA_real_,
    exposed_risk == 0 & unexposed_risk == 0 ~ 1,
    exposed_risk > 0 & unexposed_risk == 0 ~ Inf,
    TRUE ~ exposed_risk / unexposed_risk
  )
}

window_statistic <- function(seizure_day, lag_start_days, lag_end_days) {
  exposure <- make_post_seizure_exposure(seizure_day, lag_start_days, lag_end_days)
  exposed_days <- sum(exposure)
  unexposed_days <- sum(!exposure)
  exposed_seizure_days <- sum(seizure_day & exposure)
  unexposed_seizure_days <- sum(seizure_day & !exposure)
  exposed_risk <- if (exposed_days > 0) exposed_seizure_days / exposed_days else NA_real_
  unexposed_risk <- if (unexposed_days > 0) unexposed_seizure_days / unexposed_days else NA_real_

  tibble::tibble(
    exposed_days = exposed_days,
    unexposed_days = unexposed_days,
    exposed_seizure_days = exposed_seizure_days,
    unexposed_seizure_days = unexposed_seizure_days,
    exposed_risk = exposed_risk,
    unexposed_risk = unexposed_risk,
    risk_difference = exposed_risk - unexposed_risk,
    risk_ratio = safe_risk_ratio(exposed_risk, unexposed_risk)
  )
}

risk_difference_statistic <- function(seizure_day, lag_start_days, lag_end_days) {
  exposure <- make_post_seizure_exposure(seizure_day, lag_start_days, lag_end_days)
  exposed_days <- sum(exposure)
  unexposed_days <- sum(!exposure)

  if (exposed_days == 0 || unexposed_days == 0) {
    return(NA_real_)
  }

  exposed_risk <- sum(seizure_day & exposure) / exposed_days
  unexposed_risk <- sum(seizure_day & !exposure) / unexposed_days
  exposed_risk - unexposed_risk
}

all_window_statistics <- function(seizure_day) {
  POST_SEIZURE_WINDOWS %>%
    dplyr::mutate(
      statistic = purrr::map2(
        .data$lag_start_days,
        .data$lag_end_days,
        ~ window_statistic(seizure_day, .x, .y)
      )
    ) %>%
    tidyr::unnest("statistic")
}

permute_within_segments <- function(seizure_day, segment_id) {
  permuted <- seizure_day
  segment_levels <- unique(segment_id)

  for (segment_value in segment_levels) {
    index <- which(segment_id == segment_value)
    permuted[index] <- sample(seizure_day[index], length(index), replace = FALSE)
  }

  permuted
}

patient_permutation_tests <- function(patient_daily_data) {
  observed <- all_window_statistics(patient_daily_data$seizure_day)
  seizure_day <- patient_daily_data$seizure_day
  segment_id <- patient_daily_data$segment_id

  if (sum(seizure_day) == 0) {
    return(observed %>%
      dplyr::mutate(
        permutation_replicates = PERMUTATION_REPLICATES,
        valid_permutation_replicates = 0L,
        null_risk_difference_mean = NA_real_,
        null_risk_difference_lower_95 = NA_real_,
        null_risk_difference_upper_95 = NA_real_,
        p_value_increased_risk = NA_real_,
        p_value_decreased_risk = NA_real_,
        p_value_two_sided = NA_real_,
        analysis_status = "not_run_no_seizure_days"
      ))
  }

  permutation_risk_differences <- matrix(
    NA_real_,
    nrow = PERMUTATION_REPLICATES,
    ncol = nrow(POST_SEIZURE_WINDOWS),
    dimnames = list(NULL, POST_SEIZURE_WINDOWS$window)
  )

  for (replicate_id in seq_len(PERMUTATION_REPLICATES)) {
    permuted_seizure_day <- permute_within_segments(seizure_day, segment_id)

    for (window_index in seq_len(nrow(POST_SEIZURE_WINDOWS))) {
      permutation_risk_differences[replicate_id, window_index] <- risk_difference_statistic(
        permuted_seizure_day,
        POST_SEIZURE_WINDOWS$lag_start_days[[window_index]],
        POST_SEIZURE_WINDOWS$lag_end_days[[window_index]]
      )
    }
  }

  permutation_summary <- purrr::map_dfr(seq_len(nrow(POST_SEIZURE_WINDOWS)), function(window_index) {
    null_values <- permutation_risk_differences[, window_index]
    valid_values <- null_values[is.finite(null_values)]
    observed_value <- observed$risk_difference[[window_index]]
    valid_count <- length(valid_values)

    tibble::tibble(
      window = POST_SEIZURE_WINDOWS$window[[window_index]],
      valid_permutation_replicates = valid_count,
      null_risk_difference_mean = if (valid_count > 0) mean(valid_values) else NA_real_,
      null_risk_difference_lower_95 = if (valid_count > 0) as.numeric(stats::quantile(valid_values, probs = 0.025, names = FALSE)) else NA_real_,
      null_risk_difference_upper_95 = if (valid_count > 0) as.numeric(stats::quantile(valid_values, probs = 0.975, names = FALSE)) else NA_real_,
      permuted_greater_equal = if (valid_count > 0 && is.finite(observed_value)) sum(valid_values >= observed_value) else NA_integer_,
      permuted_less_equal = if (valid_count > 0 && is.finite(observed_value)) sum(valid_values <= observed_value) else NA_integer_
    )
  })

  observed %>%
    dplyr::left_join(permutation_summary, by = "window") %>%
    dplyr::mutate(
      permutation_replicates = PERMUTATION_REPLICATES,
      p_value_increased_risk = dplyr::if_else(
        is.finite(.data$risk_difference) & .data$valid_permutation_replicates > 0,
        (1 + .data$permuted_greater_equal) / (1 + .data$valid_permutation_replicates),
        NA_real_
      ),
      p_value_decreased_risk = dplyr::if_else(
        is.finite(.data$risk_difference) & .data$valid_permutation_replicates > 0,
        (1 + .data$permuted_less_equal) / (1 + .data$valid_permutation_replicates),
        NA_real_
      ),
      p_value_two_sided = pmin(1, 2 * pmin(.data$p_value_increased_risk, .data$p_value_decreased_risk)),
      analysis_status = dplyr::case_when(
        !is.finite(.data$risk_difference) ~ "not_run_insufficient_exposed_or_unexposed_days",
        .data$valid_permutation_replicates == 0 ~ "not_run_no_valid_permutations",
        TRUE ~ "ok"
      )
    ) %>%
    dplyr::select(-"permuted_greater_equal", -"permuted_less_equal")
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

daily_panel_with_exposures <- daily_panel_segmented %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::arrange(.data$date, .by_group = TRUE) %>%
  dplyr::group_modify(function(patient_data, patient_key) {
    exposure_columns <- POST_SEIZURE_WINDOWS %>%
      dplyr::mutate(
        exposure = purrr::map2(
          .data$lag_start_days,
          .data$lag_end_days,
          ~ make_post_seizure_exposure(patient_data$seizure_day, .x, .y)
        ),
        column_name = paste0("post_seizure_", .data$window)
      )

    for (row_index in seq_len(nrow(exposure_columns))) {
      patient_data[[exposure_columns$column_name[[row_index]]]] <- exposure_columns$exposure[[row_index]]
    }

    patient_data
  }) %>%
  dplyr::ungroup()

set.seed(RANDOM_SEED)

patient_results <- daily_panel_segmented %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::group_modify(function(.x, .y) {
    patient_permutation_tests(.x) %>%
      dplyr::mutate(
        observed_days = nrow(.x),
        observed_seizure_days = sum(.x$seizure_day),
        observed_seizure_events = sum(.x$seizure_count),
        n_segments = dplyr::n_distinct(.x$segment_id),
        segment_sources = paste(sort(unique(.x$segment_source)), collapse = ";")
      )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(patient_windows %>% dplyr::select("patient_id", "variant_p", "study_start_date", "study_end_date"), by = "patient_id") %>%
  dplyr::mutate(
    q_value = NA_real_,
    significant = FALSE,
    direction = dplyr::case_when(
      is.na(.data$risk_difference) ~ NA_character_,
      .data$risk_difference > 0 ~ "increased",
      .data$risk_difference < 0 ~ "decreased",
      TRUE ~ "no_difference"
    )
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
    "segment_sources",
    "window",
    "lag_start_days",
    "lag_end_days"
  )

valid_p_values <- is.finite(patient_results$p_value_two_sided)
patient_results$q_value[valid_p_values] <- stats::p.adjust(
  patient_results$p_value_two_sided[valid_p_values],
  method = "BH"
)
patient_results$significant <- !is.na(patient_results$q_value) & patient_results$q_value < ALPHA

readr::write_csv(analysis_segments, SEGMENTS_OUTPUT_PATH)
readr::write_csv(daily_panel_with_exposures, DAILY_PANEL_OUTPUT_PATH)
readr::write_csv(patient_results, PATIENT_RESULTS_OUTPUT_PATH)
