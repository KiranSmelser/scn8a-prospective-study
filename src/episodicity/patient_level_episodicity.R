# Patient-level changepoint-stratified seizure episodicity analysis.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(lubridate)
})

source("src/analysis_config.R")

PANEL_INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
SEIZURE_INPUT_PATH <- "output/tabs/seizures/seizures.csv"
CHANGEPOINT_SEGMENTS_INPUT_PATH <- "output/tabs/changepoints/change_point_segments.csv"
OUTPUT_DIR <- "output/tabs/episodicity"
PATIENT_RESULTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "patient_level_episodicity.csv")
OBSERVED_CLUMPS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "patient_level_observed_clumps.csv")
DAILY_PANEL_OUTPUT_PATH <- file.path(OUTPUT_DIR, "patient_day_episodicity_panel.csv")

EXPECTED_COHORT_SIZE <- as.integer(Sys.getenv("EPISODICITY_EXPECTED_COHORT_SIZE", "23"))
TAU_DAYS <- as.integer(strsplit(Sys.getenv("EPISODICITY_TAU_DAYS", "1,2,3"), ",")[[1]])
PRIMARY_TAU_DAYS <- as.integer(Sys.getenv("EPISODICITY_PRIMARY_TAU_DAYS", "1"))
PERMUTATION_REPLICATES <- as.integer(Sys.getenv("EPISODICITY_PERMUTATIONS", "1000"))
RANDOM_SEED <- as.integer(Sys.getenv("EPISODICITY_RANDOM_SEED", "20260708"))
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

if (length(TAU_DAYS) == 0 || any(!is.finite(TAU_DAYS)) || any(TAU_DAYS < 0L)) {
  stop("EPISODICITY_TAU_DAYS must contain non-negative integer day values.", call. = FALSE)
}

TAU_DAYS <- sort(unique(TAU_DAYS))

if (!PRIMARY_TAU_DAYS %in% TAU_DAYS) {
  stop("EPISODICITY_PRIMARY_TAU_DAYS must be one of EPISODICITY_TAU_DAYS.", call. = FALSE)
}

if (!is.finite(PERMUTATION_REPLICATES) || PERMUTATION_REPLICATES < 1L) {
  stop("EPISODICITY_PERMUTATIONS must be a positive integer.", call. = FALSE)
}

parse_event_date <- function(x) {
  parsed <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC"))
  parsed <- dplyr::if_else(is.na(parsed), suppressWarnings(lubridate::ymd(x)), parsed)
  as.Date(parsed)
}

build_clumps <- function(patient_daily_data, tau_days) {
  seizure_days <- patient_daily_data %>%
    dplyr::filter(.data$seizure_count > 0) %>%
    dplyr::arrange(.data$date) %>%
    dplyr::select("date", "seizure_count")

  if (nrow(seizure_days) == 0) {
    return(tibble::tibble(
      tau_days = integer(),
      clump_id = integer(),
      clump_start_date = as.Date(character()),
      clump_end_date = as.Date(character()),
      clump_duration_days = integer(),
      clump_seizure_days = integer(),
      clump_events = integer()
    ))
  }

  starts_new_clump <- c(TRUE, as.integer(diff(seizure_days$date)) > tau_days)

  seizure_days %>%
    dplyr::mutate(clump_id = cumsum(starts_new_clump)) %>%
    dplyr::group_by(.data$clump_id) %>%
    dplyr::summarise(
      tau_days = tau_days,
      clump_start_date = min(.data$date),
      clump_end_date = max(.data$date),
      clump_duration_days = as.integer(.data$clump_end_date - .data$clump_start_date + 1L),
      clump_seizure_days = dplyr::n(),
      clump_events = sum(.data$seizure_count),
      .groups = "drop"
    ) %>%
    dplyr::relocate("tau_days", "clump_id")
}

clump_summary_vector <- function(seizure_count, dates, tau_days) {
  seizure_count <- as.integer(seizure_count)
  seizure_day_index <- which(seizure_count > 0L)
  observed_days <- length(seizure_count)
  total_seizure_events <- sum(seizure_count)
  total_seizure_days <- length(seizure_day_index)

  if (total_seizure_days == 0L) {
    return(c(
      tau_days = tau_days,
      observed_days = observed_days,
      observed_seizure_days = 0,
      observed_seizure_events = 0,
      n_clumps = 0,
      clumps_per_28_days = 0,
      n_multi_event_clumps = 0,
      n_multi_day_clumps = 0,
      max_clump_events = 0,
      median_clump_events = NA_real_,
      mean_clump_events = NA_real_,
      max_clump_seizure_days = 0,
      max_clump_duration_days = 0,
      seizures_in_multi_event_clumps = 0,
      fraction_events_in_multi_event_clumps = NA_real_,
      seizure_days_in_multi_day_clumps = 0,
      fraction_seizure_days_in_multi_day_clumps = NA_real_,
      max_event_concentration = NA_real_,
      median_inter_clump_start_gap_days = NA_real_,
      mean_inter_clump_start_gap_days = NA_real_
    ))
  }

  seizure_dates <- dates[seizure_day_index]
  new_clump_starts <- c(1L, which(as.integer(diff(seizure_dates)) > tau_days) + 1L)
  clump_ends <- c(new_clump_starts[-1L] - 1L, length(seizure_day_index))
  n_clumps <- length(new_clump_starts)

  clump_events <- vapply(seq_len(n_clumps), function(clump_index) {
    start <- new_clump_starts[[clump_index]]
    end <- clump_ends[[clump_index]]
    sum(seizure_count[seizure_day_index[start:end]])
  }, numeric(1))
  clump_seizure_days <- clump_ends - new_clump_starts + 1L
  clump_duration_days <- as.integer(seizure_dates[clump_ends] - seizure_dates[new_clump_starts] + 1L)
  inter_clump_start_gaps <- as.integer(diff(seizure_dates[new_clump_starts]))
  seizures_in_multi_event_clumps <- sum(clump_events[clump_events >= 2L])
  seizure_days_in_multi_day_clumps <- sum(clump_seizure_days[clump_seizure_days >= 2L])

  c(
    tau_days = tau_days,
    observed_days = observed_days,
    observed_seizure_days = total_seizure_days,
    observed_seizure_events = total_seizure_events,
    n_clumps = n_clumps,
    clumps_per_28_days = ifelse(observed_days > 0, STANDARD_MONTH_DAYS * n_clumps / observed_days, NA_real_),
    n_multi_event_clumps = sum(clump_events >= 2L),
    n_multi_day_clumps = sum(clump_seizure_days >= 2L),
    max_clump_events = max(clump_events),
    median_clump_events = stats::median(clump_events),
    mean_clump_events = mean(clump_events),
    max_clump_seizure_days = max(clump_seizure_days),
    max_clump_duration_days = max(clump_duration_days),
    seizures_in_multi_event_clumps = seizures_in_multi_event_clumps,
    fraction_events_in_multi_event_clumps = ifelse(
      total_seizure_events > 0,
      seizures_in_multi_event_clumps / total_seizure_events,
      NA_real_
    ),
    seizure_days_in_multi_day_clumps = seizure_days_in_multi_day_clumps,
    fraction_seizure_days_in_multi_day_clumps = ifelse(
      total_seizure_days > 0,
      seizure_days_in_multi_day_clumps / total_seizure_days,
      NA_real_
    ),
    max_event_concentration = ifelse(
      total_seizure_events > 0,
      max(clump_events) / total_seizure_events,
      NA_real_
    ),
    median_inter_clump_start_gap_days = ifelse(
      length(inter_clump_start_gaps) > 0,
      stats::median(inter_clump_start_gaps),
      NA_real_
    ),
    mean_inter_clump_start_gap_days = ifelse(
      length(inter_clump_start_gaps) > 0,
      mean(inter_clump_start_gaps),
      NA_real_
    )
  )
}

clump_statistics <- function(patient_daily_data, tau_days) {
  tibble::as_tibble_row(clump_summary_vector(
    seizure_count = patient_daily_data$seizure_count,
    dates = patient_daily_data$date,
    tau_days = tau_days
  ))
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

summarise_null_metric <- function(null_values, observed_value) {
  valid_values <- null_values[is.finite(null_values)]
  valid_count <- length(valid_values)

  tibble::tibble(
    valid_permutation_replicates = valid_count,
    null_mean = if (valid_count > 0) mean(valid_values) else NA_real_,
    null_lower_95 = if (valid_count > 0) as.numeric(stats::quantile(valid_values, probs = 0.025, names = FALSE)) else NA_real_,
    null_upper_95 = if (valid_count > 0) as.numeric(stats::quantile(valid_values, probs = 0.975, names = FALSE)) else NA_real_,
    p_value = if (valid_count > 0 && is.finite(observed_value)) {
      (1 + sum(valid_values >= observed_value)) / (1 + valid_count)
    } else {
      NA_real_
    }
  )
}

patient_episodicity_test <- function(patient_daily_data) {
  observed <- purrr::map_dfr(TAU_DAYS, ~ clump_statistics(patient_daily_data, .x))

  base_summary <- observed %>%
    dplyr::mutate(
      n_segments = dplyr::n_distinct(patient_daily_data$segment_id),
      segment_sources = paste(sort(unique(patient_daily_data$segment_source)), collapse = ";"),
      primary_tau = .data$tau_days == PRIMARY_TAU_DAYS,
      permutation_replicates = PERMUTATION_REPLICATES
    )

  if (sum(patient_daily_data$seizure_count) == 0) {
    return(base_summary %>%
      dplyr::mutate(
        valid_permutation_replicates = 0L,
        null_max_clump_events_mean = NA_real_,
        null_max_clump_events_lower_95 = NA_real_,
        null_max_clump_events_upper_95 = NA_real_,
        p_value_max_clump_events = NA_real_,
        null_fraction_events_in_multi_event_clumps_mean = NA_real_,
        null_fraction_events_in_multi_event_clumps_lower_95 = NA_real_,
        null_fraction_events_in_multi_event_clumps_upper_95 = NA_real_,
        p_value_fraction_events_in_multi_event_clumps = NA_real_,
        null_max_clump_duration_days_mean = NA_real_,
        null_max_clump_duration_days_lower_95 = NA_real_,
        null_max_clump_duration_days_upper_95 = NA_real_,
        p_value_max_clump_duration_days = NA_real_,
        analysis_status = "not_run_no_seizure_events"
      ))
  }

  if (sum(patient_daily_data$seizure_count) < 2L) {
    return(base_summary %>%
      dplyr::mutate(
        valid_permutation_replicates = 0L,
        null_max_clump_events_mean = NA_real_,
        null_max_clump_events_lower_95 = NA_real_,
        null_max_clump_events_upper_95 = NA_real_,
        p_value_max_clump_events = NA_real_,
        null_fraction_events_in_multi_event_clumps_mean = NA_real_,
        null_fraction_events_in_multi_event_clumps_lower_95 = NA_real_,
        null_fraction_events_in_multi_event_clumps_upper_95 = NA_real_,
        p_value_fraction_events_in_multi_event_clumps = NA_real_,
        null_max_clump_duration_days_mean = NA_real_,
        null_max_clump_duration_days_lower_95 = NA_real_,
        null_max_clump_duration_days_upper_95 = NA_real_,
        p_value_max_clump_duration_days = NA_real_,
        analysis_status = "not_run_fewer_than_two_seizure_events"
      ))
  }

  null_max_clump_events <- matrix(NA_real_, nrow = PERMUTATION_REPLICATES, ncol = length(TAU_DAYS))
  null_fraction_events <- matrix(NA_real_, nrow = PERMUTATION_REPLICATES, ncol = length(TAU_DAYS))
  null_max_duration <- matrix(NA_real_, nrow = PERMUTATION_REPLICATES, ncol = length(TAU_DAYS))

  for (replicate_id in seq_len(PERMUTATION_REPLICATES)) {
    permuted_counts <- permute_counts_within_segments(
      patient_daily_data$seizure_count,
      patient_daily_data$segment_id
    )

    for (tau_index in seq_along(TAU_DAYS)) {
      summary_vector <- clump_summary_vector(
        seizure_count = permuted_counts,
        dates = patient_daily_data$date,
        tau_days = TAU_DAYS[[tau_index]]
      )
      null_max_clump_events[replicate_id, tau_index] <- summary_vector[["max_clump_events"]]
      null_fraction_events[replicate_id, tau_index] <- summary_vector[["fraction_events_in_multi_event_clumps"]]
      null_max_duration[replicate_id, tau_index] <- summary_vector[["max_clump_duration_days"]]
    }
  }

  permutation_summary <- purrr::map_dfr(TAU_DAYS, function(tau_days_value) {
    tau_index <- match(tau_days_value, TAU_DAYS)
    observed_row <- observed %>%
      dplyr::filter(.data$tau_days == tau_days_value)

    max_clump_events_summary <- summarise_null_metric(
      null_max_clump_events[, tau_index],
      observed_row$max_clump_events[[1]]
    )
    fraction_events_summary <- summarise_null_metric(
      null_fraction_events[, tau_index],
      observed_row$fraction_events_in_multi_event_clumps[[1]]
    )
    max_duration_summary <- summarise_null_metric(
      null_max_duration[, tau_index],
      observed_row$max_clump_duration_days[[1]]
    )

    tibble::tibble(
      tau_days = tau_days_value,
      valid_permutation_replicates = max_clump_events_summary$valid_permutation_replicates[[1]],
      null_max_clump_events_mean = max_clump_events_summary$null_mean[[1]],
      null_max_clump_events_lower_95 = max_clump_events_summary$null_lower_95[[1]],
      null_max_clump_events_upper_95 = max_clump_events_summary$null_upper_95[[1]],
      p_value_max_clump_events = max_clump_events_summary$p_value[[1]],
      null_fraction_events_in_multi_event_clumps_mean = fraction_events_summary$null_mean[[1]],
      null_fraction_events_in_multi_event_clumps_lower_95 = fraction_events_summary$null_lower_95[[1]],
      null_fraction_events_in_multi_event_clumps_upper_95 = fraction_events_summary$null_upper_95[[1]],
      p_value_fraction_events_in_multi_event_clumps = fraction_events_summary$p_value[[1]],
      null_max_clump_duration_days_mean = max_duration_summary$null_mean[[1]],
      null_max_clump_duration_days_lower_95 = max_duration_summary$null_lower_95[[1]],
      null_max_clump_duration_days_upper_95 = max_duration_summary$null_upper_95[[1]],
      p_value_max_clump_duration_days = max_duration_summary$p_value[[1]]
    )
  })

  base_summary %>%
    dplyr::left_join(permutation_summary, by = "tau_days") %>%
    dplyr::mutate(analysis_status = "ok")
}

patient_windows <- readr::read_csv(PANEL_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    study_start_date = as.Date(.data$study_start_date),
    study_end_date = pmin(as.Date(.data$study_end_date), ANALYSIS_CUTOFF_DATE)
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
    .data$event_date <= ANALYSIS_CUTOFF_DATE,
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
    segment_end_date = pmin(as.Date(.data$segment_end_date), ANALYSIS_CUTOFF_DATE),
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
    patient_episodicity_test(patient_data)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(patient_windows, by = "patient_id") %>%
  dplyr::mutate(
    q_value_max_clump_events = NA_real_,
    q_value_fraction_events_in_multi_event_clumps = NA_real_,
    q_value_max_clump_duration_days = NA_real_,
    significant_max_clump_events = FALSE,
    significant_fraction_events_in_multi_event_clumps = FALSE,
    significant_max_clump_duration_days = FALSE
  ) %>%
  dplyr::relocate(
    "patient_id",
    "variant_p",
    "study_start_date",
    "study_end_date",
    "tau_days",
    "primary_tau",
    "observed_days",
    "observed_seizure_days",
    "observed_seizure_events",
    "n_segments",
    "segment_sources"
  )

valid_max_event_p_values <- is.finite(patient_results$p_value_max_clump_events)
patient_results$q_value_max_clump_events[valid_max_event_p_values] <- stats::p.adjust(
  patient_results$p_value_max_clump_events[valid_max_event_p_values],
  method = "BH"
)

valid_fraction_p_values <- is.finite(patient_results$p_value_fraction_events_in_multi_event_clumps)
patient_results$q_value_fraction_events_in_multi_event_clumps[valid_fraction_p_values] <- stats::p.adjust(
  patient_results$p_value_fraction_events_in_multi_event_clumps[valid_fraction_p_values],
  method = "BH"
)

valid_duration_p_values <- is.finite(patient_results$p_value_max_clump_duration_days)
patient_results$q_value_max_clump_duration_days[valid_duration_p_values] <- stats::p.adjust(
  patient_results$p_value_max_clump_duration_days[valid_duration_p_values],
  method = "BH"
)

patient_results <- patient_results %>%
  dplyr::mutate(
    significant_max_clump_events = !is.na(.data$q_value_max_clump_events) & .data$q_value_max_clump_events < ALPHA,
    significant_fraction_events_in_multi_event_clumps = !is.na(.data$q_value_fraction_events_in_multi_event_clumps) &
      .data$q_value_fraction_events_in_multi_event_clumps < ALPHA,
    significant_max_clump_duration_days = !is.na(.data$q_value_max_clump_duration_days) &
      .data$q_value_max_clump_duration_days < ALPHA,
    significant_any_episodicity_metric = .data$significant_max_clump_events |
      .data$significant_fraction_events_in_multi_event_clumps |
      .data$significant_max_clump_duration_days
  )

observed_clumps <- daily_panel_segmented %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::group_modify(function(patient_data, patient_key) {
    purrr::map_dfr(TAU_DAYS, ~ build_clumps(patient_data, .x))
  }) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(patient_windows, by = "patient_id") %>%
  dplyr::relocate(
    "patient_id",
    "variant_p",
    "study_start_date",
    "study_end_date",
    "tau_days"
  )

readr::write_csv(daily_panel_segmented, DAILY_PANEL_OUTPUT_PATH)
readr::write_csv(patient_results, PATIENT_RESULTS_OUTPUT_PATH)
readr::write_csv(observed_clumps, OBSERVED_CLUMPS_OUTPUT_PATH)
