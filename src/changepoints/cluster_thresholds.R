suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(purrr)
  library(readr)
  library(tidyr)
})

source("src/analysis_config.R")

PANEL_INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
SEIZURE_INPUT_PATH <- "output/tabs/seizures/seizures.csv"
ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
SEGMENTS_INPUT_PATH <- "output/tabs/changepoints/change_point_segments.csv"
CHANGEPOINTS_INPUT_PATH <- "output/tabs/changepoints/patient_change_points.csv"
OUTPUT_DIR <- "output/tabs/changepoints"
OUTPUT_PATH <- file.path(OUTPUT_DIR, "cluster_thresholds.csv")

PANEL_REQUIRED_COLUMNS <- c("patient_id", "study_start_date", "study_end_date", "seizure_count")
SEIZURE_REQUIRED_COLUMNS <- c("patient_id", "date")
ASSIGNMENT_REQUIRED_COLUMNS <- c("patient_id", "pam_k3")
SEGMENT_REQUIRED_COLUMNS <- c("patient_id", "theta")
CHANGEPOINT_REQUIRED_COLUMNS <- c(
  "patient_id",
  "significant",
  "pre_observed_days",
  "post_observed_days",
  "pre_seizure_count",
  "post_seizure_count"
)

MIN_TOTAL_SEIZURES <- as.integer(Sys.getenv("CHANGEPOINT_MIN_TOTAL_SEIZURES", "1"))
MIN_SEGMENT_WEEKS <- as.integer(Sys.getenv("CHANGEPOINT_MIN_SEGMENT_WEEKS", "4"))
MIN_SEGMENT_OBSERVED_DAYS <- as.integer(Sys.getenv("CHANGEPOINT_MIN_SEGMENT_OBSERVED_DAYS", "28"))
SIMULATION_REPLICATES <- as.integer(Sys.getenv("CLUSTER_CHANGE_THRESHOLD_SIM_REPLICATES", "2000"))
THRESHOLD_PROBABILITY <- as.numeric(Sys.getenv("CLUSTER_CHANGE_THRESHOLD_PROBABILITY", "0.95"))
POISSON_LIMIT_THETA <- as.numeric(Sys.getenv("CHANGEPOINT_POISSON_LIMIT_THETA", "1000000"))
RANDOM_SEED <- 20260623L

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

parse_event_date <- function(x) {
  parsed <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC"))
  parsed <- dplyr::if_else(is.na(parsed), suppressWarnings(lubridate::ymd(x)), parsed)
  as.Date(parsed)
}

quantile_available <- function(x, probability) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }

  as.numeric(stats::quantile(x, probs = probability, names = FALSE, type = 7))
}

median_available <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }

  stats::median(x)
}

build_weekly_counts <- function(patient_window, seizure_events_by_week) {
  patient_id_value <- patient_window$patient_id[[1]]
  study_start_date <- patient_window$study_start_date[[1]]
  study_end_date <- patient_window$study_end_date[[1]]

  weekly_grid <- tibble::tibble(
    week = seq(
      lubridate::floor_date(study_start_date, unit = "week", week_start = 1),
      lubridate::floor_date(study_end_date, unit = "week", week_start = 1),
      by = "week"
    )
  ) %>%
    mutate(
      patient_id = patient_id_value,
      study_start_date = study_start_date,
      study_end_date = study_end_date,
      week_end = .data$week + 6L,
      observed_start = as.Date(pmax(.data$week, .data$study_start_date), origin = "1970-01-01"),
      observed_end = as.Date(pmin(.data$week_end, .data$study_end_date), origin = "1970-01-01"),
      observed_days = pmax(as.integer(.data$observed_end - .data$observed_start + 1L), 0L)
    )

  weekly_grid %>%
    left_join(
      seizure_events_by_week %>%
        filter(.data$patient_id == patient_id_value) %>%
        select("patient_id", "week", "seizure_count"),
      by = c("patient_id", "week")
    ) %>%
    tidyr::replace_na(list(seizure_count = 0L)) %>%
    filter(.data$observed_days > 0) %>%
    arrange(.data$week) %>%
    mutate(
      week_index = row_number(),
      seizure_count = as.integer(.data$seizure_count)
    ) %>%
    select(
      "patient_id",
      "week",
      "week_index",
      "seizure_count",
      "observed_days",
      "observed_start",
      "observed_end",
      "study_start_date",
      "study_end_date"
    )
}

valid_split_starts <- function(observed_days) {
  n <- length(observed_days)
  if (n < 2L * MIN_SEGMENT_WEEKS) {
    return(integer())
  }

  split_starts <- 2:n
  cumulative_days <- cumsum(observed_days)
  total_days <- sum(observed_days)

  split_starts[purrr::map_lgl(split_starts, function(split_start) {
    pre_weeks <- split_start - 1L
    post_weeks <- n - split_start + 1L
    pre_days <- cumulative_days[[split_start - 1L]]
    post_days <- total_days - pre_days

    pre_weeks >= MIN_SEGMENT_WEEKS &&
      post_weeks >= MIN_SEGMENT_WEEKS &&
      pre_days >= MIN_SEGMENT_OBSERVED_DAYS &&
      post_days >= MIN_SEGMENT_OBSERVED_DAYS
  })]
}

scan_split_metrics <- function(seizure_count, observed_days, split_starts) {
  if (length(split_starts) == 0) {
    return(tibble::tibble(
      max_absolute_change_per_28d = NA_real_,
      max_finite_fold_change = NA_real_,
      any_zero_transition = FALSE,
      max_zero_transition_abs_change_per_28d = 0
    ))
  }

  cumulative_counts <- cumsum(seizure_count)
  cumulative_days <- cumsum(observed_days)
  total_counts <- sum(seizure_count)
  total_days <- sum(observed_days)

  pre_counts <- cumulative_counts[split_starts - 1L]
  pre_days <- cumulative_days[split_starts - 1L]
  post_counts <- total_counts - pre_counts
  post_days <- total_days - pre_days

  pre_rates <- STANDARD_MONTH_DAYS * pre_counts / pre_days
  post_rates <- STANDARD_MONTH_DAYS * post_counts / post_days
  absolute_changes <- abs(post_rates - pre_rates)

  finite_fold_changes <- ifelse(
    pre_rates > 0 & post_rates > 0,
    pmax(post_rates / pre_rates, pre_rates / post_rates),
    NA_real_
  )
  zero_transitions <- (pre_rates == 0 & post_rates > 0) | (pre_rates > 0 & post_rates == 0)
  zero_transition_changes <- ifelse(zero_transitions, absolute_changes, 0)

  tibble::tibble(
    max_absolute_change_per_28d = max(absolute_changes, na.rm = TRUE),
    max_finite_fold_change = if (all(is.na(finite_fold_changes))) {
      NA_real_
    } else {
      max(finite_fold_changes, na.rm = TRUE)
    },
    any_zero_transition = any(zero_transitions, na.rm = TRUE),
    max_zero_transition_abs_change_per_28d = max(zero_transition_changes, na.rm = TRUE)
  )
}

scan_simulated_split_metrics <- function(simulated_counts, observed_days, split_starts) {
  if (length(split_starts) == 0 || nrow(simulated_counts) == 0) {
    return(tibble::tibble(
      max_absolute_change_per_28d = numeric(),
      max_finite_fold_change = numeric(),
      any_zero_transition = logical(),
      max_zero_transition_abs_change_per_28d = numeric()
    ))
  }

  cumulative_counts <- t(apply(simulated_counts, 1, cumsum))
  cumulative_days <- cumsum(observed_days)
  total_counts <- rowSums(simulated_counts)
  total_days <- sum(observed_days)

  pre_counts <- cumulative_counts[, split_starts - 1L, drop = FALSE]
  post_counts <- total_counts - pre_counts
  pre_days <- cumulative_days[split_starts - 1L]
  post_days <- total_days - pre_days

  pre_rates <- sweep(pre_counts, 2, pre_days, "/") * STANDARD_MONTH_DAYS
  post_rates <- sweep(post_counts, 2, post_days, "/") * STANDARD_MONTH_DAYS
  absolute_changes <- abs(post_rates - pre_rates)

  finite_fold_changes <- ifelse(
    pre_rates > 0 & post_rates > 0,
    pmax(post_rates / pre_rates, pre_rates / post_rates),
    NA_real_
  )
  zero_transitions <- (pre_rates == 0 & post_rates > 0) | (pre_rates > 0 & post_rates == 0)
  zero_transition_changes <- ifelse(zero_transitions, absolute_changes, 0)

  tibble::tibble(
    max_absolute_change_per_28d = apply(absolute_changes, 1, max, na.rm = TRUE),
    max_finite_fold_change = apply(finite_fold_changes, 1, function(x) {
      if (all(is.na(x))) {
        return(NA_real_)
      }

      max(x, na.rm = TRUE)
    }),
    any_zero_transition = apply(zero_transitions, 1, any, na.rm = TRUE),
    max_zero_transition_abs_change_per_28d = apply(zero_transition_changes, 1, max, na.rm = TRUE)
  )
}

simulate_patient_null <- function(patient_weekly_data, theta, n_replicates, patient_id_value) {
  observed_days <- as.numeric(patient_weekly_data$observed_days)
  observed_counts <- as.integer(patient_weekly_data$seizure_count)
  split_starts <- valid_split_starts(observed_days)
  total_days <- sum(observed_days)
  total_seizures <- sum(observed_counts)

  if (
    length(split_starts) == 0 ||
      total_days <= 0 ||
      total_seizures < MIN_TOTAL_SEIZURES ||
      !is.finite(theta) ||
      theta <= 0
  ) {
    return(tibble::tibble(
      replicate = integer(),
      max_absolute_change_per_28d = numeric(),
      max_finite_fold_change = numeric(),
      any_zero_transition = logical(),
      max_zero_transition_abs_change_per_28d = numeric()
    ))
  }

  stable_rate_per_day <- total_seizures / total_days
  mu <- stable_rate_per_day * observed_days
  theta <- min(theta, POISSON_LIMIT_THETA)

  simulated_counts <- matrix(
    stats::rnbinom(
      n = n_replicates * length(observed_days),
      size = theta,
      mu = rep(mu, times = n_replicates)
    ),
    nrow = n_replicates,
    byrow = TRUE
  )

  metrics <- scan_simulated_split_metrics(simulated_counts, observed_days, split_starts) %>%
    mutate(replicate = row_number(), .before = 1)

  metrics
}

panel <- check_input(PANEL_INPUT_PATH, PANEL_REQUIRED_COLUMNS, "Patient-month panel")
seizures <- check_input(SEIZURE_INPUT_PATH, SEIZURE_REQUIRED_COLUMNS, "Seizure event")
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

patient_windows <- panel %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    study_start_date = as.Date(.data$study_start_date),
    study_end_date = as.Date(.data$study_end_date),
    seizure_count = suppressWarnings(as.integer(as.numeric(.data$seizure_count)))
  ) %>%
  filter(
    !is.na(.data$patient_id),
    !is.na(.data$study_start_date),
    !is.na(.data$study_end_date),
    !is.na(.data$seizure_count)
  ) %>%
  group_by(.data$patient_id) %>%
  summarise(
    study_start_date = min(.data$study_start_date),
    study_end_date = max(.data$study_end_date),
    panel_total_seizure_events = sum(.data$seizure_count),
    .groups = "drop"
  ) %>%
  filter(.data$study_start_date <= .data$study_end_date)

seizure_events_by_week <- seizures %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    event_date = parse_event_date(.data$date)
  ) %>%
  filter(!is.na(.data$patient_id), !is.na(.data$event_date)) %>%
  inner_join(patient_windows, by = "patient_id") %>%
  filter(
    .data$event_date >= .data$study_start_date,
    .data$event_date <= .data$study_end_date
  ) %>%
  mutate(week = as.Date(lubridate::floor_date(.data$event_date, unit = "week", week_start = 1))) %>%
  count(.data$patient_id, .data$week, name = "seizure_count")

weekly_counts <- patient_windows %>%
  group_by(.data$patient_id) %>%
  group_split() %>%
  purrr::map_dfr(~ build_weekly_counts(.x, seizure_events_by_week)) %>%
  inner_join(cluster_lookup, by = "patient_id")

theta_lookup <- segments %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    theta = suppressWarnings(as.numeric(.data$theta))
  ) %>%
  filter(!is.na(.data$patient_id), is.finite(.data$theta), .data$theta > 0) %>%
  group_by(.data$patient_id) %>%
  summarise(theta = first(.data$theta), .groups = "drop")

patient_simulation_inputs <- weekly_counts %>%
  group_by(.data$patient_id, .data$pam_k3) %>%
  summarise(
    weekly_data = list(dplyr::pick("seizure_count", "observed_days")),
    total_observed_days = sum(.data$observed_days),
    total_seizures = sum(.data$seizure_count),
    total_weeks = n(),
    n_valid_split_starts = length(valid_split_starts(.data$observed_days)),
    baseline_rate_per_28d = STANDARD_MONTH_DAYS * .data$total_seizures / .data$total_observed_days,
    .groups = "drop"
  ) %>%
  inner_join(theta_lookup, by = "patient_id") %>%
  filter(
    .data$total_seizures >= MIN_TOTAL_SEIZURES,
    .data$n_valid_split_starts > 0
  )

if (nrow(patient_simulation_inputs) == 0) {
  stop("No patients were eligible for cluster change-threshold simulation.", call. = FALSE)
}

set.seed(RANDOM_SEED)

simulation_results <- patient_simulation_inputs %>%
  mutate(
    simulation = purrr::pmap(
      list(.data$weekly_data, .data$theta, .data$patient_id),
      ~ simulate_patient_null(..1, ..2, SIMULATION_REPLICATES, ..3)
    )
  ) %>%
  select("patient_id", "pam_k3", "simulation") %>%
  tidyr::unnest(cols = c("simulation"))

cluster_patient_counts <- cluster_lookup %>%
  count(.data$pam_k3, name = "n_cluster_patients")

simulation_thresholds <- simulation_results %>%
  group_by(.data$pam_k3) %>%
  summarise(
    n_simulation_patients = n_distinct(.data$patient_id),
    simulation_replicates_per_patient = SIMULATION_REPLICATES,
    n_patient_replicates = n(),
    median_baseline_rate_per_28d = median_available(
      patient_simulation_inputs$baseline_rate_per_28d[
        patient_simulation_inputs$pam_k3 == first(.data$pam_k3)
      ]
    ),
    q95_max_absolute_change_per_28d_under_null = quantile_available(
      .data$max_absolute_change_per_28d,
      THRESHOLD_PROBABILITY
    ),
    q95_max_finite_fold_change_under_null = quantile_available(
      .data$max_finite_fold_change,
      THRESHOLD_PROBABILITY
    ),
    n_replicates_with_finite_fold_change = sum(!is.na(.data$max_finite_fold_change)),
    null_probability_any_zero_transition = mean(.data$any_zero_transition, na.rm = TRUE),
    q95_max_zero_transition_abs_change_per_28d_under_null = quantile_available(
      .data$max_zero_transition_abs_change_per_28d,
      THRESHOLD_PROBABILITY
    ),
    .groups = "drop"
  )

observed_significant_changes <- change_points %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    significant = .data$significant %in% TRUE | as.character(.data$significant) == "TRUE",
    pre_observed_days = suppressWarnings(as.numeric(.data$pre_observed_days)),
    post_observed_days = suppressWarnings(as.numeric(.data$post_observed_days)),
    pre_seizure_count = suppressWarnings(as.numeric(.data$pre_seizure_count)),
    post_seizure_count = suppressWarnings(as.numeric(.data$post_seizure_count))
  ) %>%
  inner_join(cluster_lookup, by = "patient_id") %>%
  filter(.data$significant) %>%
  mutate(
    pre_rate_per_28d = STANDARD_MONTH_DAYS * .data$pre_seizure_count / .data$pre_observed_days,
    post_rate_per_28d = STANDARD_MONTH_DAYS * .data$post_seizure_count / .data$post_observed_days,
    absolute_change_per_28d = abs(.data$post_rate_per_28d - .data$pre_rate_per_28d),
    finite_fold_change = if_else(
      .data$pre_rate_per_28d > 0 & .data$post_rate_per_28d > 0,
      pmax(
        .data$post_rate_per_28d / .data$pre_rate_per_28d,
        .data$pre_rate_per_28d / .data$post_rate_per_28d
      ),
      NA_real_
    ),
    zero_transition = (.data$pre_rate_per_28d == 0 & .data$post_rate_per_28d > 0) |
      (.data$pre_rate_per_28d > 0 & .data$post_rate_per_28d == 0)
  )

observed_summary <- observed_significant_changes %>%
  group_by(.data$pam_k3) %>%
  summarise(
    observed_significant_changepoints = n(),
    observed_median_significant_abs_change_per_28d = median_available(.data$absolute_change_per_28d),
    observed_min_significant_abs_change_per_28d = min(.data$absolute_change_per_28d, na.rm = TRUE),
    observed_median_finite_fold_change = median_available(.data$finite_fold_change),
    observed_min_finite_fold_change = if (all(is.na(.data$finite_fold_change))) {
      NA_real_
    } else {
      min(.data$finite_fold_change, na.rm = TRUE)
    },
    observed_zero_transition_changepoints = sum(.data$zero_transition, na.rm = TRUE),
    .groups = "drop"
  )

threshold_table <- cluster_patient_counts %>%
  full_join(simulation_thresholds, by = "pam_k3") %>%
  full_join(observed_summary, by = "pam_k3") %>%
  mutate(
    across(
      all_of(c(
        "n_simulation_patients",
        "n_patient_replicates",
        "n_replicates_with_finite_fold_change",
        "observed_significant_changepoints",
        "observed_zero_transition_changepoints"
      )),
      ~ tidyr::replace_na(.x, 0)
    ),
    simulation_replicates_per_patient = tidyr::replace_na(
      .data$simulation_replicates_per_patient,
      SIMULATION_REPLICATES
    ),
    recommended_absolute_change_threshold_per_28d = .data$q95_max_absolute_change_per_28d_under_null,
    recommended_finite_fold_change_threshold = .data$q95_max_finite_fold_change_under_null,
    recommended_zero_transition_abs_change_threshold_per_28d =
      .data$q95_max_zero_transition_abs_change_per_28d_under_null
  ) %>%
  transmute(
    pam_k3 = .data$pam_k3,
    n_cluster_patients = as.integer(.data$n_cluster_patients),
    n_simulation_patients = as.integer(.data$n_simulation_patients),
    simulation_replicates_per_patient = as.integer(.data$simulation_replicates_per_patient),
    n_patient_replicates = as.integer(.data$n_patient_replicates),
    median_baseline_rate_per_28d = round(.data$median_baseline_rate_per_28d, 2),
    recommended_absolute_change_threshold_per_28d = round(
      .data$recommended_absolute_change_threshold_per_28d,
      2
    ),
    recommended_finite_fold_change_threshold = round(
      .data$recommended_finite_fold_change_threshold,
      2
    ),
    n_replicates_with_finite_fold_change = as.integer(.data$n_replicates_with_finite_fold_change),
    null_probability_any_zero_transition = round(.data$null_probability_any_zero_transition, 3),
    recommended_zero_transition_abs_change_threshold_per_28d = round(
      .data$recommended_zero_transition_abs_change_threshold_per_28d,
      2
    ),
    observed_significant_changepoints = as.integer(.data$observed_significant_changepoints),
    observed_median_significant_abs_change_per_28d = round(
      .data$observed_median_significant_abs_change_per_28d,
      2
    ),
    observed_min_significant_abs_change_per_28d = round(
      .data$observed_min_significant_abs_change_per_28d,
      2
    ),
    observed_median_finite_fold_change = round(.data$observed_median_finite_fold_change, 2),
    observed_min_finite_fold_change = round(.data$observed_min_finite_fold_change, 2),
    observed_zero_transition_changepoints = as.integer(.data$observed_zero_transition_changepoints)
  ) %>%
  arrange(.data$pam_k3)

readr::write_csv(threshold_table, OUTPUT_PATH)
