# Change-point analysis for patient seizure frequency.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(purrr)
  library(MASS)
})

PANEL_INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
SEIZURE_INPUT_PATH <- "output/tabs/seizures/seizures.csv"
OUTPUT_DIR <- "output/tabs/changepoints"
CANDIDATES_OUTPUT_PATH <- file.path(OUTPUT_DIR, "change_point_candidates.csv")
PATIENT_OUTPUT_PATH <- file.path(OUTPUT_DIR, "patient_change_points.csv")
SEGMENTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "change_point_segments.csv")
PANEL_REQUIRED_COLUMNS <- c("patient_id", "study_start_date", "study_end_date", "month", "seizure_count")
SEIZURE_REQUIRED_COLUMNS <- c("patient_id", "date")

MIN_TOTAL_SEIZURES <- 1L
MIN_SEGMENT_WEEKS <- 4L
MIN_SEGMENT_OBSERVED_DAYS <- 28L
BOOTSTRAP_REPLICATES <- as.integer(Sys.getenv("CHANGEPOINT_BOOTSTRAPS", "500"))
RANDOM_SEED <- 20260430L
ALPHA <- 0.05

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PANEL_INPUT_PATH)) {
  stop("Input file not found: ", PANEL_INPUT_PATH)
}

if (!file.exists(SEIZURE_INPUT_PATH)) {
  stop("Input file not found: ", SEIZURE_INPUT_PATH)
}

parse_event_date <- function(x) {
  parsed <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC"))
  parsed <- dplyr::if_else(is.na(parsed), suppressWarnings(lubridate::ymd(x)), parsed)
  as.Date(parsed)
}

safe_glm_nb <- function(formula, data) {
  tryCatch(
    suppressWarnings(MASS::glm.nb(formula = formula, data = data, control = glm.control(maxit = 100))),
    error = function(...) NULL
  )
}

safe_glm_nb_fixed_theta <- function(formula, data, theta) {
  if (!is.finite(theta) || theta <= 0) {
    return(NULL)
  }

  tryCatch(
    suppressWarnings(stats::glm(
      formula = formula,
      data = data,
      family = MASS::negative.binomial(theta = theta),
      control = glm.control(maxit = 100)
    )),
    error = function(...) NULL
  )
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
    dplyr::mutate(
      patient_id = patient_id_value,
      study_start_date = study_start_date,
      study_end_date = study_end_date,
      week_end = .data$week + 6L,
      observed_start = as.Date(pmax(.data$week, .data$study_start_date), origin = "1970-01-01"),
      observed_end = as.Date(pmin(.data$week_end, .data$study_end_date), origin = "1970-01-01"),
      observed_days = pmax(as.integer(.data$observed_end - .data$observed_start + 1L), 0L)
    )

  weekly_grid %>%
    dplyr::left_join(
      seizure_events_by_week %>%
        dplyr::filter(.data$patient_id == patient_id_value) %>%
        dplyr::select("patient_id", "week", "seizure_count"),
      by = c("patient_id", "week")
    ) %>%
    tidyr::replace_na(list(seizure_count = 0L)) %>%
    dplyr::filter(.data$observed_days > 0) %>%
    dplyr::arrange(.data$week) %>%
    dplyr::mutate(
      week_index = dplyr::row_number(),
      seizure_count = as.integer(.data$seizure_count)
    ) %>%
    dplyr::select(
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

summarise_segment <- function(data, segment_label) {
  total_seizures <- sum(data$seizure_count)
  total_days <- sum(data$observed_days)

  tibble::tibble(
    segment = segment_label,
    segment_start_date = if (nrow(data) == 0) as.Date(NA) else min(data$observed_start, na.rm = TRUE),
    segment_end_date = if (nrow(data) == 0) as.Date(NA) else max(data$observed_end, na.rm = TRUE),
    weeks = nrow(data),
    observed_days = total_days,
    seizure_count = total_seizures,
    rate_per_30_days = ifelse(total_days > 0, 30 * total_seizures / total_days, NA_real_)
  )
}

candidate_summary <- function(weekly_data, candidate_week) {
  pre_data <- weekly_data %>% dplyr::filter(.data$week < candidate_week)
  post_data <- weekly_data %>% dplyr::filter(.data$week >= candidate_week)

  pre <- summarise_segment(pre_data, "pre")
  post <- summarise_segment(post_data, "post")

  tibble::tibble(
    patient_id = weekly_data$patient_id[[1]],
    candidate_week = candidate_week,
    pre_segment_start_date = pre$segment_start_date,
    pre_segment_end_date = pre$segment_end_date,
    post_segment_start_date = post$segment_start_date,
    post_segment_end_date = post$segment_end_date,
    pre_weeks = pre$weeks,
    post_weeks = post$weeks,
    pre_observed_days = pre$observed_days,
    post_observed_days = post$observed_days,
    pre_seizure_count = pre$seizure_count,
    post_seizure_count = post$seizure_count,
    pre_rate_per_30_days = pre$rate_per_30_days,
    post_rate_per_30_days = post$rate_per_30_days,
    rate_ratio_post_vs_pre = ifelse(
      is.finite(pre$rate_per_30_days) && pre$rate_per_30_days > 0,
      post$rate_per_30_days / pre$rate_per_30_days,
      NA_real_
    ),
    direction = dplyr::case_when(
      is.na(.data$rate_ratio_post_vs_pre) ~ NA_character_,
      .data$rate_ratio_post_vs_pre > 1 ~ "increase",
      .data$rate_ratio_post_vs_pre < 1 ~ "decrease",
      TRUE ~ "no_change"
    )
  )
}

candidate_weeks_for_patient <- function(weekly_data) {
  weekly_data %>%
    dplyr::filter(
      .data$week_index > MIN_SEGMENT_WEEKS,
      .data$week_index <= max(.data$week_index) - MIN_SEGMENT_WEEKS + 1L
    ) %>%
    dplyr::pull(.data$week) %>%
    purrr::keep(function(candidate_week) {
      pre_days <- weekly_data %>%
        dplyr::filter(.data$week < candidate_week) %>%
        dplyr::summarise(observed_days = sum(.data$observed_days), .groups = "drop") %>%
        dplyr::pull(.data$observed_days)
      post_days <- weekly_data %>%
        dplyr::filter(.data$week >= candidate_week) %>%
        dplyr::summarise(observed_days = sum(.data$observed_days), .groups = "drop") %>%
        dplyr::pull(.data$observed_days)

      pre_days >= MIN_SEGMENT_OBSERVED_DAYS && post_days >= MIN_SEGMENT_OBSERVED_DAYS
    })
}

scan_candidate <- function(weekly_data, candidate_week, null_log_likelihood, theta) {
  model_data <- weekly_data %>%
    dplyr::mutate(
      post_change = as.integer(.data$week >= candidate_week),
      log_observed_days = log(.data$observed_days)
    )

  alternative_model <- safe_glm_nb_fixed_theta(
    seizure_count ~ post_change + offset(log_observed_days),
    data = model_data,
    theta = theta
  )

  candidate_stats <- candidate_summary(weekly_data, candidate_week)

  if (!is.finite(null_log_likelihood) || is.null(alternative_model)) {
    return(
      candidate_stats %>%
        dplyr::mutate(
          lrt_statistic = NA_real_,
          raw_p_value = NA_real_,
          model_status = dplyr::if_else(
            is.finite(null_log_likelihood),
            "alternative_negative_binomial_fit_failed",
            "null_negative_binomial_fit_failed"
          )
        )
    )
  }

  lrt_statistic <- max(0, 2 * (as.numeric(stats::logLik(alternative_model)) - null_log_likelihood))
  raw_p_value <- stats::pchisq(lrt_statistic, df = 1, lower.tail = FALSE)

  candidate_stats %>%
    dplyr::mutate(
      lrt_statistic = lrt_statistic,
      raw_p_value = raw_p_value,
      model_status = "ok"
    )
}

scan_patient <- function(weekly_data) {
  candidate_weeks <- candidate_weeks_for_patient(weekly_data)

  if (length(candidate_weeks) == 0) {
    return(tibble::tibble(
      patient_id = weekly_data$patient_id[[1]],
      candidate_week = as.Date(character()),
      pre_segment_start_date = as.Date(character()),
      pre_segment_end_date = as.Date(character()),
      post_segment_start_date = as.Date(character()),
      post_segment_end_date = as.Date(character()),
      pre_weeks = integer(),
      post_weeks = integer(),
      pre_observed_days = integer(),
      post_observed_days = integer(),
      pre_seizure_count = integer(),
      post_seizure_count = integer(),
      pre_rate_per_30_days = numeric(),
      post_rate_per_30_days = numeric(),
      rate_ratio_post_vs_pre = numeric(),
      direction = character(),
      lrt_statistic = numeric(),
      raw_p_value = numeric(),
      model_status = character()
    ))
  }

  model_data <- weekly_data %>%
    dplyr::mutate(log_observed_days = log(.data$observed_days))

  null_model <- safe_glm_nb(
    seizure_count ~ 1 + offset(log_observed_days),
    data = model_data
  )
  theta <- if (is.null(null_model)) NA_real_ else null_model$theta
  null_fixed_model <- safe_glm_nb_fixed_theta(
    seizure_count ~ 1 + offset(log_observed_days),
    data = model_data,
    theta = theta
  )
  null_log_likelihood <- if (is.null(null_fixed_model)) NA_real_ else as.numeric(stats::logLik(null_fixed_model))

  purrr::map_dfr(candidate_weeks, ~ scan_candidate(weekly_data, .x, null_log_likelihood, theta))
}

max_lrt_for_counts <- function(seizure_count, base_weekly_data, candidate_weeks, theta) {
  simulated_data <- base_weekly_data %>%
    dplyr::mutate(
      seizure_count = as.integer(.env$seizure_count),
      log_observed_days = log(.data$observed_days)
    )

  null_model <- safe_glm_nb_fixed_theta(
    seizure_count ~ 1 + offset(log_observed_days),
    data = simulated_data,
    theta = theta
  )

  if (is.null(null_model)) {
    return(NA_real_)
  }

  null_log_likelihood <- as.numeric(stats::logLik(null_model))
  lrt_statistics <- purrr::map_dbl(candidate_weeks, function(candidate_week) {
    candidate_data <- simulated_data %>%
      dplyr::mutate(post_change = as.integer(.data$week >= candidate_week))

    alternative_model <- safe_glm_nb_fixed_theta(
      seizure_count ~ post_change + offset(log_observed_days),
      data = candidate_data,
      theta = theta
    )

    if (is.null(alternative_model)) {
      return(NA_real_)
    }

    max(0, 2 * (as.numeric(stats::logLik(alternative_model)) - null_log_likelihood))
  })

  if (all(!is.finite(lrt_statistics))) {
    return(NA_real_)
  }

  max(lrt_statistics, na.rm = TRUE)
}

bootstrap_patient_p_value <- function(weekly_data, candidate_results, n_bootstrap) {
  valid_candidates <- candidate_results %>%
    dplyr::filter(.data$model_status == "ok", is.finite(.data$lrt_statistic))

  if (nrow(valid_candidates) == 0 || n_bootstrap <= 0) {
    return(tibble::tibble(
      observed_max_lrt = NA_real_,
      bootstrap_p_value = NA_real_,
      bootstrap_replicates = n_bootstrap,
      bootstrap_valid_replicates = 0L,
      bootstrap_status = "not_run_no_valid_candidates"
    ))
  }

  model_data <- weekly_data %>%
    dplyr::mutate(log_observed_days = log(.data$observed_days))

  null_model <- safe_glm_nb(
    seizure_count ~ 1 + offset(log_observed_days),
    data = model_data
  )

  if (is.null(null_model)) {
    return(tibble::tibble(
      observed_max_lrt = max(valid_candidates$lrt_statistic, na.rm = TRUE),
      bootstrap_p_value = NA_real_,
      bootstrap_replicates = n_bootstrap,
      bootstrap_valid_replicates = 0L,
      bootstrap_status = "not_run_null_negative_binomial_fit_failed"
    ))
  }

  observed_max_lrt <- max(valid_candidates$lrt_statistic, na.rm = TRUE)
  candidate_weeks <- valid_candidates$candidate_week
  mu <- as.numeric(stats::fitted(null_model))
  theta <- null_model$theta

  simulated_max_lrt <- replicate(
    n_bootstrap,
    {
      simulated_counts <- stats::rnbinom(n = length(mu), mu = mu, size = theta)
      max_lrt_for_counts(simulated_counts, weekly_data, candidate_weeks, theta)
    }
  )

  simulated_max_lrt <- simulated_max_lrt[is.finite(simulated_max_lrt)]
  valid_replicates <- length(simulated_max_lrt)

  if (valid_replicates == 0) {
    return(tibble::tibble(
      observed_max_lrt = observed_max_lrt,
      bootstrap_p_value = NA_real_,
      bootstrap_replicates = n_bootstrap,
      bootstrap_valid_replicates = valid_replicates,
      bootstrap_status = "failed_all_bootstrap_replicates"
    ))
  }

  tibble::tibble(
    observed_max_lrt = observed_max_lrt,
    bootstrap_p_value = (sum(simulated_max_lrt >= observed_max_lrt) + 1) / (valid_replicates + 1),
    bootstrap_replicates = n_bootstrap,
    bootstrap_valid_replicates = valid_replicates,
    bootstrap_status = "ok"
  )
}

select_patient_change_point <- function(patient_summary, candidate_results, bootstrap_result) {
  if (nrow(candidate_results) == 0) {
    return(patient_summary %>%
      dplyr::mutate(
        candidate_week = as.Date(NA),
        pre_segment_start_date = as.Date(NA),
        pre_segment_end_date = as.Date(NA),
        post_segment_start_date = as.Date(NA),
        post_segment_end_date = as.Date(NA),
        pre_weeks = NA_integer_,
        post_weeks = NA_integer_,
        pre_observed_days = NA_integer_,
        post_observed_days = NA_integer_,
        pre_seizure_count = NA_integer_,
        post_seizure_count = NA_integer_,
        pre_rate_per_30_days = NA_real_,
        post_rate_per_30_days = NA_real_,
        rate_ratio_post_vs_pre = NA_real_,
        direction = NA_character_,
        lrt_statistic = NA_real_,
        raw_p_value = NA_real_,
        observed_max_lrt = NA_real_,
        bootstrap_p_value = NA_real_,
        bootstrap_replicates = BOOTSTRAP_REPLICATES,
        bootstrap_valid_replicates = 0L,
        bootstrap_status = "not_run_no_candidate_weeks",
        q_value = NA_real_,
        significant = FALSE,
        analysis_status = "insufficient_timeline_for_candidate_scan"
      ))
  }

  best_candidate <- candidate_results %>%
    dplyr::arrange(dplyr::desc(.data$lrt_statistic), .data$candidate_week) %>%
    dplyr::slice(1)

  patient_summary %>%
    dplyr::bind_cols(best_candidate %>% dplyr::select(-dplyr::all_of("patient_id"))) %>%
    dplyr::bind_cols(bootstrap_result) %>%
    dplyr::mutate(
      q_value = NA_real_,
      significant = FALSE,
      analysis_status = dplyr::if_else(.data$model_status == "ok", "ok", .data$model_status)
    )
}

patient_month_panel <- readr::read_csv(PANEL_INPUT_PATH, show_col_types = FALSE)
seizures <- readr::read_csv(SEIZURE_INPUT_PATH, show_col_types = FALSE)

missing_panel_columns <- setdiff(PANEL_REQUIRED_COLUMNS, names(patient_month_panel))
if (length(missing_panel_columns) > 0) {
  stop("Patient-month panel is missing required columns: ", paste(missing_panel_columns, collapse = ", "))
}

missing_seizure_columns <- setdiff(SEIZURE_REQUIRED_COLUMNS, names(seizures))
if (length(missing_seizure_columns) > 0) {
  stop("Seizure event data is missing required columns: ", paste(missing_seizure_columns, collapse = ", "))
}

patient_windows <- patient_month_panel %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    study_start_date = as.Date(.data$study_start_date),
    study_end_date = as.Date(.data$study_end_date),
    seizure_count = as.integer(suppressWarnings(as.numeric(.data$seizure_count)))
  ) %>%
  dplyr::filter(
    !is.na(.data$patient_id),
    !is.na(.data$study_start_date),
    !is.na(.data$study_end_date),
    !is.na(.data$seizure_count)
  ) %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(
    study_start_date = min(.data$study_start_date),
    study_end_date = max(.data$study_end_date),
    panel_total_seizure_events = sum(.data$seizure_count),
    .groups = "drop"
  ) %>%
  dplyr::filter(.data$study_start_date <= .data$study_end_date)

if (nrow(patient_windows) == 0) {
  stop("No valid patient study windows found in input file: ", PANEL_INPUT_PATH)
}

seizure_events_by_week <- seizures %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    event_date = parse_event_date(.data$date)
  ) %>%
  dplyr::filter(!is.na(.data$patient_id), !is.na(.data$event_date)) %>%
  dplyr::inner_join(patient_windows, by = "patient_id") %>%
  dplyr::filter(
    .data$event_date >= .data$study_start_date,
    .data$event_date <= .data$study_end_date
  ) %>%
  dplyr::mutate(week = as.Date(lubridate::floor_date(.data$event_date, unit = "week", week_start = 1))) %>%
  dplyr::count(.data$patient_id, .data$week, name = "seizure_count")

weekly_counts_all <- patient_windows %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(~ build_weekly_counts(.x, seizure_events_by_week))

if (nrow(weekly_counts_all) == 0) {
  stop("No valid patient-week rows could be built from input file: ", PANEL_INPUT_PATH)
}

eligible_patients <- weekly_counts_all %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(total_seizure_events = sum(.data$seizure_count), .groups = "drop") %>%
  dplyr::filter(.data$total_seizure_events >= MIN_TOTAL_SEIZURES)

if (nrow(eligible_patients) == 0) {
  stop("No patients met the minimum seizure threshold of ", MIN_TOTAL_SEIZURES, " event(s).")
}

weekly_counts <- weekly_counts_all %>%
  dplyr::semi_join(eligible_patients, by = "patient_id") %>%
  dplyr::arrange(.data$patient_id, .data$week)

patient_summaries <- weekly_counts %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(
    total_weeks = dplyr::n(),
    total_observed_days = sum(.data$observed_days),
    total_seizure_events = sum(.data$seizure_count),
    study_start_date = min(.data$study_start_date),
    study_end_date = max(.data$study_end_date),
    first_observed_date = min(.data$observed_start),
    last_observed_date = max(.data$observed_end),
    mean_rate_per_30_days = 30 * .data$total_seizure_events / .data$total_observed_days,
    .groups = "drop"
  )

set.seed(RANDOM_SEED)

patient_results <- weekly_counts %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::group_split() %>%
  purrr::map(function(patient_weekly_data) {
    candidate_results <- scan_patient(patient_weekly_data)
    bootstrap_result <- bootstrap_patient_p_value(
      weekly_data = patient_weekly_data,
      candidate_results = candidate_results,
      n_bootstrap = BOOTSTRAP_REPLICATES
    )
    patient_summary <- patient_summaries %>%
      dplyr::filter(.data$patient_id == patient_weekly_data$patient_id[[1]])

    list(
      candidates = candidate_results,
      patient = select_patient_change_point(patient_summary, candidate_results, bootstrap_result)
    )
  })

candidate_table <- patient_results %>%
  purrr::map("candidates") %>%
  purrr::list_rbind() %>%
  dplyr::arrange(.data$patient_id, .data$candidate_week)

patient_table <- patient_results %>%
  purrr::map("patient") %>%
  purrr::list_rbind()

if (any(is.finite(patient_table$bootstrap_p_value))) {
  finite_bootstrap_p <- is.finite(patient_table$bootstrap_p_value)
  patient_table$q_value[finite_bootstrap_p] <- stats::p.adjust(
    patient_table$bootstrap_p_value[finite_bootstrap_p],
    method = "BH"
  )

  patient_table <- patient_table %>%
    dplyr::mutate(
      significant = !is.na(.data$q_value) & .data$q_value < ALPHA
    )
}

segment_table <- patient_table %>%
  dplyr::filter(!is.na(.data$candidate_week)) %>%
  dplyr::select(
    "patient_id",
    "candidate_week",
    "pre_segment_start_date",
    "pre_segment_end_date",
    "post_segment_start_date",
    "post_segment_end_date",
    "pre_weeks",
    "post_weeks",
    "pre_observed_days",
    "post_observed_days",
    "pre_seizure_count",
    "post_seizure_count",
    "pre_rate_per_30_days",
    "post_rate_per_30_days"
  ) %>%
  tidyr::pivot_longer(
    cols = -c("patient_id", "candidate_week"),
    names_to = c("segment", ".value"),
    names_pattern = "^(pre|post)_(.*)$"
  ) %>%
  dplyr::arrange(.data$patient_id, .data$candidate_week, .data$segment)

readr::write_csv(candidate_table, CANDIDATES_OUTPUT_PATH)
readr::write_csv(patient_table, PATIENT_OUTPUT_PATH)
readr::write_csv(segment_table, SEGMENTS_OUTPUT_PATH)

message("Wrote candidate-level results to: ", CANDIDATES_OUTPUT_PATH)
message("Wrote patient-level results to: ", PATIENT_OUTPUT_PATH)
message("Wrote segment summaries to: ", SEGMENTS_OUTPUT_PATH)
