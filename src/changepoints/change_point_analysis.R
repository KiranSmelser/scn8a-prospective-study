# Joint multiple-change-point analysis for patient seizure frequency.

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
LEGACY_BINARY_OUTPUT_PATHS <- file.path(
  OUTPUT_DIR,
  c("binary_segmentation_change_points.csv", "binary_segmentation_segments.csv")
)
PANEL_REQUIRED_COLUMNS <- c("patient_id", "study_start_date", "study_end_date", "month", "seizure_count")
SEIZURE_REQUIRED_COLUMNS <- c("patient_id", "date")

MIN_TOTAL_SEIZURES <- as.integer(Sys.getenv("CHANGEPOINT_MIN_TOTAL_SEIZURES", "1"))
MIN_SEGMENT_WEEKS <- as.integer(Sys.getenv("CHANGEPOINT_MIN_SEGMENT_WEEKS", "4"))
MIN_SEGMENT_OBSERVED_DAYS <- as.integer(Sys.getenv("CHANGEPOINT_MIN_SEGMENT_OBSERVED_DAYS", "28"))
MAX_CHANGEPOINTS <- as.integer(Sys.getenv("CHANGEPOINT_MAX_CHANGEPOINTS", "4"))
BOOTSTRAP_REPLICATES <- as.integer(Sys.getenv("CHANGEPOINT_BOOTSTRAPS", "500"))
SEGMENT_PENALTY_MULTIPLIER <- as.numeric(Sys.getenv("CHANGEPOINT_SEGMENT_PENALTY_MULTIPLIER", "2"))
POISSON_LIMIT_THETA <- as.numeric(Sys.getenv("CHANGEPOINT_POISSON_LIMIT_THETA", "1000000"))
RANDOM_SEED <- 20260430L
ALPHA <- 0.05

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
invisible(file.remove(LEGACY_BINARY_OUTPUT_PATHS[file.exists(LEGACY_BINARY_OUTPUT_PATHS)]))

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

rate_ratio_post_vs_pre <- function(pre_rate, post_rate) {
  dplyr::case_when(
    is.na(pre_rate) | is.na(post_rate) ~ NA_real_,
    pre_rate == 0 & post_rate == 0 ~ 1,
    pre_rate == 0 & post_rate > 0 ~ Inf,
    pre_rate > 0 ~ post_rate / pre_rate,
    TRUE ~ NA_real_
  )
}

rate_change_direction <- function(pre_rate, post_rate) {
  dplyr::case_when(
    is.na(pre_rate) | is.na(post_rate) ~ NA_character_,
    post_rate > pre_rate ~ "increase",
    post_rate < pre_rate ~ "decrease",
    TRUE ~ "no_change"
  )
}

estimate_patient_theta <- function(weekly_data) {
  model_data <- weekly_data %>%
    dplyr::mutate(log_observed_days = log(.data$observed_days))

  null_model <- safe_glm_nb(
    seizure_count ~ 1 + offset(log_observed_days),
    data = model_data
  )

  if (is.null(null_model) || !is.finite(null_model$theta) || null_model$theta <= 0) {
    return(tibble::tibble(
      theta = POISSON_LIMIT_THETA,
      theta_status = "negative_binomial_theta_fit_failed_poisson_limit"
    ))
  }

  tibble::tibble(
    theta = as.numeric(null_model$theta),
    theta_status = "ok"
  )
}

segment_log_likelihood <- function(seizure_count, observed_days, theta) {
  seizure_count <- as.integer(seizure_count)
  observed_days <- as.numeric(observed_days)

  if (!is.finite(theta) || theta <= 0 || any(!is.finite(observed_days)) || any(observed_days <= 0)) {
    return(NA_real_)
  }

  total_seizures <- sum(seizure_count)
  total_days <- sum(observed_days)
  if (total_days <= 0) {
    return(NA_real_)
  }

  if (total_seizures == 0) {
    return(0)
  }

  rate_per_day <- total_seizures / total_days
  mu <- rate_per_day * observed_days
  sum(stats::dnbinom(seizure_count, size = theta, mu = mu, log = TRUE))
}

candidate_split_indices <- function(weekly_data) {
  n <- nrow(weekly_data)
  if (n < 2 * MIN_SEGMENT_WEEKS) {
    return(integer())
  }

  starts <- 2:n
  starts[purrr::map_lgl(starts, function(split_start) {
    pre_data <- weekly_data[seq_len(split_start - 1L), , drop = FALSE]
    post_data <- weekly_data[split_start:n, , drop = FALSE]

    nrow(pre_data) >= MIN_SEGMENT_WEEKS &&
      nrow(post_data) >= MIN_SEGMENT_WEEKS &&
      sum(pre_data$observed_days) >= MIN_SEGMENT_OBSERVED_DAYS &&
      sum(post_data$observed_days) >= MIN_SEGMENT_OBSERVED_DAYS
  })]
}

candidate_summary <- function(weekly_data, split_start, theta, null_log_likelihood) {
  n <- nrow(weekly_data)
  pre_data <- weekly_data[seq_len(split_start - 1L), , drop = FALSE]
  post_data <- weekly_data[split_start:n, , drop = FALSE]
  pre <- summarise_segment(pre_data, "pre")
  post <- summarise_segment(post_data, "post")
  pre_log_likelihood <- segment_log_likelihood(pre_data$seizure_count, pre_data$observed_days, theta)
  post_log_likelihood <- segment_log_likelihood(post_data$seizure_count, post_data$observed_days, theta)
  alternative_log_likelihood <- pre_log_likelihood + post_log_likelihood
  lrt_statistic <- ifelse(
    is.finite(null_log_likelihood) && is.finite(alternative_log_likelihood),
    max(0, 2 * (alternative_log_likelihood - null_log_likelihood)),
    NA_real_
  )

  tibble::tibble(
    patient_id = weekly_data$patient_id[[1]],
    candidate_week = weekly_data$week[[split_start]],
    candidate_week_index = split_start,
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
    rate_ratio_post_vs_pre = rate_ratio_post_vs_pre(pre$rate_per_30_days, post$rate_per_30_days),
    direction = rate_change_direction(pre$rate_per_30_days, post$rate_per_30_days),
    lrt_statistic = lrt_statistic,
    raw_p_value = stats::pchisq(lrt_statistic, df = 1, lower.tail = FALSE),
    model_status = ifelse(is.finite(lrt_statistic), "ok", "likelihood_failed"),
    theta = theta,
    total_weeks = n
  )
}

scan_patient_candidates <- function(weekly_data, theta) {
  split_indices <- candidate_split_indices(weekly_data)

  if (length(split_indices) == 0) {
    return(tibble::tibble(
      patient_id = weekly_data$patient_id[[1]],
      candidate_week = as.Date(character()),
      candidate_week_index = integer(),
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
      model_status = character(),
      theta = numeric(),
      total_weeks = integer()
    ))
  }

  null_log_likelihood <- segment_log_likelihood(weekly_data$seizure_count, weekly_data$observed_days, theta)
  purrr::map_dfr(split_indices, ~ candidate_summary(weekly_data, .x, theta, null_log_likelihood))
}

build_interval_costs <- function(weekly_data, theta) {
  n <- nrow(weekly_data)
  cost <- matrix(Inf, nrow = n, ncol = n)
  log_likelihood <- matrix(NA_real_, nrow = n, ncol = n)

  for (segment_start in seq_len(n)) {
    for (segment_end in segment_start:n) {
      segment_data <- weekly_data[segment_start:segment_end, , drop = FALSE]
      if (
        nrow(segment_data) >= MIN_SEGMENT_WEEKS &&
          sum(segment_data$observed_days) >= MIN_SEGMENT_OBSERVED_DAYS
      ) {
        ll <- segment_log_likelihood(segment_data$seizure_count, segment_data$observed_days, theta)
        if (is.finite(ll)) {
          log_likelihood[segment_start, segment_end] <- ll
          cost[segment_start, segment_end] <- -2 * ll
        }
      }
    }
  }

  list(cost = cost, log_likelihood = log_likelihood)
}

select_joint_segmentation <- function(weekly_data, theta) {
  n <- nrow(weekly_data)
  interval_costs <- build_interval_costs(weekly_data, theta)
  cost <- interval_costs$cost

  if (!is.finite(cost[1, n])) {
    return(list(
      status = "insufficient_timeline_for_joint_segmentation",
      starts = 1L,
      ends = n,
      n_segments = 1L,
      n_change_points = 0L,
      null_cost = NA_real_,
      selected_cost = NA_real_,
      selected_objective = NA_real_,
      penalty_per_segment = NA_real_,
      interval_costs = interval_costs
    ))
  }

  max_segments <- min(MAX_CHANGEPOINTS + 1L, n)
  penalty_per_segment <- SEGMENT_PENALTY_MULTIPLIER * log(n)
  dp <- matrix(Inf, nrow = max_segments, ncol = n)
  backtrack <- matrix(NA_integer_, nrow = max_segments, ncol = n)

  for (segment_end in seq_len(n)) {
    if (is.finite(cost[1, segment_end])) {
      dp[1, segment_end] <- cost[1, segment_end]
      backtrack[1, segment_end] <- 1L
    }
  }

  if (max_segments >= 2L) {
    for (k in 2:max_segments) {
      for (segment_end in seq_len(n)) {
        if (segment_end < 2L) {
          next
        }

        possible_starts <- 2:segment_end
        possible_starts <- possible_starts[
          is.finite(cost[possible_starts, segment_end]) &
            is.finite(dp[k - 1L, possible_starts - 1L])
        ]

        if (length(possible_starts) > 0) {
          objectives <- dp[k - 1L, possible_starts - 1L] + cost[possible_starts, segment_end]
          best_index <- which.min(objectives)
          dp[k, segment_end] <- objectives[[best_index]]
          backtrack[k, segment_end] <- possible_starts[[best_index]]
        }
      }
    }
  }

  available_segments <- which(is.finite(dp[, n]))
  objectives <- dp[available_segments, n] + available_segments * penalty_per_segment
  selected_k <- available_segments[[which.min(objectives)]]

  starts <- integer(selected_k)
  ends <- integer(selected_k)
  current_end <- n
  for (k in seq(from = selected_k, to = 1L)) {
    current_start <- backtrack[k, current_end]
    starts[[k]] <- current_start
    ends[[k]] <- current_end
    current_end <- current_start - 1L
  }

  list(
    status = "ok",
    starts = starts,
    ends = ends,
    n_segments = selected_k,
    n_change_points = selected_k - 1L,
    null_cost = dp[1, n],
    selected_cost = dp[selected_k, n],
    selected_objective = dp[selected_k, n] + selected_k * penalty_per_segment,
    penalty_per_segment = penalty_per_segment,
    interval_costs = interval_costs
  )
}

build_segment_table <- function(weekly_data, segmentation, theta, theta_status) {
  purrr::map2_dfr(segmentation$starts, segmentation$ends, function(segment_start, segment_end) {
    segment_data <- weekly_data[segment_start:segment_end, , drop = FALSE]
    segment_summary <- summarise_segment(segment_data, "selected")
    segment_log_likelihood_value <- segment_log_likelihood(segment_data$seizure_count, segment_data$observed_days, theta)

    tibble::tibble(
      patient_id = weekly_data$patient_id[[1]],
      segment_id = which(segmentation$starts == segment_start & segmentation$ends == segment_end)[[1]],
      segment_start_week_index = segment_start,
      segment_end_week_index = segment_end,
      segment_start_date = segment_summary$segment_start_date,
      segment_end_date = segment_summary$segment_end_date,
      weeks = segment_summary$weeks,
      observed_days = segment_summary$observed_days,
      seizure_count = segment_summary$seizure_count,
      rate_per_30_days = segment_summary$rate_per_30_days,
      log_likelihood = segment_log_likelihood_value,
      theta = theta,
      theta_status = theta_status,
      n_model_segments = segmentation$n_segments,
      n_model_change_points = segmentation$n_change_points,
      null_cost = segmentation$null_cost,
      selected_cost = segmentation$selected_cost,
      selected_objective = segmentation$selected_objective,
      penalty_per_segment = segmentation$penalty_per_segment,
      segmentation_status = segmentation$status
    )
  })
}

best_local_lrt_for_counts <- function(seizure_count, observed_days, theta) {
  n <- length(seizure_count)
  null_log_likelihood <- segment_log_likelihood(seizure_count, observed_days, theta)
  if (!is.finite(null_log_likelihood)) {
    return(NA_real_)
  }

  split_starts <- 2:n
  split_starts <- split_starts[purrr::map_lgl(split_starts, function(split_start) {
    pre_weeks <- split_start - 1L
    post_weeks <- n - split_start + 1L
    pre_days <- sum(observed_days[seq_len(split_start - 1L)])
    post_days <- sum(observed_days[split_start:n])

    pre_weeks >= MIN_SEGMENT_WEEKS &&
      post_weeks >= MIN_SEGMENT_WEEKS &&
      pre_days >= MIN_SEGMENT_OBSERVED_DAYS &&
      post_days >= MIN_SEGMENT_OBSERVED_DAYS
  })]

  if (length(split_starts) == 0) {
    return(NA_real_)
  }

  lrt_values <- purrr::map_dbl(split_starts, function(split_start) {
    pre_log_likelihood <- segment_log_likelihood(
      seizure_count[seq_len(split_start - 1L)],
      observed_days[seq_len(split_start - 1L)],
      theta
    )
    post_log_likelihood <- segment_log_likelihood(
      seizure_count[split_start:n],
      observed_days[split_start:n],
      theta
    )
    if (!is.finite(pre_log_likelihood) || !is.finite(post_log_likelihood)) {
      return(NA_real_)
    }

    max(0, 2 * (pre_log_likelihood + post_log_likelihood - null_log_likelihood))
  })

  if (all(!is.finite(lrt_values))) {
    return(NA_real_)
  }

  max(lrt_values, na.rm = TRUE)
}

bootstrap_local_p_value <- function(weekly_data, segment_start, segment_end, theta, observed_lrt, n_bootstrap) {
  if (!is.finite(observed_lrt) || n_bootstrap <= 0) {
    return(tibble::tibble(
      bootstrap_p_value = NA_real_,
      bootstrap_replicates = n_bootstrap,
      bootstrap_valid_replicates = 0L,
      bootstrap_status = "not_run"
    ))
  }

  merged_data <- weekly_data[segment_start:segment_end, , drop = FALSE]
  if (length(candidate_split_indices(merged_data)) == 0) {
    return(tibble::tibble(
      bootstrap_p_value = NA_real_,
      bootstrap_replicates = n_bootstrap,
      bootstrap_valid_replicates = 0L,
      bootstrap_status = "not_run_no_valid_local_splits"
    ))
  }

  total_seizures <- sum(merged_data$seizure_count)
  total_days <- sum(merged_data$observed_days)
  mu <- total_seizures * merged_data$observed_days / total_days

  simulated_lrt <- replicate(
    n_bootstrap,
    {
      simulated_counts <- stats::rnbinom(n = nrow(merged_data), size = theta, mu = mu)
      best_local_lrt_for_counts(simulated_counts, merged_data$observed_days, theta)
    }
  )
  simulated_lrt <- simulated_lrt[is.finite(simulated_lrt)]
  valid_replicates <- length(simulated_lrt)

  if (valid_replicates == 0) {
    return(tibble::tibble(
      bootstrap_p_value = NA_real_,
      bootstrap_replicates = n_bootstrap,
      bootstrap_valid_replicates = 0L,
      bootstrap_status = "failed_all_bootstrap_replicates"
    ))
  }

  tibble::tibble(
    bootstrap_p_value = (sum(simulated_lrt >= observed_lrt) + 1) / (valid_replicates + 1),
    bootstrap_replicates = n_bootstrap,
    bootstrap_valid_replicates = valid_replicates,
    bootstrap_status = "ok"
  )
}

empty_patient_change_point_row <- function(patient_summary, theta, theta_status, segmentation) {
  patient_summary %>%
    dplyr::mutate(
      candidate_week = as.Date(NA),
      candidate_week_index = NA_integer_,
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
      bootstrap_p_value = NA_real_,
      bootstrap_replicates = BOOTSTRAP_REPLICATES,
      bootstrap_valid_replicates = 0L,
      bootstrap_status = "not_run_no_selected_change_point",
      q_value = NA_real_,
      significant = FALSE,
      theta = theta,
      theta_status = theta_status,
      n_model_segments = segmentation$n_segments,
      n_model_change_points = segmentation$n_change_points,
      null_cost = segmentation$null_cost,
      selected_cost = segmentation$selected_cost,
      selected_objective = segmentation$selected_objective,
      penalty_per_segment = segmentation$penalty_per_segment,
      analysis_status = segmentation$status
    )
}

build_change_point_table <- function(weekly_data, patient_summary, segmentation, segment_table, theta, theta_status) {
  if (segmentation$n_change_points == 0L || nrow(segment_table) < 2L) {
    return(empty_patient_change_point_row(patient_summary, theta, theta_status, segmentation))
  }

  purrr::map_dfr(seq_len(segmentation$n_change_points), function(change_index) {
    pre_segment <- segment_table[change_index, , drop = FALSE]
    post_segment <- segment_table[change_index + 1L, , drop = FALSE]
    merged_start <- pre_segment$segment_start_week_index[[1]]
    merged_end <- post_segment$segment_end_week_index[[1]]
    merged_data <- weekly_data[merged_start:merged_end, , drop = FALSE]
    merged_log_likelihood <- segment_log_likelihood(merged_data$seizure_count, merged_data$observed_days, theta)
    split_log_likelihood <- pre_segment$log_likelihood[[1]] + post_segment$log_likelihood[[1]]
    lrt_statistic <- ifelse(
      is.finite(merged_log_likelihood) && is.finite(split_log_likelihood),
      max(0, 2 * (split_log_likelihood - merged_log_likelihood)),
      NA_real_
    )
    bootstrap_result <- bootstrap_local_p_value(
      weekly_data = weekly_data,
      segment_start = merged_start,
      segment_end = merged_end,
      theta = theta,
      observed_lrt = lrt_statistic,
      n_bootstrap = BOOTSTRAP_REPLICATES
    )
    pre_rate <- pre_segment$rate_per_30_days[[1]]
    post_rate <- post_segment$rate_per_30_days[[1]]
    rate_ratio <- rate_ratio_post_vs_pre(pre_rate, post_rate)

    patient_summary %>%
      dplyr::mutate(
        candidate_week = weekly_data$week[[post_segment$segment_start_week_index[[1]]]],
        candidate_week_index = post_segment$segment_start_week_index[[1]],
        pre_segment_start_date = pre_segment$segment_start_date[[1]],
        pre_segment_end_date = pre_segment$segment_end_date[[1]],
        post_segment_start_date = post_segment$segment_start_date[[1]],
        post_segment_end_date = post_segment$segment_end_date[[1]],
        pre_weeks = pre_segment$weeks[[1]],
        post_weeks = post_segment$weeks[[1]],
        pre_observed_days = pre_segment$observed_days[[1]],
        post_observed_days = post_segment$observed_days[[1]],
        pre_seizure_count = pre_segment$seizure_count[[1]],
        post_seizure_count = post_segment$seizure_count[[1]],
        pre_rate_per_30_days = pre_rate,
        post_rate_per_30_days = post_rate,
        rate_ratio_post_vs_pre = rate_ratio,
        direction = rate_change_direction(pre_rate, post_rate),
        lrt_statistic = lrt_statistic,
        raw_p_value = stats::pchisq(lrt_statistic, df = 1, lower.tail = FALSE),
        bootstrap_p_value = bootstrap_result$bootstrap_p_value[[1]],
        bootstrap_replicates = bootstrap_result$bootstrap_replicates[[1]],
        bootstrap_valid_replicates = bootstrap_result$bootstrap_valid_replicates[[1]],
        bootstrap_status = bootstrap_result$bootstrap_status[[1]],
        q_value = NA_real_,
        significant = FALSE,
        theta = theta,
        theta_status = theta_status,
        n_model_segments = segmentation$n_segments,
        n_model_change_points = segmentation$n_change_points,
        null_cost = segmentation$null_cost,
        selected_cost = segmentation$selected_cost,
        selected_objective = segmentation$selected_objective,
        penalty_per_segment = segmentation$penalty_per_segment,
        analysis_status = segmentation$status
      )
  })
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
    theta_result <- estimate_patient_theta(patient_weekly_data)
    theta <- theta_result$theta[[1]]
    theta_status <- theta_result$theta_status[[1]]
    segmentation <- select_joint_segmentation(patient_weekly_data, theta)
    patient_summary <- patient_summaries %>%
      dplyr::filter(.data$patient_id == patient_weekly_data$patient_id[[1]])
    segment_table <- build_segment_table(patient_weekly_data, segmentation, theta, theta_status)
    change_point_table <- build_change_point_table(
      weekly_data = patient_weekly_data,
      patient_summary = patient_summary,
      segmentation = segmentation,
      segment_table = segment_table,
      theta = theta,
      theta_status = theta_status
    )
    candidate_table <- scan_patient_candidates(patient_weekly_data, theta) %>%
      dplyr::mutate(
        selected_by_joint_model = .data$candidate_week_index %in% change_point_table$candidate_week_index,
        joint_model_n_segments = segmentation$n_segments,
        joint_model_n_change_points = segmentation$n_change_points,
        joint_model_status = segmentation$status,
        penalty_per_segment = segmentation$penalty_per_segment
      )

    list(
      candidates = candidate_table,
      change_points = change_point_table,
      segments = segment_table
    )
  })

candidate_table <- patient_results %>%
  purrr::map("candidates") %>%
  purrr::list_rbind() %>%
  dplyr::arrange(.data$patient_id, .data$candidate_week)

patient_table <- patient_results %>%
  purrr::map("change_points") %>%
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

segment_table <- patient_results %>%
  purrr::map("segments") %>%
  purrr::list_rbind() %>%
  dplyr::arrange(.data$patient_id, .data$segment_start_date, .data$segment_id)

readr::write_csv(candidate_table, CANDIDATES_OUTPUT_PATH)
readr::write_csv(patient_table, PATIENT_OUTPUT_PATH)
readr::write_csv(segment_table, SEGMENTS_OUTPUT_PATH)
