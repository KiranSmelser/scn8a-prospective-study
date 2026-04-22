# Check whether daily seizure counts are consistent with a Poisson process.

suppressPackageStartupMessages({
  library(tidyverse)
})

INPUT_PATH <- "output/tabs/seizures/seizures.csv"
OUTPUT_DIR <- "output/tabs/seizures"
OVERALL_OUTPUT_PATH <- file.path(OUTPUT_DIR, "poisson_fit_overall.csv")
PATIENT_OUTPUT_PATH <- file.path(OUTPUT_DIR, "poisson_fit_by_patient.csv")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

group_bins_for_chisq <- function(observed, expected) {
  groups <- list()
  current_observed <- 0
  current_expected <- 0
  start_index <- 1

  for (i in seq_along(expected)) {
    current_observed <- current_observed + observed[i]
    current_expected <- current_expected + expected[i]

    if (current_expected >= 5) {
      groups[[length(groups) + 1]] <- list(
        start = start_index,
        end = i,
        observed = current_observed,
        expected = current_expected
      )
      start_index <- i + 1
      current_observed <- 0
      current_expected <- 0
    }
  }

  if (current_expected > 0) {
    if (length(groups) == 0) {
      groups[[1]] <- list(
        start = 1,
        end = length(expected),
        observed = current_observed,
        expected = current_expected
      )
    } else {
      last <- groups[[length(groups)]]
      groups[[length(groups)]] <- list(
        start = last$start,
        end = length(expected),
        observed = last$observed + current_observed,
        expected = last$expected + current_expected
      )
    }
  }

  tibble::tibble(
    observed = purrr::map_dbl(groups, "observed"),
    expected = purrr::map_dbl(groups, "expected")
  )
}

poisson_gof <- function(counts) {
  counts <- counts[is.finite(counts)]
  counts <- as.integer(counts[counts >= 0])
  n <- length(counts)

  if (n == 0) {
    return(tibble::tibble(
      n = 0L,
      lambda_hat = NA_real_,
      variance = NA_real_,
      variance_to_mean = NA_real_,
      chisq_statistic = NA_real_,
      chisq_df = NA_integer_,
      chisq_p_value = NA_real_,
      poisson_fit = NA_character_
    ))
  }

  lambda_hat <- mean(counts)
  variance <- stats::var(counts)
  variance_to_mean <- ifelse(lambda_hat > 0, variance / lambda_hat, NA_real_)

  if (lambda_hat == 0) {
    return(tibble::tibble(
      n = n,
      lambda_hat = lambda_hat,
      variance = variance,
      variance_to_mean = variance_to_mean,
      chisq_statistic = NA_real_,
      chisq_df = NA_integer_,
      chisq_p_value = NA_real_,
      poisson_fit = "All counts are zero; Poisson GOF test is not informative"
    ))
  }

  max_k <- max(max(counts), stats::qpois(0.999, lambda_hat))
  categories <- 0:max_k
  observed <- as.numeric(table(factor(ifelse(counts > max_k, max_k + 1L, counts), levels = 0:(max_k + 1L))))
  expected_probs <- c(stats::dpois(categories, lambda_hat), 1 - stats::ppois(max_k, lambda_hat))
  expected <- n * expected_probs

  grouped <- group_bins_for_chisq(observed = observed, expected = expected)
  chisq_df <- nrow(grouped) - 2L

  if (chisq_df < 1 || any(grouped$expected <= 0)) {
    chisq_statistic <- NA_real_
    chisq_p_value <- NA_real_
    poisson_fit <- "Insufficient grouped bins for valid chi-squared GOF"
  } else {
    chisq_statistic <- sum((grouped$observed - grouped$expected)^2 / grouped$expected)
    chisq_p_value <- stats::pchisq(chisq_statistic, df = chisq_df, lower.tail = FALSE)
    poisson_fit <- ifelse(chisq_p_value < 0.05, "No (reject at alpha=0.05)", "Yes (fail to reject at alpha=0.05)")
  }

  tibble::tibble(
    n = n,
    lambda_hat = lambda_hat,
    variance = variance,
    variance_to_mean = variance_to_mean,
    chisq_statistic = chisq_statistic,
    chisq_df = chisq_df,
    chisq_p_value = chisq_p_value,
    poisson_fit = poisson_fit
  )
}

if (!file.exists(INPUT_PATH)) {
  stop("Input file not found: ", INPUT_PATH)
}

seizures <- readr::read_csv(INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    day = as.Date(.data$date)
  ) %>%
  dplyr::filter(!is.na(.data$patient_id), !is.na(.data$day))

daily_counts <- seizures %>%
  dplyr::count(.data$patient_id, .data$day, name = "seizure_count") %>%
  dplyr::group_by(.data$patient_id) %>%
  tidyr::complete(
    day = seq(min(.data$day), max(.data$day), by = "day"),
    fill = list(seizure_count = 0L)
  ) %>%
  dplyr::ungroup()

overall_result <- poisson_gof(daily_counts$seizure_count) %>%
  dplyr::mutate(level = "overall_patient_day")

patient_results <- daily_counts %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(stats = list(poisson_gof(.data$seizure_count)), .groups = "drop") %>%
  tidyr::unnest(cols = c(stats)) %>%
  dplyr::arrange(.data$chisq_p_value)

readr::write_csv(overall_result, OVERALL_OUTPUT_PATH)
readr::write_csv(patient_results, PATIENT_OUTPUT_PATH)
