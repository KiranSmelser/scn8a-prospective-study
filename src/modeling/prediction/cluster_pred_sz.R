# Predict cluster membership from early seizure data.

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(lubridate)
  library(pROC)
  library(readr)
  library(tibble)
  library(tidyr)
})

SEIZURES_INPUT_PATH <- "output/tabs/seizures/seizures.csv"
PANEL_INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
CLUSTER_ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
OUTPUT_DIR <- "output/tabs/modeling/prediction"

FEATURES_OUTPUT_PATH <- file.path(OUTPUT_DIR, "early_seizure_features.csv")
PREDICTIONS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "ridge_multinomial_early_cluster_predictions.csv")
PERFORMANCE_OUTPUT_PATH <- file.path(OUTPUT_DIR, "ridge_multinomial_early_cluster_performance.csv")
CONFUSION_OUTPUT_PATH <- file.path(OUTPUT_DIR, "ridge_multinomial_early_cluster_confusion_matrix.csv")
AUC_OUTPUT_PATH <- file.path(OUTPUT_DIR, "ridge_multinomial_early_cluster_auc.csv")
COEFFICIENTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "ridge_multinomial_early_cluster_coefficients.csv")

HORIZON_WEEKS <- c(4L, 6L, 8L)
RANDOM_SEED <- 20260609L
RIDGE_ALPHA <- 0
INNER_CV_MAX_FOLDS <- 5L

MODEL_FEATURES <- c(
  "log1p_total_seizures",
  "log1p_seizure_rate_per_28_days",
  "log1p_mean_weekly_seizures",
  "log1p_median_weekly_seizures",
  "log1p_max_weekly_seizures",
  "log1p_sd_weekly_seizures",
  "log1p_iqr_weekly_seizures",
  "proportion_zero_seizure_weeks",
  "first_seizure_day_fraction",
  "last_seizure_day_fraction",
  "log1p_mean_interseizure_interval_days",
  "n_seizure_types",
  "log1p_mean_duration_seconds",
  "log1p_max_duration_seconds",
  "proportion_during_sleep_yes",
  "during_sleep_known_proportion"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(SEIZURES_INPUT_PATH, PANEL_INPUT_PATH, CLUSTER_ASSIGNMENTS_INPUT_PATH)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Required input file(s) not found: ", paste(missing_files, collapse = ", "), call. = FALSE)
}

parse_event_datetime <- function(x) {
  parsed <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC"))
  parsed <- dplyr::if_else(
    is.na(parsed),
    suppressWarnings(as.POSIXct(lubridate::ymd(x), tz = "UTC")),
    parsed
  )
  parsed
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }
  x[[1]]
}

safe_mean <- function(x, default = 0) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(default)
  }
  mean(x)
}

safe_max <- function(x, default = 0) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(default)
  }
  max(x)
}

safe_sd <- function(x, default = 0) {
  x <- x[is.finite(x)]
  if (length(x) < 2) {
    return(default)
  }
  stats::sd(x)
}

check_columns <- function(data, required_columns, label) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

make_stratified_foldid <- function(y, max_folds = INNER_CV_MAX_FOLDS) {
  class_counts <- table(y)
  nfolds <- min(max_folds, min(class_counts), length(y))

  if (nfolds < 2) {
    return(NULL)
  }

  foldid <- integer(length(y))
  for (class_level in names(class_counts)) {
    class_indices <- which(y == class_level)
    class_indices <- sample(class_indices, length(class_indices))
    foldid[class_indices] <- rep(seq_len(nfolds), length.out = length(class_indices))
  }

  foldid
}

fit_cv_ridge_multinomial <- function(x, y) {
  foldid <- make_stratified_foldid(y)

  if (is.null(foldid)) {
    stop("At least two observations per class are required for inner cross-validation.", call. = FALSE)
  }

  withCallingHandlers(
    glmnet::cv.glmnet(
      x = x,
      y = y,
      family = "multinomial",
      alpha = RIDGE_ALPHA,
      type.measure = "class",
      foldid = foldid,
      standardize = TRUE
    ),
    warning = function(w) {
      if (grepl("one multinomial or binomial class has fewer than 8", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

extract_probability_row <- function(prediction_array, cluster_levels) {
  prediction_dimensions <- dim(prediction_array)
  if (length(prediction_dimensions) != 3L) {
    stop("Unexpected glmnet multinomial prediction shape.", call. = FALSE)
  }

  predicted_classes <- dimnames(prediction_array)[[2]]
  probabilities <- setNames(rep(0, length(cluster_levels)), cluster_levels)
  probabilities[predicted_classes] <- as.numeric(prediction_array[1, predicted_classes, 1])
  probabilities
}

multiclass_log_loss <- function(truth, probabilities, cluster_levels) {
  eps <- 1e-15
  truth <- as.character(truth)
  row_index <- seq_along(truth)
  class_index <- match(truth, cluster_levels)
  selected_probabilities <- probabilities[cbind(row_index, class_index)]
  selected_probabilities <- pmin(pmax(selected_probabilities, eps), 1 - eps)
  -mean(log(selected_probabilities))
}

build_metrics <- function(data, cluster_levels) {
  truth <- factor(data$true_pam_k3, levels = cluster_levels)
  predicted <- factor(data$predicted_pam_k3, levels = cluster_levels)
  confusion <- table(truth, predicted)

  recalls <- diag(confusion) / rowSums(confusion)
  precisions <- diag(confusion) / colSums(confusion)
  f1 <- 2 * precisions * recalls / (precisions + recalls)
  f1[!is.finite(f1)] <- 0

  probabilities <- data %>%
    select(all_of(paste0("prob_cluster_", cluster_levels))) %>%
    as.matrix()

  tibble::tibble(
    n_patients = nrow(data),
    n_clusters = length(cluster_levels),
    min_cluster_n = min(as.integer(table(truth))),
    accuracy = mean(truth == predicted),
    balanced_accuracy = mean(recalls, na.rm = TRUE),
    macro_f1 = mean(f1),
    multiclass_log_loss = multiclass_log_loss(truth, probabilities, cluster_levels),
    majority_class = names(which.max(table(truth))),
    majority_class_accuracy = max(table(truth)) / length(truth)
  )
}

build_auc_table <- function(data, cluster_levels) {
  purrr::map_dfr(
    cluster_levels,
    function(cluster_level) {
      truth_binary <- as.integer(as.character(data$true_pam_k3) == cluster_level)
      probability <- data[[paste0("prob_cluster_", cluster_level)]]

      auc_value <- if (length(unique(truth_binary)) < 2) {
        NA_real_
      } else {
        as.numeric(pROC::auc(pROC::roc(truth_binary, probability, quiet = TRUE)))
      }

      tibble::tibble(
        class = cluster_level,
        one_vs_rest_auc = auc_value
      )
    }
  )
}

build_coefficient_table <- function(cv_fit, horizon_weeks, feature_names, cluster_levels) {
  coefficient_list <- stats::coef(cv_fit, s = "lambda.min")

  purrr::map_dfr(
    cluster_levels,
    function(cluster_level) {
      coefficient_matrix <- as.matrix(coefficient_list[[cluster_level]])
      tibble::tibble(
        horizon_weeks = horizon_weeks,
        class = cluster_level,
        term = rownames(coefficient_matrix),
        estimate = as.numeric(coefficient_matrix[, 1]),
        lambda_min = as.numeric(cv_fit$lambda.min),
        lambda_1se = as.numeric(cv_fit$lambda.1se)
      )
    }
  ) %>%
    mutate(
      feature_in_model = .data$term %in% c("(Intercept)", feature_names)
    )
}

seizures <- readr::read_csv(SEIZURES_INPUT_PATH, show_col_types = FALSE)
patient_month_panel <- readr::read_csv(PANEL_INPUT_PATH, show_col_types = FALSE)
cluster_assignments <- readr::read_csv(CLUSTER_ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE)

check_columns(
  seizures,
  c("patient_id", "type", "date", "duration", "during_sleep"),
  "Seizure input"
)
check_columns(
  patient_month_panel,
  c("patient_id", "variant_p", "study_start_date", "study_end_date"),
  "Patient-month panel"
)
check_columns(
  cluster_assignments,
  c("patient_id", "pam_k3"),
  "Cluster assignments"
)

patient_windows <- patient_month_panel %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    variant_p = as.character(.data$variant_p),
    study_start_date = as.Date(.data$study_start_date),
    study_end_date = as.Date(.data$study_end_date)
  ) %>%
  group_by(.data$patient_id) %>%
  summarise(
    variant_p = first_non_missing(.data$variant_p),
    study_start_date = min(.data$study_start_date, na.rm = TRUE),
    study_end_date = max(.data$study_end_date, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(
    !is.na(.data$patient_id),
    !is.na(.data$study_start_date),
    !is.na(.data$study_end_date),
    .data$study_start_date <= .data$study_end_date
  )

cluster_lookup <- cluster_assignments %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    pam_k3 = factor(as.character(.data$pam_k3))
  ) %>%
  filter(!is.na(.data$patient_id), !is.na(.data$pam_k3)) %>%
  distinct(.data$patient_id, .keep_all = TRUE)

cluster_levels <- levels(droplevels(cluster_lookup$pam_k3))

seizure_events <- seizures %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    seizure_type = as.character(.data$type),
    event_datetime = parse_event_datetime(.data$date),
    event_date = as.Date(.data$event_datetime),
    duration_seconds = suppressWarnings(as.numeric(.data$duration)),
    during_sleep = stringr::str_to_lower(stringr::str_squish(as.character(.data$during_sleep)))
  ) %>%
  filter(!is.na(.data$patient_id), !is.na(.data$event_date))

feature_table <- purrr::map_dfr(
  HORIZON_WEEKS,
  function(horizon_value) {
    horizon_week_count <- as.integer(horizon_value)
    horizon_days <- as.integer(horizon_week_count * 7L)

    eligible_windows <- patient_windows %>%
      mutate(
        horizon_weeks = horizon_week_count,
        horizon_days = horizon_days,
        window_start_date = .data$study_start_date,
        window_end_date = .data$study_start_date + horizon_days - 1L,
        observed_days = as.integer(.data$study_end_date - .data$study_start_date + 1L)
      ) %>%
      filter(.data$observed_days >= horizon_days)

    week_panel <- eligible_windows %>%
      select("patient_id", "horizon_weeks", "window_start_date") %>%
      tidyr::expand_grid(week_index = seq_len(horizon_week_count)) %>%
      mutate(
        week_start_date = .data$window_start_date + (.data$week_index - 1L) * 7L,
        week_end_date = .data$week_start_date + 6L
      )

    early_events <- seizure_events %>%
      inner_join(
        eligible_windows %>% select("patient_id", "window_start_date", "window_end_date"),
        by = "patient_id"
      ) %>%
      filter(
        .data$event_date >= .data$window_start_date,
        .data$event_date <= .data$window_end_date
      ) %>%
      mutate(
        horizon_weeks = horizon_week_count,
        days_since_window_start = as.integer(.data$event_date - .data$window_start_date),
        week_index = as.integer(floor(.data$days_since_window_start / 7L) + 1L),
        during_sleep_known = .data$during_sleep %in% c("yes", "no", "true", "false"),
        during_sleep_yes = .data$during_sleep %in% c("yes", "true")
      )

    weekly_counts <- week_panel %>%
      left_join(
        early_events %>%
          count(.data$patient_id, .data$horizon_weeks, .data$week_index, name = "weekly_seizure_count"),
        by = c("patient_id", "horizon_weeks", "week_index")
      ) %>%
      mutate(weekly_seizure_count = tidyr::replace_na(.data$weekly_seizure_count, 0L))

    weekly_features <- weekly_counts %>%
      group_by(.data$patient_id, .data$horizon_weeks) %>%
      summarise(
        total_seizures = sum(.data$weekly_seizure_count),
        mean_weekly_seizures = mean(.data$weekly_seizure_count),
        median_weekly_seizures = stats::median(.data$weekly_seizure_count),
        max_weekly_seizures = max(.data$weekly_seizure_count),
        sd_weekly_seizures = safe_sd(.data$weekly_seizure_count),
        iqr_weekly_seizures = stats::IQR(.data$weekly_seizure_count),
        proportion_zero_seizure_weeks = mean(.data$weekly_seizure_count == 0),
        .groups = "drop"
      )

    event_features <- early_events %>%
      arrange(.data$patient_id, .data$horizon_weeks, .data$event_datetime) %>%
      group_by(.data$patient_id, .data$horizon_weeks) %>%
      mutate(
        interseizure_interval_days = as.numeric(
          difftime(.data$event_datetime, lag(.data$event_datetime), units = "days")
        )
      ) %>%
      summarise(
        first_seizure_day = min(.data$days_since_window_start) + 1L,
        last_seizure_day = max(.data$days_since_window_start) + 1L,
        mean_interseizure_interval_days = safe_mean(.data$interseizure_interval_days, default = horizon_days),
        n_seizure_types = dplyr::n_distinct(.data$seizure_type[!is.na(.data$seizure_type)]),
        mean_duration_seconds = safe_mean(.data$duration_seconds),
        max_duration_seconds = safe_max(.data$duration_seconds),
        proportion_during_sleep_yes = if (sum(.data$during_sleep_known) == 0) {
          0
        } else {
          mean(.data$during_sleep_yes[.data$during_sleep_known])
        },
        during_sleep_known_proportion = mean(.data$during_sleep_known),
        .groups = "drop"
      )

    eligible_windows %>%
      select(
        "patient_id",
        "variant_p",
        "horizon_weeks",
        "horizon_days",
        "window_start_date",
        "window_end_date",
        "observed_days"
      ) %>%
      left_join(weekly_features, by = c("patient_id", "horizon_weeks")) %>%
      left_join(event_features, by = c("patient_id", "horizon_weeks")) %>%
      mutate(
        across(
          all_of(c(
            "first_seizure_day",
            "last_seizure_day",
            "mean_interseizure_interval_days",
            "n_seizure_types",
            "mean_duration_seconds",
            "max_duration_seconds",
            "proportion_during_sleep_yes",
            "during_sleep_known_proportion"
          )),
          ~ tidyr::replace_na(.x, 0)
        ),
        first_seizure_day = if_else(.data$total_seizures == 0, .data$horizon_days, .data$first_seizure_day),
        mean_interseizure_interval_days = if_else(
          .data$total_seizures < 2,
          .data$horizon_days,
          .data$mean_interseizure_interval_days
        ),
        seizure_rate_per_28_days = .data$total_seizures / .data$horizon_days * 28,
        first_seizure_day_fraction = .data$first_seizure_day / .data$horizon_days,
        last_seizure_day_fraction = .data$last_seizure_day / .data$horizon_days,
        log1p_total_seizures = log1p(.data$total_seizures),
        log1p_seizure_rate_per_28_days = log1p(.data$seizure_rate_per_28_days),
        log1p_mean_weekly_seizures = log1p(.data$mean_weekly_seizures),
        log1p_median_weekly_seizures = log1p(.data$median_weekly_seizures),
        log1p_max_weekly_seizures = log1p(.data$max_weekly_seizures),
        log1p_sd_weekly_seizures = log1p(.data$sd_weekly_seizures),
        log1p_iqr_weekly_seizures = log1p(.data$iqr_weekly_seizures),
        log1p_mean_interseizure_interval_days = log1p(.data$mean_interseizure_interval_days),
        log1p_mean_duration_seconds = log1p(.data$mean_duration_seconds),
        log1p_max_duration_seconds = log1p(.data$max_duration_seconds)
      )
  }
) %>%
  inner_join(cluster_lookup, by = "patient_id") %>%
  arrange(.data$horizon_weeks, .data$patient_id)

if (nrow(feature_table) == 0) {
  stop("No eligible patients had enough observed follow-up for the requested horizons.", call. = FALSE)
}

if (any(!is.finite(as.matrix(feature_table %>% select(all_of(MODEL_FEATURES)))))) {
  stop("Early seizure feature matrix contains non-finite values.", call. = FALSE)
}

readr::write_csv(feature_table, FEATURES_OUTPUT_PATH)

set.seed(RANDOM_SEED)

prediction_table <- purrr::map_dfr(
  HORIZON_WEEKS,
  function(horizon_weeks) {
    horizon_data <- feature_table %>%
      filter(.data$horizon_weeks == !!horizon_weeks) %>%
      mutate(pam_k3 = factor(as.character(.data$pam_k3), levels = cluster_levels)) %>%
      arrange(.data$patient_id)

    if (n_distinct(horizon_data$pam_k3) < length(cluster_levels)) {
      stop("Horizon ", horizon_weeks, " weeks does not contain all cluster classes.", call. = FALSE)
    }

    purrr::map_dfr(
      seq_len(nrow(horizon_data)),
      function(test_index) {
        training_data <- horizon_data[-test_index, , drop = FALSE]
        test_data <- horizon_data[test_index, , drop = FALSE]

        x_train <- training_data %>% select(all_of(MODEL_FEATURES)) %>% as.matrix()
        y_train <- factor(as.character(training_data$pam_k3), levels = cluster_levels)
        x_test <- test_data %>% select(all_of(MODEL_FEATURES)) %>% as.matrix()

        cv_fit <- fit_cv_ridge_multinomial(x_train, y_train)
        probability_row <- extract_probability_row(
          predict(cv_fit, newx = x_test, s = "lambda.min", type = "response"),
          cluster_levels
        )
        predicted_class <- names(which.max(probability_row))

        tibble::tibble(
          horizon_weeks = horizon_weeks,
          patient_id = test_data$patient_id,
          variant_p = test_data$variant_p,
          true_pam_k3 = as.character(test_data$pam_k3),
          predicted_pam_k3 = predicted_class,
          correct = predicted_class == as.character(test_data$pam_k3),
          lambda_min = as.numeric(cv_fit$lambda.min),
          lambda_1se = as.numeric(cv_fit$lambda.1se)
        ) %>%
          bind_cols(
            as_tibble_row(setNames(as.numeric(probability_row), paste0("prob_cluster_", cluster_levels)))
          )
      }
    )
  }
)

readr::write_csv(prediction_table, PREDICTIONS_OUTPUT_PATH)

performance_table <- prediction_table %>%
  group_by(.data$horizon_weeks) %>%
  group_modify(~ build_metrics(.x, cluster_levels)) %>%
  ungroup() %>%
  arrange(.data$horizon_weeks)

readr::write_csv(performance_table, PERFORMANCE_OUTPUT_PATH)

confusion_table <- prediction_table %>%
  count(.data$horizon_weeks, true_pam_k3 = .data$true_pam_k3, predicted_pam_k3 = .data$predicted_pam_k3, name = "n") %>%
  complete(
    horizon_weeks = HORIZON_WEEKS,
    true_pam_k3 = cluster_levels,
    predicted_pam_k3 = cluster_levels,
    fill = list(n = 0L)
  ) %>%
  arrange(.data$horizon_weeks, .data$true_pam_k3, .data$predicted_pam_k3)

readr::write_csv(confusion_table, CONFUSION_OUTPUT_PATH)

auc_table <- prediction_table %>%
  group_by(.data$horizon_weeks) %>%
  group_modify(~ build_auc_table(.x, cluster_levels)) %>%
  ungroup() %>%
  group_by(.data$horizon_weeks) %>%
  mutate(macro_one_vs_rest_auc = mean(.data$one_vs_rest_auc, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(.data$horizon_weeks, .data$class)

readr::write_csv(auc_table, AUC_OUTPUT_PATH)

coefficient_table <- purrr::map_dfr(
  HORIZON_WEEKS,
  function(horizon_weeks) {
    horizon_data <- feature_table %>%
      filter(.data$horizon_weeks == !!horizon_weeks) %>%
      mutate(pam_k3 = factor(as.character(.data$pam_k3), levels = cluster_levels))

    cv_fit <- fit_cv_ridge_multinomial(
      horizon_data %>% select(all_of(MODEL_FEATURES)) %>% as.matrix(),
      horizon_data$pam_k3
    )

    build_coefficient_table(cv_fit, horizon_weeks, MODEL_FEATURES, cluster_levels)
  }
)

readr::write_csv(coefficient_table, COEFFICIENTS_OUTPUT_PATH)
