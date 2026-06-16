# Predict cluster membership from seizure data and targeted registry data.

suppressPackageStartupMessages({
  library(dplyr)
  library(glmnet)
  library(pROC)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

REGISTRY_INPUT_PATH <- "data/registry.csv"
SEIZURE_FEATURES_INPUT_PATH <- "output/tabs/modeling/prediction/seizure/rolling_seizure_window_features.csv"
CLUSTER_ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
OUTPUT_DIR <- "output/tabs/modeling/prediction/combined"

ROLLING_FEATURES_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_combined_window_features.csv")
ROLLING_WINDOW_PREDICTIONS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_window_predictions.csv")
ROLLING_WINDOW_PERFORMANCE_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_window_performance.csv")
ROLLING_WINDOW_CONFUSION_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_window_confusion_matrix.csv")
ROLLING_WINDOW_AUC_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_window_auc.csv")
ROLLING_PATIENT_PREDICTIONS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_patient_predictions.csv")
ROLLING_PATIENT_PERFORMANCE_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_patient_performance.csv")
ROLLING_PATIENT_CONFUSION_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_patient_confusion_matrix.csv")
ROLLING_PERIOD_PERFORMANCE_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_period_performance.csv")
ROLLING_COEFFICIENTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_coefficients.csv")
FEATURE_MAP_OUTPUT_PATH <- file.path(OUTPUT_DIR, "combined_model_feature_map.csv")

HORIZON_WEEKS <- c(4L, 6L, 8L)
RANDOM_SEED <- 20260616L
RIDGE_ALPHA <- 0
INNER_CV_MAX_FOLDS <- 5L

SEIZURE_MODEL_FEATURES <- c(
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

TARGETED_REGISTRY_SOURCE_COLUMNS <- c(
  "age_seizure_onset_months",
  "initial_tonic",
  "initial_tonic_clonic_grand_mal",
  "seizure_type_tonic",
  "seizure_type_tonic_clonic_grand_mal",
  "seizure_type_focal_aware_simple_partial_seizure",
  "seizure_type_focal_impaired_awareness_complex_partial_seizure_limbic_psychomotor",
  "dev_skill_brush_teeth_with_no_help",
  "dev_skill_name_colors",
  "dev_skill_wash_and_dry_hands",
  "dev_skill_used_a_2_word_combination",
  "dev_skill_spoken_in_phrases",
  "dev_skill_read"
)

REGISTRY_MODEL_FEATURES <- c(
  "registry_tonic_clonic_history",
  "registry_tonic_history",
  "registry_focal_aware_or_impaired",
  "registry_log1p_age_seizure_onset_months",
  "registry_age_seizure_onset_missing",
  "registry_higher_dev_language_adl_score"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(REGISTRY_INPUT_PATH, SEIZURE_FEATURES_INPUT_PATH, CLUSTER_ASSIGNMENTS_INPUT_PATH)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Required input file(s) not found: ", paste(missing_files, collapse = ", "), call. = FALSE)
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
    n_observations = nrow(data),
    n_patients = dplyr::n_distinct(data$patient_id),
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

add_predicted_class <- function(data, cluster_levels) {
  probability_columns <- paste0("prob_cluster_", cluster_levels)
  probability_matrix <- data %>%
    select(all_of(probability_columns)) %>%
    as.matrix()

  data %>%
    mutate(
      predicted_pam_k3 = cluster_levels[max.col(probability_matrix, ties.method = "first")],
      correct = .data$predicted_pam_k3 == as.character(.data$true_pam_k3)
    )
}

build_confusion_table <- function(data, horizon_values, cluster_levels) {
  data %>%
    count(.data$horizon_weeks, true_pam_k3 = .data$true_pam_k3, predicted_pam_k3 = .data$predicted_pam_k3, name = "n") %>%
    complete(
      horizon_weeks = horizon_values,
      true_pam_k3 = cluster_levels,
      predicted_pam_k3 = cluster_levels,
      fill = list(n = 0L)
    ) %>%
    arrange(.data$horizon_weeks, .data$true_pam_k3, .data$predicted_pam_k3)
}

build_performance_table <- function(data, cluster_levels) {
  data %>%
    group_by(.data$horizon_weeks) %>%
    group_modify(~ build_metrics(.x, cluster_levels)) %>%
    ungroup() %>%
    arrange(.data$horizon_weeks)
}

build_multiclass_auc_table <- function(data, cluster_levels) {
  data %>%
    group_by(.data$horizon_weeks) %>%
    group_modify(~ build_auc_table(.x, cluster_levels)) %>%
    ungroup() %>%
    group_by(.data$horizon_weeks) %>%
    mutate(macro_one_vs_rest_auc = mean(.data$one_vs_rest_auc, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(.data$horizon_weeks, .data$class)
}

row_max_binary <- function(...) {
  values <- cbind(...)
  values <- apply(values, 2, function(x) suppressWarnings(as.numeric(x)))
  row_value <- apply(
    values,
    1,
    function(x) {
      if (all(is.na(x))) {
        return(0)
      }
      max(x, na.rm = TRUE)
    }
  )
  as.numeric(row_value > 0)
}

row_mean_score <- function(...) {
  values <- cbind(...)
  values <- apply(values, 2, function(x) suppressWarnings(as.numeric(x)))
  row_value <- rowMeans(values, na.rm = TRUE)
  row_value[is.nan(row_value)] <- 0
  row_value
}

prepare_targeted_registry_features <- function(registry) {
  onset_months <- suppressWarnings(as.numeric(registry$age_seizure_onset_months))
  onset_fill_value <- stats::median(onset_months[is.finite(onset_months)], na.rm = TRUE)
  if (!is.finite(onset_fill_value)) {
    onset_fill_value <- 0
  }
  onset_months_imputed <- onset_months
  onset_months_imputed[!is.finite(onset_months_imputed)] <- onset_fill_value

  registry %>%
    transmute(
      patient_id = as.character(.data$patient_id),
      registry_tonic_clonic_history = suppressWarnings(as.numeric(.data$seizure_type_tonic_clonic_grand_mal)),
      registry_tonic_history = row_max_binary(
        .data$initial_tonic,
        .data$seizure_type_tonic,
        .data$initial_tonic_clonic_grand_mal
      ),
      registry_focal_aware_or_impaired = row_max_binary(
        .data$seizure_type_focal_aware_simple_partial_seizure,
        .data$seizure_type_focal_impaired_awareness_complex_partial_seizure_limbic_psychomotor
      ),
      registry_log1p_age_seizure_onset_months = log1p(onset_months_imputed),
      registry_age_seizure_onset_missing = as.numeric(!is.finite(onset_months)),
      registry_higher_dev_language_adl_score = row_mean_score(
        .data$dev_skill_brush_teeth_with_no_help,
        .data$dev_skill_name_colors,
        .data$dev_skill_wash_and_dry_hands,
        .data$dev_skill_used_a_2_word_combination,
        .data$dev_skill_spoken_in_phrases,
        .data$dev_skill_read
      )
    ) %>%
    filter(!is.na(.data$patient_id)) %>%
    mutate(
      across(all_of(REGISTRY_MODEL_FEATURES), ~ tidyr::replace_na(suppressWarnings(as.numeric(.x)), 0))
    ) %>%
    distinct(.data$patient_id, .keep_all = TRUE)
}

registry_raw <- readr::read_csv(REGISTRY_INPUT_PATH, show_col_types = FALSE)
seizure_feature_table <- readr::read_csv(SEIZURE_FEATURES_INPUT_PATH, show_col_types = FALSE)
cluster_assignments <- readr::read_csv(CLUSTER_ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE)

check_columns(registry_raw, c("patient_id", TARGETED_REGISTRY_SOURCE_COLUMNS), "Registry input")
check_columns(
  seizure_feature_table,
  c(
    "window_id", "patient_id", "variant_p", "horizon_weeks", "period_index",
    "window_start_date", "window_end_date", "pam_k3", SEIZURE_MODEL_FEATURES
  ),
  "Rolling seizure features"
)
check_columns(cluster_assignments, c("patient_id", "pam_k3"), "Cluster assignments")

cluster_lookup <- cluster_assignments %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    pam_k3 = factor(as.character(.data$pam_k3))
  ) %>%
  filter(!is.na(.data$patient_id), !is.na(.data$pam_k3)) %>%
  distinct(.data$patient_id, .keep_all = TRUE)

cluster_levels <- levels(droplevels(cluster_lookup$pam_k3))

registry_design_table <- prepare_targeted_registry_features(registry_raw)
COMBINED_MODEL_FEATURES <- c(SEIZURE_MODEL_FEATURES, REGISTRY_MODEL_FEATURES)

feature_map <- bind_rows(
  tibble::tibble(
    model_feature = SEIZURE_MODEL_FEATURES,
    feature_source = "seizure_window",
    source_column = SEIZURE_MODEL_FEATURES
  ),
  tibble::tibble(
    model_feature = REGISTRY_MODEL_FEATURES,
    feature_source = "targeted_registry",
    source_column = c(
      "seizure_type_tonic_clonic_grand_mal",
      "initial_tonic|seizure_type_tonic|initial_tonic_clonic_grand_mal",
      "seizure_type_focal_aware_simple_partial_seizure|seizure_type_focal_impaired_awareness_complex_partial_seizure_limbic_psychomotor",
      "age_seizure_onset_months",
      "age_seizure_onset_months",
      paste(
        c(
          "dev_skill_brush_teeth_with_no_help",
          "dev_skill_name_colors",
          "dev_skill_wash_and_dry_hands",
          "dev_skill_used_a_2_word_combination",
          "dev_skill_spoken_in_phrases",
          "dev_skill_read"
        ),
        collapse = "|"
      )
    )
  )
)
readr::write_csv(feature_map, FEATURE_MAP_OUTPUT_PATH)

rolling_feature_table <- seizure_feature_table %>%
  mutate(
    patient_id = as.character(.data$patient_id),
    pam_k3 = factor(as.character(.data$pam_k3), levels = cluster_levels)
  ) %>%
  inner_join(registry_design_table, by = "patient_id") %>%
  arrange(.data$horizon_weeks, .data$patient_id, .data$period_index)

if (nrow(rolling_feature_table) == 0) {
  stop("No combined rolling windows were available for the requested horizons.", call. = FALSE)
}

if (any(!is.finite(as.matrix(rolling_feature_table %>% select(all_of(COMBINED_MODEL_FEATURES)))))) {
  stop("Combined feature matrix contains non-finite values.", call. = FALSE)
}

readr::write_csv(rolling_feature_table, ROLLING_FEATURES_OUTPUT_PATH)

set.seed(RANDOM_SEED)

rolling_prediction_table <- purrr::map_dfr(
  HORIZON_WEEKS,
  function(horizon_weeks) {
    horizon_data <- rolling_feature_table %>%
      filter(.data$horizon_weeks == !!horizon_weeks) %>%
      mutate(pam_k3 = factor(as.character(.data$pam_k3), levels = cluster_levels)) %>%
      arrange(.data$patient_id, .data$period_index)

    if (n_distinct(horizon_data$pam_k3) < length(cluster_levels)) {
      stop("Rolling horizon ", horizon_weeks, " weeks does not contain all cluster classes.", call. = FALSE)
    }

    purrr::map_dfr(
      sort(unique(horizon_data$patient_id)),
      function(test_patient_id) {
        training_data <- horizon_data %>%
          filter(.data$patient_id != test_patient_id)
        test_data <- horizon_data %>%
          filter(.data$patient_id == test_patient_id)

        if (n_distinct(training_data$pam_k3) < length(cluster_levels)) {
          stop(
            "Training data for held-out patient ",
            test_patient_id,
            " at horizon ",
            horizon_weeks,
            " weeks does not contain all cluster classes.",
            call. = FALSE
          )
        }

        x_train <- training_data %>% select(all_of(COMBINED_MODEL_FEATURES)) %>% as.matrix()
        y_train <- factor(as.character(training_data$pam_k3), levels = cluster_levels)
        x_test <- test_data %>% select(all_of(COMBINED_MODEL_FEATURES)) %>% as.matrix()

        cv_fit <- fit_cv_ridge_multinomial(x_train, y_train)
        probability_array <- predict(cv_fit, newx = x_test, s = "lambda.min", type = "response")
        predicted_classes <- dimnames(probability_array)[[2]]
        probability_matrix <- matrix(
          0,
          nrow = nrow(test_data),
          ncol = length(cluster_levels),
          dimnames = list(NULL, cluster_levels)
        )
        probability_matrix[, predicted_classes] <- probability_array[, predicted_classes, 1, drop = FALSE]

        test_data %>%
          transmute(
            horizon_weeks = .data$horizon_weeks,
            window_id = .data$window_id,
            patient_id = .data$patient_id,
            variant_p = .data$variant_p,
            period_index = .data$period_index,
            window_start_date = .data$window_start_date,
            window_end_date = .data$window_end_date,
            true_pam_k3 = as.character(.data$pam_k3),
            lambda_min = as.numeric(cv_fit$lambda.min),
            lambda_1se = as.numeric(cv_fit$lambda.1se)
          ) %>%
          bind_cols(
            as_tibble(probability_matrix, .name_repair = ~ paste0("prob_cluster_", .x))
          ) %>%
          add_predicted_class(cluster_levels)
      }
    )
  }
)

readr::write_csv(rolling_prediction_table, ROLLING_WINDOW_PREDICTIONS_OUTPUT_PATH)

rolling_window_performance_table <- build_performance_table(rolling_prediction_table, cluster_levels)
readr::write_csv(rolling_window_performance_table, ROLLING_WINDOW_PERFORMANCE_OUTPUT_PATH)

rolling_window_confusion_table <- build_confusion_table(rolling_prediction_table, HORIZON_WEEKS, cluster_levels)
readr::write_csv(rolling_window_confusion_table, ROLLING_WINDOW_CONFUSION_OUTPUT_PATH)

rolling_window_auc_table <- build_multiclass_auc_table(rolling_prediction_table, cluster_levels)
readr::write_csv(rolling_window_auc_table, ROLLING_WINDOW_AUC_OUTPUT_PATH)

probability_columns <- paste0("prob_cluster_", cluster_levels)

rolling_patient_prediction_table <- rolling_prediction_table %>%
  group_by(.data$horizon_weeks, .data$patient_id, .data$variant_p, .data$true_pam_k3) %>%
  summarise(
    n_windows = n(),
    first_window_start_date = min(.data$window_start_date),
    last_window_end_date = max(.data$window_end_date),
    across(all_of(probability_columns), mean),
    .groups = "drop"
  ) %>%
  add_predicted_class(cluster_levels) %>%
  arrange(.data$horizon_weeks, .data$patient_id)

readr::write_csv(rolling_patient_prediction_table, ROLLING_PATIENT_PREDICTIONS_OUTPUT_PATH)

rolling_patient_performance_table <- build_performance_table(rolling_patient_prediction_table, cluster_levels)
readr::write_csv(rolling_patient_performance_table, ROLLING_PATIENT_PERFORMANCE_OUTPUT_PATH)

rolling_patient_confusion_table <- build_confusion_table(rolling_patient_prediction_table, HORIZON_WEEKS, cluster_levels)
readr::write_csv(rolling_patient_confusion_table, ROLLING_PATIENT_CONFUSION_OUTPUT_PATH)

rolling_period_performance_table <- rolling_prediction_table %>%
  group_by(.data$horizon_weeks, .data$period_index) %>%
  filter(dplyr::n_distinct(.data$true_pam_k3) == length(cluster_levels)) %>%
  group_modify(~ build_metrics(.x, cluster_levels)) %>%
  ungroup() %>%
  arrange(.data$horizon_weeks, .data$period_index)

readr::write_csv(rolling_period_performance_table, ROLLING_PERIOD_PERFORMANCE_OUTPUT_PATH)

rolling_coefficient_table <- purrr::map_dfr(
  HORIZON_WEEKS,
  function(horizon_weeks) {
    horizon_data <- rolling_feature_table %>%
      filter(.data$horizon_weeks == !!horizon_weeks) %>%
      mutate(pam_k3 = factor(as.character(.data$pam_k3), levels = cluster_levels))

    cv_fit <- fit_cv_ridge_multinomial(
      horizon_data %>% select(all_of(COMBINED_MODEL_FEATURES)) %>% as.matrix(),
      horizon_data$pam_k3
    )

    build_coefficient_table(cv_fit, horizon_weeks, COMBINED_MODEL_FEATURES, cluster_levels)
  }
) %>%
  left_join(
    feature_map,
    by = c("term" = "model_feature")
  ) %>%
  mutate(
    feature_source = if_else(.data$term == "(Intercept)", "intercept", .data$feature_source),
    source_column = if_else(.data$term == "(Intercept)", "(Intercept)", .data$source_column)
  )

readr::write_csv(rolling_coefficient_table, ROLLING_COEFFICIENTS_OUTPUT_PATH)
