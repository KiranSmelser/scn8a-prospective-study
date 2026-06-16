# Predict cluster membership from ONLY registry data.

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
PANEL_INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
CLUSTER_ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
OUTPUT_DIR <- "output/tabs/modeling/prediction/registry"

ROLLING_FEATURES_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_registry_window_features.csv")
ROLLING_WINDOW_PREDICTIONS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_window_predictions.csv")
ROLLING_WINDOW_PERFORMANCE_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_window_performance.csv")
ROLLING_WINDOW_CONFUSION_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_window_confusion_matrix.csv")
ROLLING_WINDOW_AUC_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_window_auc.csv")
ROLLING_PATIENT_PREDICTIONS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_patient_predictions.csv")
ROLLING_PATIENT_PERFORMANCE_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_patient_performance.csv")
ROLLING_PATIENT_CONFUSION_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_patient_confusion_matrix.csv")
ROLLING_PERIOD_PERFORMANCE_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_period_performance.csv")
ROLLING_COEFFICIENTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "rolling_ridge_multinomial_coefficients.csv")
FEATURE_MAP_OUTPUT_PATH <- file.path(OUTPUT_DIR, "registry_model_feature_map.csv")

HORIZON_WEEKS <- c(4L, 6L, 8L)
ROLLING_STEP_DAYS <- 7L
RANDOM_SEED <- 20260616L
RIDGE_ALPHA <- 0
INNER_CV_MAX_FOLDS <- 5L

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(REGISTRY_INPUT_PATH, PANEL_INPUT_PATH, CLUSTER_ASSIGNMENTS_INPUT_PATH)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Required input file(s) not found: ", paste(missing_files, collapse = ", "), call. = FALSE)
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }
  x[[1]]
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

clean_text_vector <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  tidyr::replace_na(x, "Unknown")
}

impute_numeric_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  finite_values <- x[is.finite(x)]
  fill_value <- if (length(finite_values) == 0) {
    0
  } else {
    stats::median(finite_values)
  }
  x[!is.finite(x)] <- fill_value
  x
}

prepare_registry_predictors <- function(registry) {
  registry_predictors <- registry %>%
    mutate(patient_id = as.character(.data$patient_id)) %>%
    filter(!is.na(.data$patient_id)) %>%
    distinct(.data$patient_id, .keep_all = TRUE)

  predictor_columns <- setdiff(names(registry_predictors), "patient_id")
  for (column_name in predictor_columns) {
    column_value <- registry_predictors[[column_name]]

    registry_predictors[[column_name]] <- if (inherits(column_value, "Date") || inherits(column_value, "POSIXt")) {
      impute_numeric_vector(column_value)
    } else if (is.logical(column_value)) {
      impute_numeric_vector(as.integer(column_value))
    } else if (is.numeric(column_value) || is.integer(column_value)) {
      impute_numeric_vector(column_value)
    } else {
      clean_text_vector(column_value)
    }
  }

  registry_predictors
}

registry_raw <- readr::read_csv(REGISTRY_INPUT_PATH, show_col_types = FALSE)
patient_month_panel <- readr::read_csv(PANEL_INPUT_PATH, show_col_types = FALSE)
cluster_assignments <- readr::read_csv(CLUSTER_ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE)

check_columns(registry_raw, c("patient_id"), "Registry input")
check_columns(patient_month_panel, c("patient_id", "variant_p", "study_start_date", "study_end_date"), "Patient-month panel")
check_columns(cluster_assignments, c("patient_id", "pam_k3"), "Cluster assignments")

registry_predictors <- prepare_registry_predictors(registry_raw)

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

registry_model_input <- registry_predictors %>%
  inner_join(cluster_lookup, by = "patient_id") %>%
  arrange(.data$patient_id)

if (nrow(registry_model_input) == 0) {
  stop("No registry rows matched cluster assignments.", call. = FALSE)
}

registry_design_matrix <- stats::model.matrix(pam_k3 ~ . - patient_id, data = registry_model_input) %>%
  as.matrix()
registry_design_matrix <- registry_design_matrix[, colnames(registry_design_matrix) != "(Intercept)", drop = FALSE]

if (ncol(registry_design_matrix) == 0) {
  stop("Registry model matrix contains no predictor columns.", call. = FALSE)
}

if (any(!is.finite(registry_design_matrix))) {
  stop("Registry model matrix contains non-finite values.", call. = FALSE)
}

registry_design_table <- tibble::as_tibble(
  registry_design_matrix,
  .name_repair = "unique"
) %>%
  mutate(patient_id = registry_model_input$patient_id, .before = 1)

MODEL_FEATURES <- setdiff(names(registry_design_table), "patient_id")

feature_map <- tibble::tibble(
  model_feature = MODEL_FEATURES,
  source_column = purrr::map_chr(
    MODEL_FEATURES,
    function(feature_name) {
      source_matches <- names(registry_predictors)[names(registry_predictors) != "patient_id"]
      source_matches <- source_matches[startsWith(feature_name, source_matches)]
      if (length(source_matches) == 0) {
        return(NA_character_)
      }
      source_matches[[which.max(nchar(source_matches))]]
    }
  )
)
readr::write_csv(feature_map, FEATURE_MAP_OUTPUT_PATH)

rolling_windows <- purrr::map_dfr(
  HORIZON_WEEKS,
  function(horizon_value) {
    horizon_week_count <- as.integer(horizon_value)
    horizon_days <- as.integer(horizon_week_count * 7L)

    patient_windows %>%
      mutate(
        horizon_weeks = horizon_week_count,
        horizon_days = horizon_days,
        observed_days = as.integer(.data$study_end_date - .data$study_start_date + 1L),
        last_window_start_date = .data$study_end_date - horizon_days + 1L
      ) %>%
      filter(
        .data$observed_days >= horizon_days,
        .data$last_window_start_date >= .data$study_start_date
      ) %>%
      rowwise() %>%
      mutate(
        window_start_date = list(seq(
          from = .data$study_start_date,
          to = .data$last_window_start_date,
          by = paste(ROLLING_STEP_DAYS, "days")
        ))
      ) %>%
      ungroup() %>%
      select(
        "patient_id",
        "variant_p",
        "horizon_weeks",
        "horizon_days",
        "observed_days",
        "window_start_date"
      ) %>%
      tidyr::unnest("window_start_date") %>%
      group_by(.data$patient_id, .data$horizon_weeks) %>%
      arrange(.data$window_start_date, .by_group = TRUE) %>%
      mutate(
        period_index = row_number(),
        window_end_date = .data$window_start_date + .data$horizon_days - 1L,
        window_id = paste(.data$patient_id, .data$horizon_weeks, .data$period_index, sep = "__")
      ) %>%
      ungroup()
  }
) %>%
  select(
    "window_id",
    "patient_id",
    "variant_p",
    "horizon_weeks",
    "horizon_days",
    "period_index",
    "window_start_date",
    "window_end_date",
    "observed_days"
  ) %>%
  arrange(.data$horizon_weeks, .data$patient_id, .data$period_index)

rolling_feature_table <- rolling_windows %>%
  inner_join(registry_design_table, by = "patient_id") %>%
  inner_join(cluster_lookup, by = "patient_id") %>%
  arrange(.data$horizon_weeks, .data$patient_id, .data$period_index)

if (nrow(rolling_feature_table) == 0) {
  stop("No rolling registry windows were available for the requested horizons.", call. = FALSE)
}

if (any(!is.finite(as.matrix(rolling_feature_table %>% select(all_of(MODEL_FEATURES)))))) {
  stop("Rolling registry feature matrix contains non-finite values.", call. = FALSE)
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

        x_train <- training_data %>% select(all_of(MODEL_FEATURES)) %>% as.matrix()
        y_train <- factor(as.character(training_data$pam_k3), levels = cluster_levels)
        x_test <- test_data %>% select(all_of(MODEL_FEATURES)) %>% as.matrix()

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
      horizon_data %>% select(all_of(MODEL_FEATURES)) %>% as.matrix(),
      horizon_data$pam_k3
    )

    build_coefficient_table(cv_fit, horizon_weeks, MODEL_FEATURES, cluster_levels)
  }
)

readr::write_csv(rolling_coefficient_table, ROLLING_COEFFICIENTS_OUTPUT_PATH)
