suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(glmmTMB)
  library(broom.mixed)
  library(performance)
})

source("src/analysis_config.R")

INPUT_PATH <- "output/tabs/modeling/patient_month_panel.csv"
CLUSTER_ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
OUTPUT_DIR <- "output/tabs/modeling"
COEFFICIENTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "mixed_effects_model_coefficients.csv")
FIT_SUMMARY_OUTPUT_PATH <- file.path(OUTPUT_DIR, "mixed_effects_model_fit_summary.csv")
DIAGNOSTICS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "mixed_effects_model_diagnostics.csv")
RESIDUALS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "mixed_effects_model_residuals.csv")
ZERO_CALIBRATION_BY_PATIENT_OUTPUT_PATH <- file.path(OUTPUT_DIR, "mixed_effects_model_zero_calibration_by_patient.csv")
ZERO_CALIBRATION_BY_CLUSTER_OUTPUT_PATH <- file.path(OUTPUT_DIR, "mixed_effects_model_zero_calibration_by_cluster.csv")
ZI_COEFFICIENTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "zero_inflated_mixed_effects_model_coefficients.csv")
ZI_FIT_SUMMARY_OUTPUT_PATH <- file.path(OUTPUT_DIR, "zero_inflated_mixed_effects_model_fit_summary.csv")
ZI_DIAGNOSTICS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "zero_inflated_mixed_effects_model_diagnostics.csv")
ZI_RESIDUALS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "zero_inflated_mixed_effects_model_residuals.csv")
ZI_ZERO_CALIBRATION_BY_PATIENT_OUTPUT_PATH <- file.path(OUTPUT_DIR, "zero_inflated_mixed_effects_model_zero_calibration_by_patient.csv")
ZI_ZERO_CALIBRATION_BY_CLUSTER_OUTPUT_PATH <- file.path(OUTPUT_DIR, "zero_inflated_mixed_effects_model_zero_calibration_by_cluster.csv")
MODEL_COMPARISON_OUTPUT_PATH <- file.path(OUTPUT_DIR, "mixed_effects_model_comparison.csv")

MODEL_NAME <- "patient_month_cluster_adjusted_mixed_effects_nb_main_effects"
ZI_MODEL_NAME <- "patient_month_cluster_adjusted_covariate_informed_zero_inflated_mixed_effects_nb_main_effects"
BASE_TERMS <- c(
  "month_index",
  "active_med_count",
  "med_started_flag",
  "med_stopped_flag"
)
MEDICATION_TERMS <- c(
  "med_exposed_carbamazepine",
  "med_exposed_clobazam",
  "med_exposed_oxcarbazepine",
  "med_exposed_valproic_acid",
  "med_exposed_cannabidiol",
  "med_exposed_lacosamide",
  "med_exposed_zonisamide",
  "med_exposed_lamotrigine"
)
CLUSTER_TERMS <- c("pam_k3")
FIXED_EFFECT_TERMS <- c(BASE_TERMS, MEDICATION_TERMS, CLUSTER_TERMS)
NUMERIC_FIXED_EFFECT_TERMS <- c(BASE_TERMS, MEDICATION_TERMS)
REQUIRED_COLUMNS <- c(
  "patient_id",
  "month",
  "study_start_date",
  "study_end_date",
  "seizure_count",
  FIXED_EFFECT_TERMS
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_PATH)) {
  stop("Input file not found: ", INPUT_PATH)
}

if (!file.exists(CLUSTER_ASSIGNMENTS_INPUT_PATH)) {
  stop("Cluster assignments file not found: ", CLUSTER_ASSIGNMENTS_INPUT_PATH)
}

patient_month_panel <- readr::read_csv(INPUT_PATH, show_col_types = FALSE)
cluster_assignments <- readr::read_csv(CLUSTER_ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE)

cluster_required_columns <- c("patient_id", "pam_k3")
missing_cluster_columns <- setdiff(cluster_required_columns, names(cluster_assignments))
if (length(missing_cluster_columns) > 0) {
  stop("Cluster assignments are missing required columns: ", paste(missing_cluster_columns, collapse = ", "))
}

cluster_lookup <- cluster_assignments %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    pam_k3 = factor(as.character(.data$pam_k3), levels = c("1", "2", "3"))
  ) %>%
  dplyr::distinct(.data$patient_id, .keep_all = TRUE)

patient_month_panel <- patient_month_panel %>%
  dplyr::mutate(patient_id = as.character(.data$patient_id)) %>%
  dplyr::left_join(cluster_lookup, by = "patient_id")

patient_zero_rates <- patient_month_panel %>%
  dplyr::mutate(seizure_count = suppressWarnings(as.numeric(.data$seizure_count))) %>%
  dplyr::filter(!is.na(.data$seizure_count)) %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(
    patient_zero_rate = mean(.data$seizure_count == 0),
    .groups = "drop"
  )

patient_month_panel <- patient_month_panel %>%
  dplyr::left_join(patient_zero_rates, by = "patient_id")

missing_columns <- setdiff(REQUIRED_COLUMNS, names(patient_month_panel))
if (length(missing_columns) > 0) {
  stop("Input data is missing required columns: ", paste(missing_columns, collapse = ", "))
}

modeling_data <- patient_month_panel %>%
  dplyr::mutate(
    patient_id = factor(.data$patient_id),
    month = as.Date(.data$month),
    study_start_date = as.Date(.data$study_start_date),
    study_end_date = pmin(as.Date(.data$study_end_date), ANALYSIS_CUTOFF_DATE),
    patient_zero_rate = suppressWarnings(as.numeric(.data$patient_zero_rate)),
    next_month_start = as.Date(format(.data$month + 32, "%Y-%m-01")),
    month_end = .data$next_month_start - 1,
    observed_start = as.Date(pmax(.data$study_start_date, .data$month), origin = "1970-01-01"),
    observed_end = as.Date(pmin(.data$study_end_date, .data$month_end), origin = "1970-01-01"),
    observed_days_in_month = pmax(as.integer(.data$observed_end - .data$observed_start + 1), 0L),
    pam_k3 = factor(as.character(.data$pam_k3), levels = c("1", "2", "3")),
    across(
      all_of(c("seizure_count", NUMERIC_FIXED_EFFECT_TERMS)),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) %>%
  dplyr::select(all_of(c(REQUIRED_COLUMNS, "patient_zero_rate", "observed_days_in_month"))) %>%
  dplyr::filter(
    .data$month <= ANALYSIS_CUTOFF_DATE,
    if_all(all_of(REQUIRED_COLUMNS), ~ !is.na(.x)),
    !is.na(.data$patient_zero_rate),
    !is.na(.data$observed_days_in_month),
    .data$observed_days_in_month > 0
  ) %>%
  dplyr::mutate(
    zi_month_index_scaled = as.numeric(scale(.data$month_index)),
    zi_active_med_count_scaled = as.numeric(scale(.data$active_med_count)),
    zi_patient_zero_rate_scaled = as.numeric(scale(.data$patient_zero_rate))
  )

if (nrow(modeling_data) == 0) {
  stop("No complete cases available for model fitting after filtering required columns.")
}

count_formula <- as.formula(
  paste(
    "seizure_count ~",
    paste(FIXED_EFFECT_TERMS, collapse = " + "),
    "+ offset(log(observed_days_in_month))",
    "+ (1 | patient_id)"
  )
)

zi_formula <- ~ zi_month_index_scaled + zi_active_med_count_scaled + pam_k3 + zi_patient_zero_rate_scaled

build_coefficient_table <- function(model, include_zero_inflation = FALSE) {
  cond_table <- broom.mixed::tidy(
    model,
    effects = "fixed",
    component = "cond"
  ) %>%
    dplyr::mutate(
      component = "conditional",
      exp_estimate = exp(.data$estimate),
      exp_conf.low = exp(.data$estimate - z_critical * .data$std.error),
      exp_conf.high = exp(.data$estimate + z_critical * .data$std.error),
      effect_scale = "incidence_rate_ratio",
      conf.low = .data$estimate - z_critical * .data$std.error,
      conf.high = .data$estimate + z_critical * .data$std.error
    )

  if (!include_zero_inflation) {
    return(
      cond_table %>%
        dplyr::select(
          all_of(c(
            "component",
            "term",
            "estimate",
            "std.error",
            "statistic",
            "p.value",
            "conf.low",
            "conf.high",
            "effect_scale",
            "exp_estimate",
            "exp_conf.low",
            "exp_conf.high"
          ))
        )
    )
  }

  zi_table <- broom.mixed::tidy(
    model,
    effects = "fixed",
    component = "zi"
  ) %>%
    dplyr::mutate(
      component = "zero_inflation",
      exp_estimate = exp(.data$estimate),
      exp_conf.low = exp(.data$estimate - z_critical * .data$std.error),
      exp_conf.high = exp(.data$estimate + z_critical * .data$std.error),
      effect_scale = "zero_inflation_odds_ratio",
      conf.low = .data$estimate - z_critical * .data$std.error,
      conf.high = .data$estimate + z_critical * .data$std.error
    )

  dplyr::bind_rows(cond_table, zi_table) %>%
    dplyr::select(
      all_of(c(
        "component",
        "term",
        "estimate",
        "std.error",
        "statistic",
        "p.value",
        "conf.low",
        "conf.high",
        "effect_scale",
        "exp_estimate",
        "exp_conf.low",
        "exp_conf.high"
      ))
    )
}

build_fit_summary <- function(model, model_name, model_engine, family_label, modeling_data, zi_formula_label = "none") {
  tibble::tibble(
    model_name = model_name,
    model_engine = model_engine,
    family = family_label,
    zero_inflation_formula = zi_formula_label,
    nobs = stats::nobs(model),
    number_of_patients = dplyr::n_distinct(modeling_data$patient_id),
    total_person_days = sum(modeling_data$observed_days_in_month),
    min_observed_days_in_month = min(modeling_data$observed_days_in_month),
    median_observed_days_in_month = stats::median(modeling_data$observed_days_in_month),
    max_observed_days_in_month = max(modeling_data$observed_days_in_month),
    partial_month_rows = sum(modeling_data$observed_days_in_month < 28),
    AIC = stats::AIC(model),
    BIC = stats::BIC(model),
    logLik = as.numeric(stats::logLik(model)),
    deviance = tryCatch(
      as.numeric(stats::deviance(model)),
      error = function(...) NA_real_
    )
  )
}

build_diagnostics_table <- function(model) {
  fitted_count <- predict(model, type = "response")
  pearson_residual <- residuals(model, type = "pearson")
  diagnostic_overdispersion <- suppressWarnings(performance::check_overdispersion(model))
  diagnostic_zero_inflation <- suppressWarnings(performance::check_zeroinflation(model))
  diagnostic_convergence <- performance::check_convergence(model)
  diagnostic_singularity <- performance::check_singularity(model)
  diagnostic_r2 <- suppressWarnings(performance::r2_nakagawa(model))
  diagnostic_icc <- suppressWarnings(performance::icc(model))
  diagnostic_performance <- suppressWarnings(performance::model_performance(model))

  tibble::tibble(
    metric = c(
      "converged",
      "singular_fit",
      "overdispersion_dispersion_ratio",
      "overdispersion_chisq",
      "overdispersion_p_value",
      "observed_zeros",
      "predicted_zeros",
      "zero_ratio_observed_to_predicted",
      "zero_inflation_flag",
      "r2_conditional",
      "r2_marginal",
      "icc_adjusted",
      "icc_conditional",
      "rmse",
      "sigma",
      "random_intercept_sd",
      "dispersion_parameter",
      "min_fitted_count",
      "mean_fitted_count",
      "median_fitted_count",
      "max_abs_pearson_residual"
    ),
    value = c(
      as.character(isTRUE(diagnostic_convergence)),
      as.character(isTRUE(diagnostic_singularity)),
      as.character(unclass(diagnostic_overdispersion)$dispersion_ratio),
      as.character(unclass(diagnostic_overdispersion)$chisq_statistic),
      as.character(unclass(diagnostic_overdispersion)$p_value),
      as.character(unclass(diagnostic_zero_inflation)$observed.zeros),
      as.character(unclass(diagnostic_zero_inflation)$predicted.zeros),
      as.character(unclass(diagnostic_zero_inflation)$ratio),
      as.character(unclass(diagnostic_zero_inflation)$ratio < (1 - unclass(diagnostic_zero_inflation)$tolerance)),
      as.character(unclass(diagnostic_r2)$R2_conditional),
      as.character(unclass(diagnostic_r2)$R2_marginal),
      as.character(unclass(diagnostic_icc)$ICC_adjusted),
      as.character(unclass(diagnostic_icc)$ICC_conditional),
      as.character(unclass(diagnostic_performance)$RMSE),
      as.character(unclass(diagnostic_performance)$Sigma),
      as.character(as.numeric(attr(VarCorr(model)$cond$patient_id, "stddev"))),
      as.character(sigma(model)),
      as.character(min(fitted_count)),
      as.character(mean(fitted_count)),
      as.character(stats::median(fitted_count)),
      as.character(max(abs(pearson_residual)))
    )
  )
}

build_residuals_table <- function(model, modeling_data) {
  fitted_count <- predict(model, type = "response")
  predicted_zero_probability <- predict_zero_probability(model)
  pearson_residual <- residuals(model, type = "pearson")
  model_family <- stats::family(model)$family
  deviance_residual <- if (grepl("^truncated_", model_family)) {
    rep(NA_real_, nrow(modeling_data))
  } else {
    residuals(model, type = "deviance")
  }

  modeling_data %>%
    dplyr::transmute(
      patient_id = as.character(.data$patient_id),
      month = .data$month,
      seizure_count = .data$seizure_count,
      observed_days_in_month = .data$observed_days_in_month,
      observed_seizures_per_day = .data$seizure_count / .data$observed_days_in_month,
      fitted_count = fitted_count,
      fitted_seizures_per_day = fitted_count / .data$observed_days_in_month,
      predicted_zero_probability = predicted_zero_probability,
      pearson_residual = pearson_residual,
      deviance_residual = deviance_residual
    ) %>%
    dplyr::arrange(dplyr::desc(abs(.data$pearson_residual)))
}

predict_zero_probability <- function(model) {
  conditional_mean <- predict(model, type = "conditional")
  zero_inflation_probability <- predict(model, type = "zprob")
  dispersion_parameter <- sigma(model)
  conditional_zero_probability <- (dispersion_parameter / (dispersion_parameter + conditional_mean))^dispersion_parameter

  zero_inflation_probability + (1 - zero_inflation_probability) * conditional_zero_probability
}

build_zero_calibration_by_patient <- function(model, modeling_data) {
  predicted_zero_probability <- predict_zero_probability(model)

  modeling_data %>%
    dplyr::transmute(
      patient_id = as.character(.data$patient_id),
      pam_k3 = as.character(.data$pam_k3),
      observed_zero = as.integer(.data$seizure_count == 0),
      predicted_zero_probability = predicted_zero_probability
    ) %>%
    dplyr::group_by(.data$patient_id, .data$pam_k3) %>%
    dplyr::summarise(
      patient_months = dplyr::n(),
      observed_zero_months = sum(.data$observed_zero),
      expected_zero_months = sum(.data$predicted_zero_probability),
      observed_zero_rate = mean(.data$observed_zero),
      expected_zero_rate = mean(.data$predicted_zero_probability),
      zero_month_difference = .data$observed_zero_months - .data$expected_zero_months,
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(abs(.data$zero_month_difference)))
}

build_zero_calibration_by_cluster <- function(model, modeling_data) {
  predicted_zero_probability <- predict_zero_probability(model)

  modeling_data %>%
    dplyr::transmute(
      pam_k3 = as.character(.data$pam_k3),
      observed_zero = as.integer(.data$seizure_count == 0),
      predicted_zero_probability = predicted_zero_probability
    ) %>%
    dplyr::group_by(.data$pam_k3) %>%
    dplyr::summarise(
      patient_months = dplyr::n(),
      observed_zero_months = sum(.data$observed_zero),
      expected_zero_months = sum(.data$predicted_zero_probability),
      observed_zero_rate = mean(.data$observed_zero),
      expected_zero_rate = mean(.data$predicted_zero_probability),
      zero_month_difference = .data$observed_zero_months - .data$expected_zero_months,
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$pam_k3)
}

z_critical <- stats::qnorm(0.975)

mixed_effects_model <- glmmTMB::glmmTMB(
  formula = count_formula,
  data = modeling_data,
  family = glmmTMB::nbinom2(link = "log")
)

zi_mixed_effects_model <- glmmTMB::glmmTMB(
  formula = count_formula,
  ziformula = zi_formula,
  data = modeling_data,
  family = glmmTMB::nbinom2(link = "log")
)

coefficient_table <- build_coefficient_table(mixed_effects_model, include_zero_inflation = FALSE)
zi_coefficient_table <- build_coefficient_table(zi_mixed_effects_model, include_zero_inflation = TRUE)

model_fit_summary <- build_fit_summary(
  mixed_effects_model,
  model_name = MODEL_NAME,
  model_engine = "glmmTMB",
  family_label = "nbinom2",
  modeling_data = modeling_data,
  zi_formula_label = "none"
)

zi_model_fit_summary <- build_fit_summary(
  zi_mixed_effects_model,
  model_name = ZI_MODEL_NAME,
  model_engine = "glmmTMB",
  family_label = "nbinom2",
  modeling_data = modeling_data,
  zi_formula_label = "~ scale(month_index) + scale(active_med_count) + pam_k3 + scale(patient_zero_rate)"
)

diagnostics_table <- build_diagnostics_table(mixed_effects_model)
zi_diagnostics_table <- build_diagnostics_table(zi_mixed_effects_model)

residuals_table <- build_residuals_table(mixed_effects_model, modeling_data)
zi_residuals_table <- build_residuals_table(zi_mixed_effects_model, modeling_data)
zero_calibration_by_patient <- build_zero_calibration_by_patient(mixed_effects_model, modeling_data)
zero_calibration_by_cluster <- build_zero_calibration_by_cluster(mixed_effects_model, modeling_data)
zi_zero_calibration_by_patient <- build_zero_calibration_by_patient(zi_mixed_effects_model, modeling_data)
zi_zero_calibration_by_cluster <- build_zero_calibration_by_cluster(zi_mixed_effects_model, modeling_data)

model_comparison_table <- dplyr::bind_rows(model_fit_summary, zi_model_fit_summary) %>%
  dplyr::left_join(
    dplyr::bind_rows(
      diagnostics_table %>%
        dplyr::filter(.data$metric %in% c("observed_zeros", "predicted_zeros", "zero_ratio_observed_to_predicted", "zero_inflation_flag")) %>%
        tidyr::pivot_wider(names_from = "metric", values_from = "value") %>%
        dplyr::mutate(model_name = MODEL_NAME),
      zi_diagnostics_table %>%
        dplyr::filter(.data$metric %in% c("observed_zeros", "predicted_zeros", "zero_ratio_observed_to_predicted", "zero_inflation_flag")) %>%
        tidyr::pivot_wider(names_from = "metric", values_from = "value") %>%
        dplyr::mutate(model_name = ZI_MODEL_NAME)
    ),
    by = "model_name"
  )

readr::write_csv(coefficient_table, COEFFICIENTS_OUTPUT_PATH)
readr::write_csv(model_fit_summary, FIT_SUMMARY_OUTPUT_PATH)
readr::write_csv(diagnostics_table, DIAGNOSTICS_OUTPUT_PATH)
readr::write_csv(residuals_table, RESIDUALS_OUTPUT_PATH)
readr::write_csv(zero_calibration_by_patient, ZERO_CALIBRATION_BY_PATIENT_OUTPUT_PATH)
readr::write_csv(zero_calibration_by_cluster, ZERO_CALIBRATION_BY_CLUSTER_OUTPUT_PATH)
readr::write_csv(zi_coefficient_table, ZI_COEFFICIENTS_OUTPUT_PATH)
readr::write_csv(zi_model_fit_summary, ZI_FIT_SUMMARY_OUTPUT_PATH)
readr::write_csv(zi_diagnostics_table, ZI_DIAGNOSTICS_OUTPUT_PATH)
readr::write_csv(zi_residuals_table, ZI_RESIDUALS_OUTPUT_PATH)
readr::write_csv(zi_zero_calibration_by_patient, ZI_ZERO_CALIBRATION_BY_PATIENT_OUTPUT_PATH)
readr::write_csv(zi_zero_calibration_by_cluster, ZI_ZERO_CALIBRATION_BY_CLUSTER_OUTPUT_PATH)
readr::write_csv(model_comparison_table, MODEL_COMPARISON_OUTPUT_PATH)
