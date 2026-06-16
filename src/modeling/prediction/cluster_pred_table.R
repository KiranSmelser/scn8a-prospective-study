# Format patient-level cluster predictions as confusion tables.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

PREDICTION_DIRS <- c(
  "output/tabs/modeling/prediction/seizure",
  "output/tabs/modeling/prediction/registry",
  "output/tabs/modeling/prediction/combined"
)

format_confusion_table <- function(predictions) {
  required_columns <- c("horizon_weeks", "true_pam_k3", "predicted_pam_k3")
  missing_columns <- setdiff(required_columns, names(predictions))

  if (length(missing_columns) > 0) {
    stop(
      "Input file is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  cluster_levels <- sort(unique(c(
    as.character(predictions$true_pam_k3),
    as.character(predictions$predicted_pam_k3)
  )))

  confusion_long <- predictions %>%
    transmute(
      horizon_weeks = as.integer(.data$horizon_weeks),
      true_cluster = as.character(.data$true_pam_k3),
      predicted_cluster = as.character(.data$predicted_pam_k3)
    ) %>%
    count(.data$horizon_weeks, .data$true_cluster, .data$predicted_cluster, name = "n") %>%
    complete(
      horizon_weeks = sort(unique(predictions$horizon_weeks)),
      true_cluster = cluster_levels,
      predicted_cluster = cluster_levels,
      fill = list(n = 0L)
    ) %>%
    group_by(.data$horizon_weeks, .data$true_cluster) %>%
    mutate(
      row_total = sum(.data$n),
      row_percent = if_else(.data$row_total > 0, .data$n / .data$row_total, NA_real_),
      cell = if_else(
        is.na(.data$row_percent),
        "0 (NA)",
        sprintf("%d (%.1f%%)", .data$n, 100 * .data$row_percent)
      )
    ) %>%
    ungroup()

  confusion_wide <- confusion_long %>%
    transmute(
      horizon_weeks = .data$horizon_weeks,
      true_cluster = .data$true_cluster,
      predicted_cluster = paste0("predicted_cluster_", .data$predicted_cluster),
      cell = .data$cell
    ) %>%
    pivot_wider(
      names_from = "predicted_cluster",
      values_from = "cell"
    )

  cluster_accuracy <- confusion_long %>%
    filter(.data$true_cluster == .data$predicted_cluster) %>%
    transmute(
      horizon_weeks = .data$horizon_weeks,
      true_cluster = .data$true_cluster,
      row_total = .data$row_total,
      cluster_accuracy = .data$row_percent,
      cluster_accuracy_label = if_else(
        is.na(.data$row_percent),
        NA_character_,
        sprintf("%.1f%%", 100 * .data$row_percent)
      )
    )

  confusion_wide %>%
    left_join(cluster_accuracy, by = c("horizon_weeks", "true_cluster")) %>%
    arrange(.data$horizon_weeks, .data$true_cluster)
}

for (output_dir in PREDICTION_DIRS) {
  input_path <- file.path(output_dir, "rolling_ridge_multinomial_patient_predictions.csv")
  output_path <- file.path(output_dir, "rolling_patient_confusion_matrix_formatted.csv")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(input_path)) {
    stop("Input file not found: ", input_path, call. = FALSE)
  }

  predictions <- readr::read_csv(input_path, show_col_types = FALSE)
  formatted_confusion_table <- format_confusion_table(predictions)
  readr::write_csv(formatted_confusion_table, output_path)
}
