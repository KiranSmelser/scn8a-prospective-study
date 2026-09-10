# Patient-level seizure burden distribution by seizure type.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(jsonlite)
})

source("src/desc/seizure_type_standardization.R")

OUTPUT_FIG_DIR <- "output/figs/seizure_patterns"
OUTPUT_TAB_DIR <- "output/tabs/seizure_patterns"
dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
png_output_path <- file.path(OUTPUT_FIG_DIR, "patient_level_seizure_burden_by_type_violin_box.png")
if (file.exists(png_output_path)) {
  invisible(file.remove(png_output_path))
}

analysis_end <- analysis_end_date()

events <- read_events_corrected() %>%
  standardize_seizure_events(filter_to_seizure = TRUE) %>%
  dplyr::filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    !is.na(.data$patient_id),
    !is.na(.data$month),
    !is.na(.data$event_date),
    .data$event_date <= analysis_end,
    !.data$requires_review
  ) %>%
  dplyr::mutate(seizure_type_plot = .data$seizure_type_primary) %>%
  dplyr::filter(
    !is.na(.data$seizure_type_plot),
    !.data$seizure_type_plot %in% c(
      "Unknown/Unmapped",
      "Unspecified",
      "Device-detected event",
      "Mixed/Multiple"
    )
  )

usage_months <- if (file.exists("output/tabs/app_usage.json")) {
  usage_records <- jsonlite::fromJSON("output/tabs/app_usage.json", simplifyVector = FALSE)
  purrr::map_dfr(
    usage_records,
    function(rec) {
      patient <- rec$patient_id
      usage <- rec$app_usage_per_month
      if (is.null(patient) || is.null(usage) || length(usage) == 0) {
        return(tibble::tibble())
      }
      usage_vec <- unlist(usage, use.names = TRUE)
      if (is.null(usage_vec) || length(usage_vec) == 0) {
        return(tibble::tibble())
      }
      tibble::tibble(
        patient_id = patient,
        month = as.Date(paste0(names(usage_vec), "-01")),
        app_usage = as.numeric(usage_vec)
      )
    }
  ) %>%
    dplyr::mutate(month = as.Date(lubridate::floor_date(.data$month, "month"))) %>%
    dplyr::filter(
      .data$patient_id %in% TARGET_PATIENT_IDS,
      !is.na(.data$patient_id),
      !is.na(.data$month),
      .data$month <= as.Date(lubridate::floor_date(analysis_end, "month")),
      .data$app_usage >= 1
    ) %>%
    dplyr::distinct(.data$patient_id, .data$month)
} else {
  tibble::tibble()
}

if (nrow(usage_months) == 0) {
  usage_months <- events %>%
    dplyr::distinct(.data$patient_id, .data$month)
}

event_counts <- events %>%
  dplyr::count(.data$patient_id, .data$month, .data$seizure_type_plot, name = "seizure_count")

patient_type_pairs <- events %>%
  dplyr::distinct(.data$patient_id, .data$seizure_type_plot)

monthly_complete <- usage_months %>%
  dplyr::inner_join(patient_type_pairs, by = "patient_id", relationship = "many-to-many") %>%
  dplyr::left_join(event_counts, by = c("patient_id", "month", "seizure_type_plot")) %>%
  tidyr::replace_na(list(seizure_count = 0L))

patient_type_burden <- monthly_complete %>%
  dplyr::group_by(.data$patient_id, .data$seizure_type_plot) %>%
  dplyr::summarise(
    months_recorded = dplyr::n(),
    months_with_type = sum(.data$seizure_count > 0),
    total_type_events = sum(.data$seizure_count),
    avg_events_per_month = mean(.data$seizure_count),
    median_events_per_month = median(.data$seizure_count),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$avg_events_per_month), .data$patient_id, .data$seizure_type_plot)

type_summary <- patient_type_burden %>%
  dplyr::group_by(.data$seizure_type_plot) %>%
  dplyr::summarise(
    n_patients = dplyr::n(),
    median_avg_events_per_month = median(.data$avg_events_per_month),
    mean_avg_events_per_month = mean(.data$avg_events_per_month),
    p25_avg_events_per_month = stats::quantile(.data$avg_events_per_month, probs = 0.25),
    p75_avg_events_per_month = stats::quantile(.data$avg_events_per_month, probs = 0.75),
    max_avg_events_per_month = max(.data$avg_events_per_month),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$median_avg_events_per_month), dplyr::desc(.data$n_patients), .data$seizure_type_plot)

plot_data <- patient_type_burden %>%
  dplyr::inner_join(type_summary %>% dplyr::select(seizure_type_plot), by = "seizure_type_plot") %>%
  dplyr::mutate(
    seizure_type_plot = factor(.data$seizure_type_plot, levels = type_summary$seizure_type_plot)
  )

if (nrow(plot_data) > 0) {
  plot_data_stats <- plot_data %>%
    dplyr::add_count(.data$seizure_type_plot, name = "n_group") %>%
    dplyr::filter(.data$n_group >= 2) %>%
    dplyr::select(-n_group)

  if (dplyr::n_distinct(plot_data_stats$seizure_type_plot) >= 2) {
    p <- ggstatsplot::ggbetweenstats(
      data = plot_data_stats,
      x = seizure_type_plot,
      y = avg_events_per_month,
      type = "np",
      pairwise.comparisons = TRUE,
      pairwise.display = "significant",
      p.adjust.method = "holm",
      package = "RColorBrewer",
      palette = "Set3",
      title = "Seizure Burden by Type",
      xlab = "Seizure type",
      ylab = "Seizures per month",
      ggtheme = ggplot2::theme_classic(base_size = 11),
      messages = FALSE
    ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1),
        plot.title = ggplot2::element_text(face = "bold")
      )

    ggplot2::ggsave(
      filename = file.path(OUTPUT_FIG_DIR, "patient_level_seizure_burden_by_type_violin_box.pdf"),
      plot = p,
      width = max(9, 0.75 * dplyr::n_distinct(plot_data_stats$seizure_type_plot)),
      height = 6.5
    )
  } else {
    warning("Not enough seizure type groups with >=2 observations for pairwise comparisons; plot not generated.")
  }
}

readr::write_csv(
  patient_type_burden,
  file.path(OUTPUT_TAB_DIR, "patient_level_seizure_burden_by_type.csv")
)
readr::write_csv(
  type_summary,
  file.path(OUTPUT_TAB_DIR, "seizure_burden_by_type_summary.csv")
)
readr::write_csv(
  monthly_complete,
  file.path(OUTPUT_TAB_DIR, "seizure_monthly_counts_by_patient_type.csv")
)
