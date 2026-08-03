# Patient-level seizure episodicity figure and summary table.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggh4x)
  library(ggplot2)
  library(readr)
  library(scales)
  library(stringr)
  library(tibble)
})

EPISODICITY_RESULTS_INPUT_PATH <- "output/tabs/episodicity/patient_level_episodicity.csv"
CLUSTER_ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
OUTPUT_FIG_DIR <- "output/figs/episodicity"
OUTPUT_TAB_DIR <- "output/tabs/episodicity"
PLOT_DATA_OUTPUT_PATH <- file.path(OUTPUT_TAB_DIR, "patient_level_episodicity_plot_data.csv")
SUMMARY_TABLE_OUTPUT_PATH <- file.path(OUTPUT_TAB_DIR, "episodicity_summary.csv")
EPISODICITY_PDF_OUTPUT_PATH <- file.path(
  OUTPUT_FIG_DIR,
  "patient_level_episodicity.pdf"
)
PRIMARY_TAU_DAYS <- 1L

dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(EPISODICITY_RESULTS_INPUT_PATH)) {
  stop("Episodicity results input file not found: ", EPISODICITY_RESULTS_INPUT_PATH, call. = FALSE)
}

if (!file.exists(CLUSTER_ASSIGNMENTS_INPUT_PATH)) {
  stop("Cluster assignments input file not found: ", CLUSTER_ASSIGNMENTS_INPUT_PATH, call. = FALSE)
}

wrap_variant_label <- function(variant_p, patient_id) {
  variant <- dplyr::if_else(
    is.na(variant_p) | variant_p == "",
    "Unknown variant",
    variant_p
  )
  short_patient_id <- stringr::str_sub(patient_id, 1L, 6L)
  stringr::str_wrap(paste0(variant, " (", short_patient_id, ")"), width = 28)
}

format_cluster_label <- function(pam_cluster) {
  dplyr::if_else(
    is.na(pam_cluster),
    "Cluster Unknown",
    paste("Cluster", pam_cluster)
  )
}

scope_labels <- c(
  "all seizure types" = "All seizure types",
  "focal-only" = "Focal only",
  "tonic-clonic-only" = "Tonic-clonic only"
)
scope_row_labels <- c(
  "all seizure types" = "All",
  "focal-only" = "Focal",
  "tonic-clonic-only" = "Tonic-clonic"
)

cluster_assignments <- readr::read_csv(CLUSTER_ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE)

missing_cluster_columns <- setdiff(c("patient_id", "pam_k3"), names(cluster_assignments))
if (length(missing_cluster_columns) > 0) {
  stop(
    "Cluster assignments are missing required columns: ",
    paste(missing_cluster_columns, collapse = ", "),
    call. = FALSE
  )
}

cluster_lookup <- cluster_assignments %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    pam_cluster = suppressWarnings(as.integer(.data$pam_k3))
  )

cluster_label_levels <- cluster_lookup %>%
  dplyr::filter(!is.na(.data$pam_cluster)) %>%
  dplyr::distinct(.data$pam_cluster) %>%
  dplyr::arrange(.data$pam_cluster) %>%
  dplyr::pull(.data$pam_cluster) %>%
  format_cluster_label() %>%
  c("Cluster Unknown")

results <- readr::read_csv(EPISODICITY_RESULTS_INPUT_PATH, show_col_types = FALSE)

required_result_columns <- c(
  "seizure_scope",
  "patient_id",
  "variant_p",
  "tau_days",
  "analysis_status",
  "observed_seizure_events",
  "max_clump_events",
  "null_max_clump_events_lower_95",
  "null_max_clump_events_upper_95",
  "q_value_max_clump_events",
  "significant_max_clump_events",
  "fraction_events_in_multi_event_clumps",
  "null_fraction_events_in_multi_event_clumps_lower_95",
  "null_fraction_events_in_multi_event_clumps_upper_95",
  "q_value_fraction_events_in_multi_event_clumps",
  "significant_fraction_events_in_multi_event_clumps",
  "max_clump_duration_days",
  "null_max_clump_duration_days_lower_95",
  "null_max_clump_duration_days_upper_95",
  "q_value_max_clump_duration_days",
  "significant_max_clump_duration_days",
  "significant_any_episodicity_metric"
)

missing_result_columns <- setdiff(required_result_columns, names(results))
if (length(missing_result_columns) > 0) {
  stop(
    "Episodicity results are missing required columns: ",
    paste(missing_result_columns, collapse = ", "),
    call. = FALSE
  )
}

primary_results <- results %>%
  dplyr::filter(.data$tau_days == PRIMARY_TAU_DAYS) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    variant_p = as.character(.data$variant_p),
    seizure_scope_label = factor(
      unname(scope_labels[.data$seizure_scope]),
      levels = unname(scope_labels)
    ),
    patient_label = wrap_variant_label(.data$variant_p, .data$patient_id)
  ) %>%
  dplyr::left_join(cluster_lookup, by = "patient_id") %>%
  dplyr::mutate(
    pam_cluster_sort = dplyr::if_else(
      is.na(.data$pam_cluster),
      Inf,
      as.double(.data$pam_cluster)
    ),
    pam_cluster_label = factor(format_cluster_label(.data$pam_cluster), levels = cluster_label_levels)
  )

if (nrow(primary_results) == 0) {
  stop("No results were found for the primary tau of 1 day.", call. = FALSE)
}

patient_order <- primary_results %>%
  dplyr::filter(.data$seizure_scope == "all seizure types") %>%
  dplyr::mutate(
    max_clump_excess = .data$max_clump_events - .data$null_max_clump_events_lower_95
  ) %>%
  dplyr::arrange(
    .data$pam_cluster_sort,
    dplyr::desc(dplyr::coalesce(.data$max_clump_excess, -Inf)),
    .data$patient_label
  ) %>%
  dplyr::pull(.data$patient_label)

multiple_seizure_type_patients <- primary_results %>%
  dplyr::filter(
    .data$seizure_scope %in% c("focal-only", "tonic-clonic-only"),
    .data$observed_seizure_events > 0
  ) %>%
  dplyr::distinct(.data$patient_id, .data$seizure_scope) %>%
  dplyr::count(.data$patient_id, name = "observed_seizure_types") %>%
  dplyr::filter(.data$observed_seizure_types > 1L) %>%
  dplyr::pull(.data$patient_id)

plot_results <- primary_results %>%
  dplyr::mutate(
    has_multiple_seizure_types = .data$patient_id %in% multiple_seizure_type_patients,
    scope_row_label = unname(scope_row_labels[.data$seizure_scope])
  ) %>%
  dplyr::filter(
    .data$seizure_scope == "all seizure types" |
      (.data$has_multiple_seizure_types & .data$significant_any_episodicity_metric)
  ) %>%
  dplyr::mutate(
    patient_scope_key = paste(.data$patient_id, .data$seizure_scope, sep = "__"),
    patient_sort = match(.data$patient_label, patient_order),
    scope_sort = match(.data$seizure_scope, names(scope_row_labels))
  )

patient_scope_rows <- plot_results %>%
  dplyr::distinct(
    .data$patient_scope_key,
    .data$patient_label,
    .data$scope_row_label,
    .data$has_multiple_seizure_types,
    .data$patient_sort,
    .data$scope_sort
  ) %>%
  dplyr::arrange(.data$patient_sort, .data$scope_sort) %>%
  dplyr::mutate(
    axis_label = dplyr::case_when(
      !.data$has_multiple_seizure_types ~ as.character(.data$patient_label),
      .data$scope_row_label == "All" ~ as.character(.data$patient_label),
      TRUE ~ paste0("  ", .data$scope_row_label)
    )
  )

patient_scope_order <- patient_scope_rows %>%
  dplyr::pull(.data$patient_scope_key)
patient_scope_axis_labels <- stats::setNames(
  patient_scope_rows$axis_label,
  patient_scope_rows$patient_scope_key
)

make_metric_rows <- function(
    data,
    metric,
    metric_label,
    observed,
    null_lower,
    null_upper,
    q_value,
    significant) {
  data %>%
    dplyr::transmute(
      seizure_scope = .data$seizure_scope,
      seizure_scope_label = .data$seizure_scope_label,
      patient_id = .data$patient_id,
      variant_p = .data$variant_p,
      patient_label = .data$patient_label,
      patient_scope_key = .data$patient_scope_key,
      scope_row_label = .data$scope_row_label,
      has_multiple_seizure_types = .data$has_multiple_seizure_types,
      pam_cluster = .data$pam_cluster,
      pam_cluster_sort = .data$pam_cluster_sort,
      pam_cluster_label = .data$pam_cluster_label,
      tau_days = .data$tau_days,
      analysis_status = .data$analysis_status,
      metric = metric,
      metric_label = metric_label,
      observed_value = {{ observed }},
      null_lower_95 = {{ null_lower }},
      null_upper_95 = {{ null_upper }},
      q_value = {{ q_value }},
      significant = {{ significant }}
    )
}

plot_data <- dplyr::bind_rows(
  make_metric_rows(
    plot_results,
    metric = "max_clump_events",
    metric_label = "Maximum events per clump",
    observed = .data$max_clump_events,
    null_lower = .data$null_max_clump_events_lower_95,
    null_upper = .data$null_max_clump_events_upper_95,
    q_value = .data$q_value_max_clump_events,
    significant = .data$significant_max_clump_events
  ),
  make_metric_rows(
    plot_results,
    metric = "fraction_events_in_multi_event_clumps",
    metric_label = "Events in multi-event clumps",
    observed = .data$fraction_events_in_multi_event_clumps,
    null_lower = .data$null_fraction_events_in_multi_event_clumps_lower_95,
    null_upper = .data$null_fraction_events_in_multi_event_clumps_upper_95,
    q_value = .data$q_value_fraction_events_in_multi_event_clumps,
    significant = .data$significant_fraction_events_in_multi_event_clumps
  ),
  make_metric_rows(
    plot_results,
    metric = "max_clump_duration_days",
    metric_label = "Maximum clump duration",
    observed = .data$max_clump_duration_days,
    null_lower = .data$null_max_clump_duration_days_lower_95,
    null_upper = .data$null_max_clump_duration_days_upper_95,
    q_value = .data$q_value_max_clump_duration_days,
    significant = .data$significant_max_clump_duration_days
  )
) %>%
  dplyr::filter(
    .data$analysis_status == "ok",
    !is.na(.data$seizure_scope_label),
    .data$seizure_scope == "all seizure types" | .data$significant
  ) %>%
  dplyr::mutate(
    patient_scope_key = factor(.data$patient_scope_key, levels = rev(patient_scope_order)),
    metric_facet_label = factor(
      .data$metric_label,
      levels = c(
        "Maximum events per clump",
        "Events in multi-event clumps",
        "Maximum clump duration"
      ),
      labels = c(
        "Maximum events\nper clump",
        "Events in multi-event\nclumps (%)",
        "Maximum clump\nduration (days)"
      )
    ),
    effect_category = factor(
      dplyr::if_else(
        .data$significant,
        "FDR significant episodicity",
        "Not significant"
      ),
      levels = c("FDR significant episodicity", "Not significant")
    )
  )

readr::write_csv(plot_data, PLOT_DATA_OUTPUT_PATH)

format_count_percent <- function(count, denominator) {
  ifelse(
    denominator > 0,
    sprintf("%d (%.1f%%)", count, 100 * count / denominator),
    NA_character_
  )
}

summary_table <- results %>%
  dplyr::mutate(
    seizure_scope_label = factor(
      unname(scope_labels[.data$seizure_scope]),
      levels = unname(scope_labels)
    )
  ) %>%
  dplyr::filter(!is.na(.data$seizure_scope_label)) %>%
  dplyr::group_by(.data$seizure_scope, .data$seizure_scope_label, .data$tau_days) %>%
  dplyr::summarise(
    analyzable_n = sum(.data$analysis_status == "ok"),
    max_clump_events_n = sum(.data$significant_max_clump_events & .data$analysis_status == "ok", na.rm = TRUE),
    multi_event_fraction_n = sum(
      .data$significant_fraction_events_in_multi_event_clumps & .data$analysis_status == "ok",
      na.rm = TRUE
    ),
    max_clump_duration_n = sum(
      .data$significant_max_clump_duration_days & .data$analysis_status == "ok",
      na.rm = TRUE
    ),
    any_metric_n = sum(.data$significant_any_episodicity_metric & .data$analysis_status == "ok", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    max_clump_events_q_lt_0_05 = format_count_percent(.data$max_clump_events_n, .data$analyzable_n),
    multi_event_fraction_q_lt_0_05 = format_count_percent(.data$multi_event_fraction_n, .data$analyzable_n),
    max_clump_duration_q_lt_0_05 = format_count_percent(.data$max_clump_duration_n, .data$analyzable_n),
    any_metric_q_lt_0_05 = format_count_percent(.data$any_metric_n, .data$analyzable_n)
  ) %>%
  dplyr::arrange(.data$seizure_scope_label, .data$tau_days) %>%
  dplyr::select(
    seizure_scope = "seizure_scope_label",
    "tau_days",
    "analyzable_n",
    "max_clump_events_q_lt_0_05",
    "multi_event_fraction_q_lt_0_05",
    "max_clump_duration_q_lt_0_05",
    "any_metric_q_lt_0_05"
  )

readr::write_csv(summary_table, SUMMARY_TABLE_OUTPUT_PATH)

effect_colors <- c(
  "FDR significant episodicity" = "#B22222",
  "Not significant" = "#6B7280"
)

p_episodicity <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(y = .data$patient_scope_key)
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$null_lower_95,
      xend = .data$null_upper_95,
      yend = .data$patient_scope_key
    ),
    linewidth = 0.45,
    color = "grey78",
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = .data$observed_value, color = .data$effect_category),
    size = 2.4,
    alpha = 0.95,
    na.rm = TRUE
  ) +
  ggh4x::facet_nested(
    rows = ggplot2::vars(pam_cluster_label),
    cols = ggplot2::vars(metric_facet_label),
    scales = "free",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::scale_y_discrete(labels = patient_scope_axis_labels) +
  ggh4x::facetted_pos_scales(
    x = list(
      metric_facet_label == "Maximum events\nper clump" ~ ggplot2::scale_x_continuous(
        trans = scales::pseudo_log_trans(base = 10),
        breaks = c(1, 3, 10, 30, 100, 300, 1000),
        labels = scales::label_number(accuracy = 1)
      ),
      metric_facet_label == "Events in multi-event\nclumps (%)" ~ ggplot2::scale_x_continuous(
        labels = scales::label_percent(accuracy = 1),
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.25)
      ),
      metric_facet_label == "Maximum clump\nduration (days)" ~ ggplot2::scale_x_continuous(
        trans = scales::pseudo_log_trans(base = 10),
        breaks = c(1, 3, 10, 30, 100, 300),
        labels = scales::label_number(accuracy = 1)
      )
    )
  ) +
  ggplot2::scale_color_manual(values = effect_colors, drop = FALSE) +
  ggplot2::labs(
    title = "Patient-Level Seizure Episodicity",
    x = NULL,
    y = "Patient",
    color = NULL,
    caption = "Grey bars show the central 95% of each patient-scope permutation null distribution. Points show observed values.",
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(size = 10, color = "grey25"),
    plot.caption = ggplot2::element_text(size = 8, color = "grey35", hjust = 0),
    axis.text.y = ggplot2::element_text(size = 7),
    axis.text.x = ggplot2::element_text(size = 8),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(face = "bold", angle = 0, size = 8),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 8),
    panel.spacing.x = grid::unit(0.8, "lines"),
    panel.spacing.y = grid::unit(0.25, "lines")
  )

ggplot2::ggsave(
  filename = EPISODICITY_PDF_OUTPUT_PATH,
  plot = p_episodicity,
  width = 13,
  height = max(
    14,
    0.30 * dplyr::n_distinct(plot_data$patient_scope_key) + 5
  ),
  limitsize = FALSE
)
