# Patient-level post-seizure risk figure.

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(readr)
  library(stringr)
})

PERMUTATION_RESULTS_INPUT_PATH <- "output/tabs/periodicity/patient_level_changepoint_stratified_permutation.csv"
CLUSTER_ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
OUTPUT_FIG_DIR <- "output/figs/periodicity"
OUTPUT_TAB_DIR <- "output/tabs/periodicity"
PLOT_DATA_OUTPUT_PATH <- file.path(OUTPUT_TAB_DIR, "patient_level_permutation_plot_data.csv")
RISK_DIFFERENCE_PDF_OUTPUT_PATH <- file.path(
  OUTPUT_FIG_DIR,
  "patient_level_post_seizure_risk_difference.pdf"
)

dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PERMUTATION_RESULTS_INPUT_PATH)) {
  stop("Permutation results input file not found: ", PERMUTATION_RESULTS_INPUT_PATH, call. = FALSE)
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
  variant <- dplyr::recode(
    variant,
    "K1473K + Pro1428_Lys1473del [predicted inframe exon skipping]" =
      "c.4419+1A>G"
  )
  stringr::str_wrap(variant, width = 28)
}

format_cluster_label <- function(pam_cluster) {
  dplyr::if_else(
    is.na(pam_cluster),
    "Cluster Unknown",
    paste("Cluster", pam_cluster)
  )
}

window_labels <- c(
  lag_1_day = "1 day",
  lag_2_3_days = "2-3 days",
  lag_4_7_days = "4-7 days"
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

results <- readr::read_csv(PERMUTATION_RESULTS_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    variant_p = as.character(.data$variant_p),
    window = factor(.data$window, levels = names(window_labels), labels = unname(window_labels)),
    patient_label = wrap_variant_label(.data$variant_p, .data$patient_id),
    effect_category = dplyr::case_when(
      .data$significant & .data$direction == "increased" ~ "FDR significant increased risk",
      .data$significant & .data$direction == "decreased" ~ "FDR significant decreased risk",
      TRUE ~ "Not significant"
    ),
    plot_risk_difference = .data$risk_difference
  ) %>%
  dplyr::filter(!is.na(.data$window)) %>%
  dplyr::left_join(cluster_lookup, by = "patient_id") %>%
  dplyr::mutate(
    pam_cluster_sort = dplyr::if_else(
      is.na(.data$pam_cluster),
      Inf,
      as.double(.data$pam_cluster)
    ),
    pam_cluster_label = factor(format_cluster_label(.data$pam_cluster), levels = cluster_label_levels)
  ) %>%
  dplyr::filter(.data$analysis_status == "ok")

patient_order <- results %>%
  dplyr::filter(.data$window == "1 day") %>%
  dplyr::arrange(
    .data$pam_cluster_sort,
    dplyr::desc(dplyr::coalesce(.data$plot_risk_difference, -Inf)),
    .data$patient_label
  ) %>%
  dplyr::pull(.data$patient_id) %>%
  unique()

patient_axis_labels <- results %>%
  dplyr::distinct(.data$patient_id, .data$patient_label) %>%
  tibble::deframe()

plot_data <- results %>%
  dplyr::mutate(
    patient_key = factor(.data$patient_id, levels = rev(patient_order)),
    patient_order_rank = match(.data$patient_id, patient_order),
    effect_category = factor(
      .data$effect_category,
      levels = c(
        "FDR significant increased risk",
        "FDR significant decreased risk",
        "Not significant"
      )
    )
  )

readr::write_csv(plot_data, PLOT_DATA_OUTPUT_PATH)

effect_colors <- c(
  "FDR significant increased risk" = "#B22222",
  "FDR significant decreased risk" = "#2B6CB0",
  "Not significant" = "#6B7280"
)

p_risk_difference <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(y = .data$patient_key)
) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55") +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$null_risk_difference_lower_95,
      xend = .data$null_risk_difference_upper_95,
      yend = .data$patient_key
    ),
    linewidth = 0.45,
    color = "grey78",
    na.rm = TRUE
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = .data$plot_risk_difference,
      color = .data$effect_category
    ),
    size = 4,
    alpha = 0.95,
    na.rm = TRUE
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(pam_cluster_label),
    cols = ggplot2::vars(window),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::scale_y_discrete(labels = patient_axis_labels) +
  ggplot2::scale_color_manual(values = effect_colors) +
  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(size = 5))
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(-0.45, 0.75),
    breaks = seq(-0.4, 0.7, by = 0.2)
  ) +
  ggplot2::labs(
    x = "Seizure-risk difference",
    y = NULL,
    color = NULL,
    caption = "Grey bars show the central 95% of the patient-window permutation null distribution. Points show observed risk differences."
  ) +
  ggplot2::theme_classic(base_size = 16) +
  ggplot2::theme(
    plot.caption = ggplot2::element_text(size = 11, color = "grey35", hjust = 0),
    axis.text.y = ggplot2::element_text(size = 12, color = "#222222"),
    axis.text.x = ggplot2::element_text(size = 13, color = "#222222"),
    axis.title.x = ggplot2::element_text(size = 14, margin = ggplot2::margin(t = 8)),
    axis.title.y = ggplot2::element_text(size = 14),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 14),
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(face = "bold", angle = 0, size = 12),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 12),
    legend.key.width = grid::unit(0.6, "cm"),
    panel.spacing.x = grid::unit(0.8, "lines")
  )

ggplot2::ggsave(
  filename = RISK_DIFFERENCE_PDF_OUTPUT_PATH,
  plot = p_risk_difference,
  width = 13,
  height = max(8, 0.32 * dplyr::n_distinct(plot_data$patient_key) + 3)
)
