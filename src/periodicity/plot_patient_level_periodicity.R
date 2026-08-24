# Patient-level seizure-count autocorrelation heatmap.

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(readr)
  library(stringr)
})

PERIODICITY_INPUT_PATH <- "output/tabs/periodicity/patient_level_periodicity.csv"
LAG_INPUT_PATH <- "output/tabs/periodicity/patient_level_periodicity_lags.csv"
CLUSTER_ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
OUTPUT_FIG_DIR <- "output/figs/periodicity"
HEATMAP_PDF_OUTPUT_PATH <- file.path(OUTPUT_FIG_DIR, "periodicity_heatmap.pdf")

dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PERIODICITY_INPUT_PATH)) {
  stop("Periodicity input file not found: ", PERIODICITY_INPUT_PATH, call. = FALSE)
}

if (!file.exists(LAG_INPUT_PATH)) {
  stop("Lag diagnostics input file not found: ", LAG_INPUT_PATH, call. = FALSE)
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
      "K1473K + Pro1428_Lys1473del"
  )
  short_patient_id <- stringr::str_sub(patient_id, 1L, 6L)
  stringr::str_wrap(paste0(variant, " (", short_patient_id, ")"), width = 30)
}

format_cluster_label <- function(pam_cluster) {
  dplyr::if_else(
    is.na(pam_cluster),
    "Cluster Unknown",
    paste("Cluster", pam_cluster)
  )
}

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

periodicity_results <- readr::read_csv(PERIODICITY_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    variant_p = as.character(.data$variant_p),
    patient_label = wrap_variant_label(.data$variant_p, .data$patient_id),
    patient_status = dplyr::case_when(
      .data$significant ~ "FDR significant",
      .data$analysis_status == "ok" ~ "Not significant",
      TRUE ~ "Not tested"
    )
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

patient_order <- periodicity_results %>%
  dplyr::arrange(
    .data$pam_cluster_sort,
    dplyr::desc(.data$significant),
    .data$strongest_positive_lag_days,
    dplyr::desc(.data$strongest_positive_autocorrelation),
    .data$patient_label
  ) %>%
  dplyr::pull(.data$patient_label)

heatmap_data <- readr::read_csv(LAG_INPUT_PATH, show_col_types = FALSE) %>%
  dplyr::mutate(
    patient_id = as.character(.data$patient_id),
    variant_p = as.character(.data$variant_p),
    patient_label = wrap_variant_label(.data$variant_p, .data$patient_id)
  ) %>%
  dplyr::left_join(
    periodicity_results %>%
      dplyr::select(
        "patient_id",
        "patient_status",
        "observed_seizure_days",
        "observed_seizure_events",
        "pam_cluster",
        "pam_cluster_label"
      ),
    by = "patient_id"
  ) %>%
  dplyr::filter(.data$analysis_status == "ok") %>%
  dplyr::mutate(
    patient_label = factor(.data$patient_label, levels = rev(patient_order)),
    selected_strongest_positive_lag = as.logical(.data$selected_strongest_positive_lag)
  )

selected_lags <- heatmap_data %>%
  dplyr::filter(.data$selected_strongest_positive_lag)

reference_lags <- c(1, 7, 14, 30, 60, 90)
reference_lags <- reference_lags[
  reference_lags >= min(heatmap_data$lag_days, na.rm = TRUE) &
    reference_lags <= max(heatmap_data$lag_days, na.rm = TRUE)
]

plot_width <- 12
plot_height <- max(6, 0.32 * dplyr::n_distinct(heatmap_data$patient_label) + 2.5)

p_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(x = .data$lag_days, y = .data$patient_label, fill = .data$autocorrelation)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.25) +
  ggplot2::geom_vline(
    xintercept = reference_lags,
    color = "grey35",
    linewidth = 0.25,
    linetype = "dashed",
    alpha = 0.55
  ) +
  ggplot2::geom_point(
    data = selected_lags,
    ggplot2::aes(x = .data$lag_days, y = .data$patient_label),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.1,
    stroke = 0.65,
    color = "black",
    fill = "white"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#C62828",
    midpoint = 0,
    name = "Autocorrelation"
  ) +
  ggplot2::scale_x_continuous(
    breaks = reference_lags,
    expand = c(0, 0)
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(pam_cluster_label),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::labs(
    title = "Patient-Level Seizure Periodicity",
    x = "Days",
    y = "Patient"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(hjust = 1, vjust = 1),
    axis.text.y = ggplot2::element_text(size = 7),
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    legend.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 10, color = "grey25"),
    strip.background = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.text.y.left = ggplot2::element_text(face = "bold", angle = 0, size = 8),
    legend.position = "right"
  )

ggplot2::ggsave(
  filename = HEATMAP_PDF_OUTPUT_PATH,
  plot = p_heatmap,
  width = plot_width,
  height = plot_height,
  bg = "white"
)
