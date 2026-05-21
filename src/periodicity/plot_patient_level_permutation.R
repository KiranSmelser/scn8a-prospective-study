# Patient-level post-seizure risk figure.

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(readr)
  library(stringr)
})

PERMUTATION_RESULTS_INPUT_PATH <- "output/tabs/periodicity/patient_level_changepoint_stratified_permutation.csv"
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

wrap_variant_label <- function(variant_p, patient_id) {
  variant <- dplyr::if_else(
    is.na(variant_p) | variant_p == "",
    "Unknown variant",
    variant_p
  )
  short_patient_id <- stringr::str_sub(patient_id, 1L, 6L)
  stringr::str_wrap(paste0(variant, " (", short_patient_id, ")"), width = 28)
}

window_labels <- c(
  lag_1_day = "1 day",
  lag_2_3_days = "2-3 days",
  lag_4_7_days = "4-7 days",
  lag_8_14_days = "8-14 days"
)

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
  dplyr::filter(.data$analysis_status == "ok")

patient_order <- results %>%
  dplyr::filter(.data$window == "1 day") %>%
  dplyr::arrange(
    dplyr::desc(dplyr::coalesce(.data$plot_risk_difference, -Inf)),
    .data$patient_label
  ) %>%
  dplyr::pull(.data$patient_label)

plot_data <- results %>%
  dplyr::mutate(
    patient_label = factor(.data$patient_label, levels = rev(patient_order)),
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
  ggplot2::aes(y = .data$patient_label)
) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55") +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = .data$null_risk_difference_lower_95,
      xend = .data$null_risk_difference_upper_95,
      yend = .data$patient_label
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
    size = 2.4,
    alpha = 0.95,
    na.rm = TRUE
  ) +
  ggplot2::facet_wrap(~window, nrow = 1) +
  ggplot2::scale_color_manual(values = effect_colors) +
  ggplot2::scale_x_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(-0.45, 0.75),
    breaks = seq(-0.4, 0.7, by = 0.2)
  ) +
  ggplot2::labs(
    title = "Patient-Level Post-Seizure Risk",
    x = "Seizure-risk difference",
    y = NULL,
    color = NULL,
    caption = "Grey bars show the central 95% of the patient-window permutation null distribution. Points show observed risk differences."
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(size = 10, color = "grey25"),
    plot.caption = ggplot2::element_text(size = 8, color = "grey35", hjust = 0),
    axis.text.y = ggplot2::element_text(size = 7),
    axis.text.x = ggplot2::element_text(size = 8),
    axis.title.x = ggplot2::element_text(size = 9, margin = ggplot2::margin(t = 8)),
    strip.background = ggplot2::element_rect(fill = "grey92", color = NA),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 8),
    panel.spacing.x = grid::unit(0.8, "lines")
  )

ggplot2::ggsave(
  filename = RISK_DIFFERENCE_PDF_OUTPUT_PATH,
  plot = p_risk_difference,
  width = 13,
  height = 8
)
