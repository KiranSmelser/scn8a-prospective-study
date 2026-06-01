# Forest plot for fixed-effect coefficients from the zero-inflated model.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

COEFFICIENTS_INPUT_PATH <- "output/tabs/modeling/zero_inflated_mixed_effects_model_coefficients.csv"
OUTPUT_FIG_DIR <- "output/figs/modeling"
OUTPUT_TAB_DIR <- "output/tabs/modeling"
PLOT_DATA_OUTPUT_PATH <- file.path(OUTPUT_TAB_DIR, "zero_inflated_mixed_effects_model_coefficient_plot_data.csv")
PDF_OUTPUT_PATH <- file.path(OUTPUT_FIG_DIR, "zero_inflated_mixed_effects_model_coefficients.pdf")

dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(COEFFICIENTS_INPUT_PATH)) {
  stop("Coefficient input file not found: ", COEFFICIENTS_INPUT_PATH, call. = FALSE)
}

required_columns <- c(
  "component",
  "term",
  "exp_estimate",
  "exp_conf.low",
  "exp_conf.high",
  "p.value"
)

coefficient_data <- readr::read_csv(COEFFICIENTS_INPUT_PATH, show_col_types = FALSE)
missing_columns <- setdiff(required_columns, names(coefficient_data))

if (length(missing_columns) > 0) {
  stop(
    "Coefficient input file is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

term_labels <- c(
  month_index = "Study month",
  active_med_count = "Active medication count",
  med_started_flag = "Medication started this month",
  med_stopped_flag = "Medication stopped this month",
  med_exposed_carbamazepine = "Carbamazepine exposure",
  med_exposed_clobazam = "Clobazam exposure",
  med_exposed_oxcarbazepine = "Oxcarbazepine exposure",
  med_exposed_valproic_acid = "Valproic acid exposure",
  med_exposed_cannabidiol = "Cannabidiol exposure",
  med_exposed_lacosamide = "Lacosamide exposure",
  med_exposed_zonisamide = "Zonisamide exposure",
  med_exposed_lamotrigine = "Lamotrigine exposure",
  pam_k32 = "Cluster 2 vs 1",
  pam_k33 = "Cluster 3 vs 1",
  zi_month_index_scaled = "Study month (per 1 SD increase)",
  zi_active_med_count_scaled = "Active medication count (per 1 SD increase)",
  zi_patient_zero_rate_scaled = "Patient zero-seizure month rate (per 1 SD increase)"
)

plot_data <- coefficient_data %>%
  dplyr::filter(.data$term != "(Intercept)") %>%
  dplyr::mutate(
    component = as.character(.data$component),
    term = as.character(.data$term),
    term_label = dplyr::recode(.data$term, !!!term_labels, .default = .data$term),
    exp_estimate = suppressWarnings(as.numeric(.data$exp_estimate)),
    exp_conf.low = suppressWarnings(as.numeric(.data$exp_conf.low)),
    exp_conf.high = suppressWarnings(as.numeric(.data$exp_conf.high)),
    p.value = suppressWarnings(as.numeric(.data$p.value)),
    statistically_significant = !is.na(.data$p.value) & .data$p.value < 0.05
  )

if (nrow(plot_data) == 0) {
  stop("No non-intercept coefficients are available to plot.", call. = FALSE)
}

if (any(!is.finite(plot_data$exp_estimate)) ||
    any(!is.finite(plot_data$exp_conf.low)) ||
    any(!is.finite(plot_data$exp_conf.high)) ||
    any(plot_data$exp_conf.low <= 0) ||
    any(plot_data$exp_conf.high <= 0)) {
  stop("Exponentiated coefficient estimates and confidence limits must be finite and greater than zero.", call. = FALSE)
}

readr::write_csv(
  plot_data %>%
    dplyr::select(
      "component",
      "term",
      "term_label",
      "exp_estimate",
      "exp_conf.low",
      "exp_conf.high",
      "p.value",
      "statistically_significant"
    ),
  PLOT_DATA_OUTPUT_PATH
)

make_component_plot <- function(data, title, x_label) {
  term_order <- data %>%
    dplyr::pull(.data$term_label) %>%
    rev()

  data <- data %>%
    dplyr::mutate(term_label = factor(.data$term_label, levels = term_order))

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$exp_estimate,
      y = .data$term_label,
      xmin = .data$exp_conf.low,
      xmax = .data$exp_conf.high
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linewidth = 0.45,
      color = "grey45",
      linetype = "dashed"
    ) +
    ggplot2::geom_errorbarh(
      height = 0.16,
      linewidth = 0.55,
      color = "#1F4E79"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$statistically_significant),
      shape = 21,
      size = 2.8,
      stroke = 0.65,
      color = "#1F4E79"
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "#1F4E79", `FALSE` = "white"),
      breaks = c(TRUE, FALSE),
      labels = c("p < 0.05", "p >= 0.05"),
      name = NULL
    ) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 10.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = "none",
      plot.margin = ggplot2::margin(7, 10, 7, 7)
    )

  p
}

count_data <- plot_data %>%
  dplyr::filter(.data$component == "conditional")

zero_inflation_data <- plot_data %>%
  dplyr::filter(.data$component == "zero_inflation")

if (nrow(count_data) == 0 || nrow(zero_inflation_data) == 0) {
  stop("Both conditional and zero-inflation coefficients are required.", call. = FALSE)
}

p_count <- make_component_plot(
  data = count_data,
  title = "A. Count component",
  x_label = "Incidence rate ratio (95% CI)"
)

p_zero_inflation <- make_component_plot(
  data = zero_inflation_data,
  title = "B. Zero-inflation component",
  x_label = "Odds ratio for an excess zero (95% CI)"
)

draw_figure <- function() {
  grid::grid.newpage()
  figure_layout <- grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit(c(nrow(count_data), nrow(zero_inflation_data) + 1), "null")
  )
  grid::pushViewport(grid::viewport(layout = figure_layout))
  print(p_count, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_zero_inflation, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid::popViewport()
}

grDevices::pdf(PDF_OUTPUT_PATH, width = 8.5, height = 10, useDingbats = FALSE)
draw_figure()
grDevices::dev.off()
