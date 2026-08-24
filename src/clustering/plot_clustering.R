suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(lubridate)
  library(readr)
  library(tidyr)
})

source("src/analysis_config.R")

FEATURES_INPUT_PATH <- "output/tabs/clustering/seizure_freq_features.csv"
ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
SEIZURES_INPUT_PATH <- "output/tabs/seizures/seizures.csv"
OUTPUT_FIG_DIR <- "output/figs/clustering"

RATE_FEATURES <- c(
  "mean_monthly_seizure_rate",
  "iqr_monthly_seizure_rate"
)
DIRECT_SCALE_FEATURES <- c("proportion_zero_seizure_months")
PCA_COMPONENT_COUNT <- 2L
publication_base_size <- 16
publication_boxplot_base_size <- 18

dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggstatsplot", quietly = TRUE)) {
  stop("Package 'ggstatsplot' is required for boxplots.")
}

if (!file.exists(FEATURES_INPUT_PATH)) {
  stop("Feature input file not found: ", FEATURES_INPUT_PATH)
}

if (!file.exists(ASSIGNMENTS_INPUT_PATH)) {
  stop("Assignment input file not found: ", ASSIGNMENTS_INPUT_PATH)
}

if (!file.exists(SEIZURES_INPUT_PATH)) {
  stop("Seizure input file not found: ", SEIZURES_INPUT_PATH)
}

robust_scale <- function(x) {
  center <- stats::median(x, na.rm = TRUE)
  spread <- stats::IQR(x, na.rm = TRUE)

  if (is.na(spread) || spread == 0) {
    spread <- stats::mad(x, center = center, constant = 1, na.rm = TRUE)
  }

  if (is.na(spread) || spread == 0) {
    spread <- 1
  }

  (x - center) / spread
}

parse_event_datetime <- function(x) {
  parsed <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC"))
  parsed <- dplyr::if_else(
    is.na(parsed),
    suppressWarnings(as.POSIXct(lubridate::ymd(x), tz = "UTC")),
    parsed
  )
  parsed
}

features <- readr::read_csv(FEATURES_INPUT_PATH, show_col_types = FALSE) %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    variant_p = as.character(.data$variant_p),
    across(all_of(RATE_FEATURES), ~ suppressWarnings(as.numeric(.x))),
    across(all_of(DIRECT_SCALE_FEATURES), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  filter(
    !is.na(.data$patient_id),
    if_all(all_of(c(RATE_FEATURES, DIRECT_SCALE_FEATURES)), ~ !is.na(.x))
  ) %>%
  mutate(
    across(all_of(RATE_FEATURES), log1p, .names = "{.col}_log1p"),
    across(all_of(c(paste0(RATE_FEATURES, "_log1p"), DIRECT_SCALE_FEATURES)), robust_scale, .names = "{.col}_scaled")
  )

scaled_feature_names <- c(
  paste0(RATE_FEATURES, "_log1p_scaled"),
  paste0(DIRECT_SCALE_FEATURES, "_scaled")
)

pca_input_matrix <- features %>%
  select(all_of(scaled_feature_names)) %>%
  as.matrix()

if (any(!is.finite(pca_input_matrix))) {
  stop("Scaled feature matrix contains non-finite values.")
}

pca_fit <- stats::prcomp(pca_input_matrix, center = FALSE, scale. = FALSE)

if (ncol(pca_fit$x) < PCA_COMPONENT_COUNT) {
  stop(
    "PCA produced fewer than ",
    PCA_COMPONENT_COUNT,
    " components. Components available: ",
    ncol(pca_fit$x)
  )
}

explained_variance <- summary(pca_fit)$importance[2, seq_len(PCA_COMPONENT_COUNT)]

assignments <- readr::read_csv(ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE) %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    variant_p = as.character(.data$variant_p),
    pam_k3 = factor(
      suppressWarnings(as.integer(.data$pam_k3)),
      levels = c(1, 2, 3),
      labels = c("1", "2", "3")
    )
  )

pca_plot_data <- features %>%
  select("patient_id", "variant_p") %>%
  bind_cols(
    tibble::tibble(
      pc1 = pca_fit$x[, 1],
      pc2 = pca_fit$x[, 2]
    )
  ) %>%
  left_join(assignments, by = c("patient_id", "variant_p")) %>%
  mutate(
    variant_label = dplyr::recode(
      .data$variant_p,
      "K1473K + Pro1428_Lys1473del [predicted inframe exon skipping]" =
        "K1473K + Pro1428_Lys1473del"
    )
  )

p_pca <- ggplot2::ggplot(
  pca_plot_data,
  ggplot2::aes(x = .data$pc1, y = .data$pc2, color = .data$pam_k3)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  ggplot2::geom_point(size = 3.6, alpha = 0.9) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = .data$variant_label),
    size = 4.2,
    box.padding = 0.4,
    point.padding = 0.35,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = EPILEPSIA_CLUSTER_COLORS,
    breaks = names(EPILEPSIA_CLUSTER_COLORS),
    na.translate = FALSE
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(override.aes = list(size = 5))
  ) +
  ggplot2::labs(
    x = "PC1",
    y = "PC2",
    color = "Cluster"
  ) +
  ggplot2::theme_classic(base_size = publication_base_size) +
  ggplot2::theme(
    legend.position = "top",
    legend.title = ggplot2::element_text(size = 14),
    legend.text = ggplot2::element_text(size = 14),
    legend.key.width = grid::unit(0.6, "cm"),
    axis.text = ggplot2::element_text(size = 14, color = "#222222"),
    axis.title = ggplot2::element_text(size = 15)
  )

ggplot2::ggsave(
  filename = file.path(OUTPUT_FIG_DIR, "pca_scatter.png"),
  plot = p_pca,
  width = 8.5,
  height = 6.5,
  dpi = 300
)

boxplot_data <- features %>%
  left_join(assignments %>% select("patient_id", "pam_k3"), by = "patient_id") %>%
  filter(!is.na(.data$pam_k3)) %>%
  mutate(pam_k3 = factor(.data$pam_k3, levels = c("1", "2", "3")))

interseizure_intervals <- readr::read_csv(SEIZURES_INPUT_PATH, show_col_types = FALSE) %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    event_datetime = parse_event_datetime(.data$date)
  ) %>%
  filter(
    !is.na(.data$patient_id),
    !is.na(.data$event_datetime),
    as.Date(.data$event_datetime) <= ANALYSIS_CUTOFF_DATE
  ) %>%
  arrange(.data$patient_id, .data$event_datetime) %>%
  group_by(.data$patient_id) %>%
  mutate(
    days_since_previous_seizure = as.numeric(
      difftime(.data$event_datetime, lag(.data$event_datetime), units = "days")
    )
  ) %>%
  summarise(
    mean_interseizure_interval_days = mean(.data$days_since_previous_seizure, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(.data$mean_interseizure_interval_days))

boxplot_data <- boxplot_data %>%
  left_join(interseizure_intervals, by = "patient_id")

feature_plot_specs <- tibble::tribble(
  ~feature_name, ~plot_label, ~y_label, ~filename,
  "mean_monthly_seizure_rate", "Mean Monthly Seizure Rate", "Mean seizures per month", "boxplot_mean_seizure_rate.png",
  "iqr_monthly_seizure_rate", "IQR of Monthly Seizure Rate", "IQR of seizures per month", "boxplot_iqr.png",
  "proportion_zero_seizure_months", "Proportion of Zero-Seizure Months", "Proportion of months with zero seizures", "boxplot_zero_seizure_months.png",
  "mean_interseizure_interval_days", "Time Between Seizures", "Mean days between seizures", "boxplot_mean_seizure_gaps.png"
)

for (i in seq_len(nrow(feature_plot_specs))) {
  spec <- feature_plot_specs[i, ]
  feature_name <- as.character(spec$feature_name[[1]])
  feature_sym <- rlang::sym(feature_name)
  use_figure2_format <- feature_name != "iqr_monthly_seizure_rate"
  plot_data <- boxplot_data %>%
    filter(!is.na(.data[[feature_name]]))

  if (dplyr::n_distinct(plot_data$pam_k3) < 2) {
    warning("Skipping ", feature_name, ": fewer than 2 clusters with non-missing values.")
    next
  }

  p_box <- ggstatsplot::ggbetweenstats(
    data = plot_data,
    x = pam_k3,
    y = !!feature_sym,
    type = "np",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    p.adjust.method = "holm",
    package = "RColorBrewer",
    palette = "Set2",
    title = if (use_figure2_format) NULL else as.character(spec$plot_label[[1]]),
    xlab = "Cluster",
    ylab = as.character(spec$y_label[[1]]),
    centrality.point.args = if (use_figure2_format) {
      list(size = 6, color = "darkred")
    } else {
      list(size = 5, color = "darkred")
    },
    centrality.label.args = if (use_figure2_format) {
      list(size = 4, nudge_x = 0.4, segment.linetype = 4, min.segment.length = 0)
    } else {
      list(size = 3, nudge_x = 0.4, segment.linetype = 4, min.segment.length = 0)
    },
    point.args = list(
      position = ggplot2::position_jitterdodge(dodge.width = 0.6),
      alpha = 0.4,
      size = if (use_figure2_format) 3.5 else 3,
      stroke = 0,
      na.rm = TRUE
    ),
    ggsignif.args = list(
      textsize = if (use_figure2_format) 4 else 3,
      tip_length = 0.01,
      na.rm = TRUE
    ),
    ggtheme = ggplot2::theme_classic(
      base_size = if (use_figure2_format) publication_boxplot_base_size else 11
    ),
    messages = FALSE
  )

  if (use_figure2_format) {
    p_box <- p_box + ggplot2::theme(
      legend.position = "none",
      plot.subtitle = ggplot2::element_text(size = 13, lineheight = 0.95),
      plot.caption = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 14, color = "#222222"),
      axis.title = ggplot2::element_text(size = 15),
      axis.title.y.right = ggplot2::element_text(size = 13)
    )
  } else {
    p_box <- p_box + ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
  }

  ggplot2::ggsave(
    filename = file.path(OUTPUT_FIG_DIR, as.character(spec$filename[[1]])),
    plot = p_box,
    width = 8,
    height = 6.5,
    dpi = 300
  )
}
