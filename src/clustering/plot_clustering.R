suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(tidyr)
})

FEATURES_INPUT_PATH <- "output/tabs/clustering/seizure_freq_features.csv"
ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
OUTPUT_FIG_DIR <- "output/figs/clustering"

RATE_FEATURES <- c(
  "mean_monthly_seizure_rate",
  "iqr_monthly_seizure_rate"
)
DIRECT_SCALE_FEATURES <- c("proportion_zero_seizure_months")
PCA_COMPONENT_COUNT <- 2L

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
  left_join(assignments, by = c("patient_id", "variant_p"))

p_pca <- ggplot2::ggplot(
  pca_plot_data,
  ggplot2::aes(x = .data$pc1, y = .data$pc2, color = .data$pam_k3)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
  ggplot2::geom_point(size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = .data$variant_p),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.25,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_brewer(palette = "Set2", na.translate = FALSE) +
  ggplot2::labs(
    title = "Patient Seizure-Frequency Clusters",
    x = "PC1",
    y = "PC2",
    color = "Cluster"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "top"
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

feature_plot_specs <- tibble::tribble(
  ~feature_name, ~plot_label, ~y_label, ~filename,
  "mean_monthly_seizure_rate", "Mean Monthly Seizure Rate", "Mean seizures per month", "boxplot_mean_seizure_rate.png",
  "iqr_monthly_seizure_rate", "IQR of Monthly Seizure Rate", "IQR of seizures per month", "boxplot_iqr.png",
  "proportion_zero_seizure_months", "Proportion of Zero-Seizure Months", "Proportion of months with zero seizures", "boxplot_zero_seizure_months.png"
)

for (i in seq_len(nrow(feature_plot_specs))) {
  spec <- feature_plot_specs[i, ]
  feature_name <- as.character(spec$feature_name[[1]])
  feature_sym <- rlang::sym(feature_name)

  p_box <- ggstatsplot::ggbetweenstats(
    data = boxplot_data,
    x = pam_k3,
    y = !!feature_sym,
    type = "np",
    pairwise.comparisons = TRUE,
    pairwise.display = "significant",
    p.adjust.method = "holm",
    package = "RColorBrewer",
    palette = "Set2",
    title = as.character(spec$plot_label[[1]]),
    xlab = "Cluster",
    ylab = as.character(spec$y_label[[1]]),
    ggtheme = ggplot2::theme_classic(base_size = 11),
    messages = FALSE
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )

  ggplot2::ggsave(
    filename = file.path(OUTPUT_FIG_DIR, as.character(spec$filename[[1]])),
    plot = p_box,
    width = 8,
    height = 6.5,
    dpi = 300
  )
}
