suppressPackageStartupMessages({
  library(dplyr)
  library(ggalluvial)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

source("src/analysis_config.R")

ASSIGNMENTS_INPUT_PATH <- "output/tabs/clustering/cluster_assignments.csv"
OUTPUT_FIG_DIR <- "output/figs/clustering"
OUTPUT_PDF_PATH <- file.path(OUTPUT_FIG_DIR, "clusters_alluvial.pdf")
REQUIRED_COLUMNS <- c("patient_id", "pam_k3", "pam_k3_tonic_clonic", "pam_k3_focal")
CLUSTER_LEVELS <- c("1", "2", "3")
NO_TONIC_CLONIC_LABEL <- "No\ntonic-\nclonic"
NO_FOCAL_LABEL <- "No\nfocal"
TONIC_CLONIC_LEVELS <- c(CLUSTER_LEVELS, NO_TONIC_CLONIC_LABEL)
FOCAL_LEVELS <- c(CLUSTER_LEVELS, NO_FOCAL_LABEL)

dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggalluvial", quietly = TRUE)) {
  stop("Package 'ggalluvial' is required for the cluster alluvial plot.")
}

if (!file.exists(ASSIGNMENTS_INPUT_PATH)) {
  stop("Assignment input file not found: ", ASSIGNMENTS_INPUT_PATH)
}

assignments <- readr::read_csv(ASSIGNMENTS_INPUT_PATH, show_col_types = FALSE)

missing_columns <- setdiff(REQUIRED_COLUMNS, names(assignments))
if (length(missing_columns) > 0) {
  stop("Assignment input is missing required columns: ", paste(missing_columns, collapse = ", "))
}

alluvial_data <- assignments %>%
  transmute(
    patient_id = as.character(.data$patient_id),
    all = as.character(suppressWarnings(as.integer(.data$pam_k3))),
    `tonic-clonic` = as.character(suppressWarnings(as.integer(.data$pam_k3_tonic_clonic))),
    focal = as.character(suppressWarnings(as.integer(.data$pam_k3_focal)))
  ) %>%
  filter(
    !is.na(.data$patient_id),
    !is.na(.data$all)
  ) %>%
  mutate(
    `tonic-clonic` = if_else(is.na(.data$`tonic-clonic`), NO_TONIC_CLONIC_LABEL, .data$`tonic-clonic`),
    focal = if_else(is.na(.data$focal), NO_FOCAL_LABEL, .data$focal),
    all = factor(.data$all, levels = CLUSTER_LEVELS),
    `tonic-clonic` = factor(.data$`tonic-clonic`, levels = TONIC_CLONIC_LEVELS),
    focal = factor(.data$focal, levels = FOCAL_LEVELS),
    all_cluster = factor(
      as.character(.data$all),
      levels = CLUSTER_LEVELS,
      labels = CLUSTER_LEVELS
    )
  )

if (nrow(alluvial_data) == 0) {
  stop("No patients have non-missing all-cluster assignments.", call. = FALSE)
}

plot_data <- alluvial_data %>%
  tidyr::pivot_longer(
    cols = c("all", "tonic-clonic", "focal"),
    names_to = "cluster_scope",
    values_to = "cluster"
  ) %>%
  mutate(
    cluster_scope = factor(.data$cluster_scope, levels = c("all", "tonic-clonic", "focal"))
  )

y_breaks <- seq(0, ceiling(nrow(alluvial_data) / 5) * 5, by = 5)

p_alluvial <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = .data$cluster_scope,
    stratum = .data$cluster,
    alluvium = .data$patient_id,
    fill = .data$all_cluster
  )
) +
  ggalluvial::geom_flow(width = 0.25, alpha = 0.75, color = "white", linewidth = 0.05) +
  ggalluvial::geom_stratum(width = 0.25, fill = "white", color = "black", linewidth = 0.45) +
  ggplot2::geom_text(
    stat = "stratum",
    ggplot2::aes(label = ggplot2::after_stat(stratum)),
    size = 3.5,
    color = "black"
  ) +
  ggplot2::scale_x_discrete(
    labels = c(
      all = "All",
      `tonic-clonic` = "Tonic-clonic",
      focal = "Focal"
    ),
    expand = c(0.08, 0.08)
  ) +
  ggplot2::scale_y_continuous(breaks = y_breaks, expand = c(0, 0.15)) +
  ggplot2::scale_fill_manual(
    values = EPILEPSIA_CLUSTER_COLORS,
    breaks = names(EPILEPSIA_CLUSTER_COLORS),
    na.translate = FALSE,
    name = "Cluster"
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Number of Patients"
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(color = "grey25"),
    axis.text.y = ggplot2::element_text(color = "grey25"),
    axis.title.y = ggplot2::element_text(color = "black"),
    axis.line.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 13),
    legend.text = ggplot2::element_text(size = 11),
    legend.key.size = grid::unit(0.7, "cm"),
    plot.margin = ggplot2::margin(t = 8, r = 12, b = 8, l = 8)
  )

ggplot2::ggsave(
  filename = OUTPUT_PDF_PATH,
  plot = p_alluvial,
  width = 8,
  height = 5
)
