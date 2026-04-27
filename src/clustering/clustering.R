suppressPackageStartupMessages({
  library(cluster)
  library(dplyr)
  library(readr)
})

INPUT_PATH <- "output/tabs/clustering/seizure_freq_features.csv"
OUTPUT_DIR <- "output/tabs/clustering"
ASSIGNMENTS_OUTPUT_PATH <- file.path(OUTPUT_DIR, "cluster_assignments.csv")
CLUSTER_COUNT <- 3L
RATE_FEATURES <- c(
  "mean_monthly_seizure_rate",
  "iqr_monthly_seizure_rate"
)
DIRECT_SCALE_FEATURES <- c("proportion_zero_seizure_months")
REQUIRED_COLUMNS <- c("patient_id", "variant_p", RATE_FEATURES, DIRECT_SCALE_FEATURES)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_PATH)) {
  stop("Input file not found: ", INPUT_PATH)
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

feature_table <- readr::read_csv(INPUT_PATH, show_col_types = FALSE)

missing_columns <- setdiff(REQUIRED_COLUMNS, names(feature_table))
if (length(missing_columns) > 0) {
  stop("Input data is missing required columns: ", paste(missing_columns, collapse = ", "))
}

clustering_input <- feature_table %>%
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
PCA_COMPONENT_COUNT <- 2L

if (nrow(clustering_input) < CLUSTER_COUNT) {
  stop(
    "Not enough patients for requested cluster count. Patients available: ",
    nrow(clustering_input),
    "; requested k: ",
    CLUSTER_COUNT
  )
}

pca_input_matrix <- clustering_input %>%
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

pc_feature_names <- paste0("pc", seq_len(PCA_COMPONENT_COUNT))
pc_feature_matrix <- pca_fit$x[, seq_len(PCA_COMPONENT_COUNT), drop = FALSE]
colnames(pc_feature_matrix) <- pc_feature_names

ward_fit <- stats::hclust(stats::dist(pc_feature_matrix), method = "ward.D2")

ward_assignments <- tibble::tibble(
  ward_k3 = stats::cutree(ward_fit, k = CLUSTER_COUNT)
)

pam_fit <- cluster::pam(pc_feature_matrix, k = CLUSTER_COUNT, metric = "euclidean", stand = FALSE)
pam_assignments <- tibble::tibble(
  pam_k3 = pam_fit$clustering
)

cluster_assignments <- clustering_input %>%
  select("patient_id", "variant_p", all_of(RATE_FEATURES), all_of(DIRECT_SCALE_FEATURES)) %>%
  bind_cols(ward_assignments, pam_assignments) %>%
  arrange(.data$patient_id)

readr::write_csv(cluster_assignments, ASSIGNMENTS_OUTPUT_PATH)
