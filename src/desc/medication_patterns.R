# Active/weaned medication summaries and combination patterns.

suppressPackageStartupMessages({
  library(tidyverse)
})

source("src/analysis_config.R")

INPUT_PATH <- "output/tabs/medication_standardization/non_rescue_epilepsy_medications_standardized.csv"
WHATSAPP_STATUS_PATH <- "data/whatsapp_status.csv"
OUTPUT_TAB_DIR <- "output/tabs/medication_patterns"
OUTPUT_FIG_DIR <- "output/figs/medication_patterns"

dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

legacy_png_paths <- file.path(
  OUTPUT_FIG_DIR,
  c(
    "active_vs_weaned_medications_all.png",
    "active_medication_pairwise_cooccurrence_heatmap.png",
    "active_medication_combinations_upset.png",
    "patient_level_medication_burden_histogram.png"
  )
)
invisible(file.remove(legacy_png_paths[file.exists(legacy_png_paths)]))

if (!file.exists(INPUT_PATH)) {
  stop("Missing input file: ", INPUT_PATH, call. = FALSE)
}

pct <- function(numerator, denominator) {
  if (!is.finite(denominator) || denominator <= 0) {
    return(rep(NA_real_, length(numerator)))
  }
  100 * numerator / denominator
}

build_combinations <- function(meds, min_size = 2L, max_size = 4L) {
  med_vec <- sort(unique(as.character(meds)))
  med_vec <- med_vec[!is.na(med_vec) & nzchar(med_vec)]
  n_meds <- length(med_vec)
  if (n_meds < min_size) {
    return(character())
  }

  sizes <- seq(min_size, min(max_size, n_meds))
  combos <- purrr::map(
    sizes,
    function(k) {
      combn(med_vec, k, FUN = function(x) paste(x, collapse = " + "))
    }
  ) %>%
    unlist(use.names = FALSE)

  sort(unique(combos))
}

load_patient_mutation_lookup <- function(path) {
  if (!file.exists(path)) {
    return(setNames(character(), character()))
  }

  whatsapp_status <- readr::read_csv(path, show_col_types = FALSE)
  required_cols <- c("patient_id", "variant_p")
  if (!all(required_cols %in% names(whatsapp_status))) {
    return(setNames(character(), character()))
  }

  lookup <- whatsapp_status %>%
    transmute(
      patient_id = stringr::str_squish(as.character(patient_id)),
      mutation = stringr::str_squish(as.character(variant_p)),
      mutation = na_if(mutation, ""),
      mutation = if_else(
        !is.na(mutation) & stringr::str_to_lower(mutation) %in% c("na", "n/a", "nan", "unknown"),
        NA_character_,
        mutation
      )
    ) %>%
    filter(!is.na(patient_id) & nzchar(patient_id)) %>%
    group_by(patient_id) %>%
    summarise(
      mutation_label = {
        values <- sort(unique(na.omit(mutation)))
        if (length(values) == 0) "NA" else paste(values, collapse = " / ")
      },
      .groups = "drop"
    )

  setNames(lookup$mutation_label, lookup$patient_id)
}

intersection_variant_label <- function(patient_ids, mutation_lookup) {
  ids <- stringr::str_split(as.character(patient_ids), ";", simplify = FALSE)[[1]]
  ids <- stringr::str_squish(ids)
  ids <- ids[nzchar(ids)]

  if (length(ids) == 0) {
    return("NA")
  }

  variants <- unname(mutation_lookup[ids])
  variants[is.na(variants) | !nzchar(variants)] <- "NA"
  variant_counts <- table(variants)
  variant_labels <- names(variant_counts)
  duplicate_variants <- variant_counts > 1
  variant_labels[duplicate_variants] <- paste0(
    variant_labels[duplicate_variants],
    " (n = ",
    variant_counts[duplicate_variants],
    ")"
  )
  paste(variant_labels, collapse = ", ")
}

meds <- readr::read_csv(INPUT_PATH, show_col_types = FALSE) %>%
  mutate(
    patient_id = stringr::str_squish(as.character(patient_id)),
    name_standardized = if_else(
      is.na(name_standardized) | !nzchar(trimws(name_standardized)),
      "Unknown medication",
      trimws(name_standardized)
    ),
    is_active_on_analysis_date = as.logical(is_active_on_analysis_date),
    requires_review = as.logical(requires_review),
    is_active_on_analysis_date = tidyr::replace_na(is_active_on_analysis_date, FALSE),
    requires_review = tidyr::replace_na(requires_review, FALSE)
  ) %>%
  filter(.data$patient_id %in% TARGET_PATIENT_IDS)

# Exclude medications still requiring review from all downstream analyses.
meds <- meds %>%
  filter(!requires_review)

patient_mutation_lookup <- load_patient_mutation_lookup(WHATSAPP_STATUS_PATH)

patients_total <- meds %>% distinct(patient_id) %>% nrow()

# Collapse to patient-medication status to avoid double counting medication IDs.
patient_med_status <- meds %>%
  group_by(patient_id, name_standardized) %>%
  summarise(
    n_medication_records = n_distinct(medication_id),
    is_active = any(is_active_on_analysis_date),
    is_weaned = !any(is_active_on_analysis_date),
    any_review_flag = any(requires_review),
    .groups = "drop"
  )

# 1) Active Medications Frequency Table
active_medication_frequency <- patient_med_status %>%
  filter(is_active) %>%
  group_by(name_standardized) %>%
  summarise(
    n_patients_active = n_distinct(patient_id),
    pct_patients_active = round(pct(n_patients_active, patients_total), 1),
    n_active_patient_med_rows = n(),
    n_active_rows_with_review_flag = sum(any_review_flag),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients_active), name_standardized)

readr::write_csv(
  active_medication_frequency,
  file.path(OUTPUT_TAB_DIR, "active_medication_frequency.csv")
)

# 2) Weaned Medications Frequency Table
weaned_medication_frequency <- patient_med_status %>%
  filter(is_weaned) %>%
  group_by(name_standardized) %>%
  summarise(
    n_patients_weaned = n_distinct(patient_id),
    pct_patients_weaned = round(pct(n_patients_weaned, patients_total), 1),
    n_weaned_patient_med_rows = n(),
    n_weaned_rows_with_review_flag = sum(any_review_flag),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients_weaned), name_standardized)

readr::write_csv(
  weaned_medication_frequency,
  file.path(OUTPUT_TAB_DIR, "weaned_medication_frequency.csv")
)

# 3) Active vs Weaned Comparison Plot
comparison_long <- bind_rows(
  active_medication_frequency %>%
    transmute(
      name_standardized,
      status = "Active",
      n_patients = n_patients_active
    ),
  weaned_medication_frequency %>%
    transmute(
      name_standardized,
      status = "Weaned",
      n_patients = n_patients_weaned
    )
) %>%
  group_by(name_standardized, status) %>%
  summarise(n_patients = sum(n_patients), .groups = "drop")

comparison_med_totals <- comparison_long %>%
  group_by(name_standardized) %>%
  summarise(total_n_patients = sum(n_patients), .groups = "drop") %>%
  arrange(desc(total_n_patients), name_standardized)

comparison_plot_data <- comparison_long %>%
  tidyr::complete(name_standardized, status = c("Active", "Weaned"), fill = list(n_patients = 0)) %>%
  left_join(comparison_med_totals, by = "name_standardized") %>%
  mutate(name_standardized = forcats::fct_reorder(name_standardized, total_n_patients))

comparison_plot <- ggplot(
  comparison_plot_data,
  aes(x = n_patients, y = name_standardized, fill = status)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_fill_manual(values = c("Active" = "#1f78b4", "Weaned" = "#e31a1c")) +
  labs(
    title = "Active vs Weaned Medications",
    x = "Number of Patients",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "top"
  )

ggplot2::ggsave(
  file.path(OUTPUT_FIG_DIR, "active_vs_weaned_medications_all.pdf"),
  plot = comparison_plot,
  width = 10,
  height = 7
)

readr::write_csv(
  comparison_plot_data %>%
    mutate(name_standardized = as.character(name_standardized)) %>%
    arrange(desc(total_n_patients), name_standardized, status),
  file.path(OUTPUT_TAB_DIR, "active_vs_weaned_plot_data_all.csv")
)

# 4) Most Common Active Medication Combinations
active_med_list_by_patient <- patient_med_status %>%
  filter(is_active) %>%
  group_by(patient_id) %>%
  summarise(
    active_medications = list(sort(unique(name_standardized))),
    n_active_medications = lengths(active_medications),
    .groups = "drop"
  )

n_patients_with_2plus_active <- active_med_list_by_patient %>%
  filter(n_active_medications >= 2) %>%
  nrow()

combination_table <- active_med_list_by_patient %>%
  mutate(combination = purrr::map(active_medications, build_combinations, min_size = 2L, max_size = 4L)) %>%
  select(patient_id, n_active_medications, combination) %>%
  tidyr::unnest(combination) %>%
  distinct(patient_id, combination, .keep_all = TRUE) %>%
  mutate(combination_size = stringr::str_count(combination, " \\+ ") + 1L) %>%
  group_by(combination, combination_size) %>%
  summarise(
    n_patients = n_distinct(patient_id),
    pct_patients_all = round(pct(n_patients, patients_total), 1),
    pct_patients_with_2plus_active = round(pct(n_patients, n_patients_with_2plus_active), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients), combination_size, combination) %>%
  mutate(rank = row_number()) %>%
  select(rank, combination, combination_size, n_patients, pct_patients_all, pct_patients_with_2plus_active)

readr::write_csv(
  combination_table,
  file.path(OUTPUT_TAB_DIR, "most_common_active_medication_combinations_size2to4.csv")
)

readr::write_csv(
  combination_table %>% slice_head(n = 30),
  file.path(OUTPUT_TAB_DIR, "most_common_active_medication_combinations_top30.csv")
)

# 5) Pairwise Co-Occurrence Heatmap (active medications)
active_medications_order <- active_medication_frequency %>%
  arrange(desc(n_patients_active), name_standardized) %>%
  pull(name_standardized)

incidence <- patient_med_status %>%
  filter(is_active) %>%
  transmute(
    patient_id,
    medication = factor(name_standardized, levels = active_medications_order),
    present = 1L
  ) %>%
  tidyr::pivot_wider(
    names_from = medication,
    values_from = present,
    values_fill = 0
  )

if (nrow(incidence) > 0 && ncol(incidence) > 2) {
  incidence_matrix <- incidence %>%
    select(-patient_id) %>%
    as.matrix()
  storage.mode(incidence_matrix) <- "integer"

  co_counts <- t(incidence_matrix) %*% incidence_matrix
  med_names <- colnames(incidence_matrix)
  n_patients_active_any <- nrow(incidence_matrix)
  med_totals <- colSums(incidence_matrix)

  pairwise_full_raw <- as.data.frame(as.table(co_counts), stringsAsFactors = FALSE)
  colnames(pairwise_full_raw) <- c("med_a", "med_b", "n_patients_both")
  pairwise_full <- as_tibble(pairwise_full_raw) %>%
    mutate(
      med_a = as.character(med_a),
      med_b = as.character(med_b),
      n_patients_a = med_totals[med_a],
      n_patients_b = med_totals[med_b],
      pct_patients_all = pct(n_patients_both, n_patients_active_any),
      expected_count = (n_patients_a * n_patients_b) / n_patients_active_any,
      lift = if_else(expected_count > 0, n_patients_both / expected_count, NA_real_),
      is_diagonal = med_a == med_b
    )

  pairwise_unique <- pairwise_full %>%
    filter(med_a < med_b) %>%
    arrange(desc(n_patients_both), desc(lift), med_a, med_b)

  readr::write_csv(
    pairwise_unique %>%
      mutate(
        pct_patients_all = round(pct_patients_all, 1),
        lift = round(lift, 3)
      ),
    file.path(OUTPUT_TAB_DIR, "active_medication_pairwise_cooccurrence.csv")
  )

  pairwise_heatmap_data <- pairwise_full %>%
    mutate(
      med_a_idx = match(med_a, active_medications_order),
      med_b_idx = match(med_b, active_medications_order)
    ) %>%
    filter(med_a_idx >= med_b_idx) %>%
    mutate(
      med_a = factor(med_a, levels = active_medications_order),
      med_b = factor(med_b, levels = rev(active_medications_order))
    )

  pairwise_heatmap <- ggplot(
    pairwise_heatmap_data,
    aes(x = med_a, y = med_b, fill = n_patients_both)
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = n_patients_both), size = 2.6, color = "black") +
    scale_fill_gradient(
      low = "#f7fbff",
      high = "#08519c",
      name = "Patients\n(both active)"
    ) +
    labs(
      title = "Active Medication Pairwise Co-Occurrences",
      x = NULL,
      y = NULL
    ) +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "#d9d9d9", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 9),
      legend.position = "none"
    )

  ggplot2::ggsave(
    file.path(OUTPUT_FIG_DIR, "active_medication_pairwise_cooccurrence_heatmap.pdf"),
    plot = pairwise_heatmap,
    width = 11,
    height = 10
  )

  readr::write_csv(
    pairwise_heatmap_data %>%
      mutate(
        med_a = as.character(med_a),
        med_b = as.character(med_b),
        pct_patients_all = round(pct_patients_all, 1),
        lift = round(lift, 3)
      ) %>%
      arrange(med_a_idx, med_b_idx),
    file.path(OUTPUT_TAB_DIR, "active_medication_pairwise_heatmap_data.csv")
  )
}

# 6) Active Medication Combination UpSet Plot
top_upset_meds <- active_medication_frequency %>%
  arrange(desc(n_patients_active), name_standardized) %>%
  pull(name_standardized)

medication_guide_category_order <- c(
  "Sodium Channel Blockers",
  "GABAergic",
  "Calcium Channel Blockers",
  "SVP2A",
  "Other/Multiple",
  "Steroids"
)

upset_category_palette <- c(
  "Sodium Channel Blockers" = "#1b9e77",
  "GABAergic" = "#d95f02",
  "Calcium Channel Blockers" = "#7570b3",
  "SVP2A" = "#e7298a",
  "Other/Multiple" = "#66a61e",
  "Steroids" = "#e6ab02",
  "Uncategorized" = "#666666"
)

medication_guide_map <- tribble(
  ~medication_normalized, ~guide_category, ~guide_med_rank, ~ordering_source,
  "carbamazepine", "Sodium Channel Blockers", 1L, "SCN8A medication guide",
  "lacosamide", "Sodium Channel Blockers", 2L, "SCN8A medication guide",
  "lamotrigine", "Sodium Channel Blockers", 3L, "SCN8A medication guide",
  "oxcarbazepine", "Sodium Channel Blockers", 4L, "SCN8A medication guide",
  "phenytoin", "Sodium Channel Blockers", 5L, "SCN8A medication guide",
  "rufinamide", "Sodium Channel Blockers", 6L, "SCN8A medication guide",
  "valproate", "Sodium Channel Blockers", 7L, "SCN8A medication guide",
  "valproic acid", "Sodium Channel Blockers", 7L, "SCN8A medication guide",
  "clobazam", "GABAergic", 1L, "SCN8A medication guide",
  "clonazepam", "GABAergic", 2L, "SCN8A medication guide",
  "felbamate", "GABAergic", 3L, "SCN8A medication guide",
  "phenobarbital", "GABAergic", 5L, "SCN8A medication guide",
  "vigabatrin", "GABAergic", 9L, "SCN8A medication guide",
  "ethosuximide", "Calcium Channel Blockers", 1L, "SCN8A medication guide",
  "gabapentin", "Calcium Channel Blockers", 2L, "SCN8A medication guide",
  "zonisamide", "Calcium Channel Blockers", 3L, "SCN8A medication guide",
  "brivaracetam", "SVP2A", 1L, "SCN8A medication guide",
  "levetiracetam", "SVP2A", 2L, "SCN8A medication guide",
  "topiramate", "Other/Multiple", 5L, "SCN8A medication guide",
  "cenobamate", "Other/Multiple", 6L, "SCN8A medication guide"
)

# Fallback for medications not listed in the SCN8A guide.
mechanism_fallback_map <- tribble(
  ~medication_normalized, ~guide_category, ~guide_med_rank, ~ordering_source,
  "praxis", "Sodium Channel Blockers", 8L, "Mechanism fallback: Relutrigine (PRAX-562) persistent sodium current inhibitor",
  "cannabidiol", "Other/Multiple", 7L, "Mechanism fallback: EPIDIOLEX label reports anticonvulsant mechanism is unknown"
)

  upset_medication_order <- tibble(
  medication = top_upset_meds,
  medication_normalized = stringr::str_to_lower(stringr::str_squish(top_upset_meds))
) %>%
  left_join(medication_guide_map, by = "medication_normalized") %>%
  left_join(mechanism_fallback_map, by = "medication_normalized", suffix = c("_guide", "_fallback")) %>%
  mutate(
    guide_category = coalesce(guide_category_guide, guide_category_fallback, "Other/Multiple"),
    guide_med_rank = coalesce(guide_med_rank_guide, guide_med_rank_fallback, 999L),
    ordering_source = coalesce(
      ordering_source_guide,
      ordering_source_fallback,
      "Mechanism fallback: uncategorized"
    ),
    guide_category_rank = match(guide_category, medication_guide_category_order),
    guide_category_rank = if_else(
      is.na(guide_category_rank),
      length(medication_guide_category_order) + 1L,
      guide_category_rank
    )
  ) %>%
  arrange(guide_category_rank, guide_med_rank, medication) %>%
  mutate(y_axis_order = row_number()) %>%
  select(y_axis_order, medication, guide_category, guide_med_rank, ordering_source)

  top_upset_meds <- upset_medication_order$medication
  n_upset_meds <- length(top_upset_meds)
  upset_category_levels_present <- upset_medication_order %>%
    distinct(guide_category) %>%
    pull(guide_category)

  upset_set_table <- upset_medication_order %>%
  transmute(
    med_key = paste0("m", row_number()),
    medication,
    guide_category,
    guide_med_rank,
    ordering_source,
    y_axis_order
  )

if (length(top_upset_meds) >= 2 && nrow(active_med_list_by_patient) > 0) {
  upset_membership <- tibble(patient_id = active_med_list_by_patient$patient_id)
  for (i in seq_along(top_upset_meds)) {
    med_name <- top_upset_meds[[i]]
    med_key <- upset_set_table$med_key[[i]]
    upset_membership[[med_key]] <- purrr::map_lgl(
      active_med_list_by_patient$active_medications,
      ~ med_name %in% .x
    )
  }

  upset_keys <- upset_set_table$med_key

  upset_intersections <- upset_membership %>%
    group_by(across(all_of(upset_keys))) %>%
    summarise(
      n_patients = n(),
      patient_ids = paste(sort(unique(patient_id)), collapse = ";"),
      .groups = "drop"
    ) %>%
    mutate(
      n_sets = rowSums(as.data.frame(across(all_of(upset_keys)))),
      combination = purrr::pmap_chr(
        across(all_of(upset_keys)),
        function(...) {
          flags <- as.logical(c(...))
          meds <- upset_set_table$medication[flags]
          if (length(meds) == 0) "None" else paste(meds, collapse = " + ")
        }
      )
    ) %>%
    filter(n_sets > 0) %>%
    arrange(desc(n_patients), desc(n_sets), combination) %>%
    mutate(
      intersection_id = paste0("I", row_number()),
      variant_label = purrr::map_chr(
        patient_ids,
        intersection_variant_label,
        mutation_lookup = patient_mutation_lookup
      )
    )

  upset_n_intersections <- nrow(upset_intersections)

  upset_plot_intersections <- upset_intersections %>%
    slice_head(n = upset_n_intersections) %>%
    mutate(
      plot_order = row_number(),
      plot_intersection_id = paste0("P", plot_order)
    )
  plot_levels <- upset_plot_intersections$plot_order
  upset_plot_intersections <- upset_plot_intersections %>%
    mutate(plot_x = factor(plot_order, levels = plot_levels))

  upset_matrix <- upset_plot_intersections %>%
    select(plot_order, plot_x, intersection_id, all_of(upset_keys)) %>%
    pivot_longer(
      cols = all_of(upset_keys),
      names_to = "med_key",
      values_to = "present"
    ) %>%
    left_join(upset_set_table, by = "med_key") %>%
    mutate(
      y_idx = n_upset_meds - match(medication, top_upset_meds) + 1L,
      present = as.logical(present)
    )

  upset_segments <- upset_matrix %>%
    filter(present) %>%
    group_by(plot_order, plot_x, intersection_id) %>%
    summarise(
      y_min = min(y_idx),
      y_max = max(y_idx),
      n_present = n(),
      .groups = "drop"
    ) %>%
    filter(n_present >= 2)

  upset_matrix_plot <- ggplot(
    upset_matrix,
    aes(x = plot_x, y = y_idx)
  ) +
    geom_segment(
      data = upset_segments,
      aes(x = plot_x, xend = plot_x, y = y_min, yend = y_max),
      inherit.aes = FALSE,
      linewidth = 0.5,
      color = "black"
    ) +
    geom_point(color = "grey85", size = 2.2) +
    geom_point(
      data = upset_matrix %>% filter(present),
      aes(color = guide_category),
      size = 2.4
    ) +
    scale_y_continuous(
      breaks = seq_len(n_upset_meds),
      labels = rev(top_upset_meds)
    ) +
    scale_color_manual(
      values = upset_category_palette,
      breaks = upset_category_levels_present,
      drop = FALSE
    ) +
    scale_x_discrete(
      limits = as.character(plot_levels),
      labels = setNames(
        upset_plot_intersections$variant_label,
        as.character(upset_plot_intersections$plot_x)
      )
    ) +
    labs(
      title = "Active Medication Combinations UpSet Plot",
      x = "SCN8A Variant(s)",
      y = NULL,
      color = "Category"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

  upset_plot_path <- file.path(OUTPUT_FIG_DIR, "active_medication_combinations_upset.pdf")
  ggplot2::ggsave(
    filename = upset_plot_path,
    plot = upset_matrix_plot,
    width = max(10, 0.3 * upset_n_intersections + 3),
    height = max(6, 0.25 * length(top_upset_meds) + 3)
  )

  readr::write_csv(
    upset_set_table,
    file.path(OUTPUT_TAB_DIR, "active_medication_upset_sets.csv")
  )
  readr::write_csv(
    upset_medication_order,
    file.path(OUTPUT_TAB_DIR, "active_medication_upset_yaxis_order.csv")
  )
  readr::write_csv(
    upset_intersections %>%
      select(intersection_id, n_patients, n_sets, combination, variant_label, all_of(upset_keys), patient_ids),
    file.path(OUTPUT_TAB_DIR, "active_medication_upset_intersections_all.csv")
  )
  readr::write_csv(
    upset_plot_intersections %>%
      select(
        plot_order, plot_intersection_id, intersection_id, n_patients, n_sets,
        combination, variant_label, all_of(upset_keys), patient_ids
      ),
    file.path(OUTPUT_TAB_DIR, "active_medication_upset_intersections.csv")
  )
  readr::write_csv(
    upset_matrix %>%
      left_join(
        upset_plot_intersections %>%
          select(plot_order, plot_intersection_id, intersection_id),
        by = c("plot_order", "intersection_id")
      ) %>%
      select(plot_order, plot_intersection_id, intersection_id, med_key, medication, y_idx, present) %>%
      arrange(plot_order, y_idx),
    file.path(OUTPUT_TAB_DIR, "active_medication_upset_matrix.csv")
  )
}

# Remove outputs from previous runs.
burden_tab_path <- file.path(OUTPUT_TAB_DIR, "patient_level_medication_burden.csv")
burden_fig_path <- file.path(OUTPUT_FIG_DIR, "patient_level_medication_burden_histogram.pdf")
if (file.exists(burden_tab_path)) {
  invisible(file.remove(burden_tab_path))
}
if (file.exists(burden_fig_path)) {
  invisible(file.remove(burden_fig_path))
}
