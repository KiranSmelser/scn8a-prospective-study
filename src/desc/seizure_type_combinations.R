# Top seizure type combinations plot.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

source("src/desc/seizure_type_standardization.R")

OUTPUT_FIG_DIR <- "output/figs/seizure_patterns"
OUTPUT_TAB_DIR <- "output/tabs/seizure_patterns"
dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
png_output_path <- file.path(OUTPUT_FIG_DIR, "seizure_type_combinations.png")
if (file.exists(png_output_path)) {
  invisible(file.remove(png_output_path))
}
old_membership_output_path <- file.path(
  OUTPUT_TAB_DIR,
  "seizure_type_combinations_membership_patient_month.csv"
)
if (file.exists(old_membership_output_path)) {
  invisible(file.remove(old_membership_output_path))
}

analysis_end <- analysis_end_date()
max_sets <- 8L
max_intersections_to_plot <- 25L
epilepsia_teal <- "#5698A3"
publication_base_size <- 16

events <- readr::read_csv("data/events.csv", show_col_types = FALSE) %>%
  standardize_seizure_events(filter_to_seizure = TRUE) %>%
  dplyr::filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    !is.na(.data$patient_id),
    !is.na(.data$month),
    !is.na(.data$event_date),
    .data$event_date <= analysis_end,
    !.data$requires_review
  ) %>%
  dplyr::mutate(seizure_type_combo = .data$seizure_type_primary) %>%
  dplyr::filter(
    !is.na(.data$seizure_type_combo),
    !.data$seizure_type_combo %in% c(
      "Unknown/Unmapped",
      "Unspecified",
      "Device-detected event",
      "Mixed/Multiple"
    )
  ) %>%
  dplyr::distinct(.data$patient_id, .data$month, .data$seizure_type_combo)

if (nrow(events) == 0) {
  warning("No seizure patient-month combinations available after filtering. No outputs generated.")
  quit(save = "no", status = 0)
}

type_patient_month_frequency <- events %>%
  dplyr::count(.data$seizure_type_combo, name = "n_patient_months") %>%
  dplyr::arrange(dplyr::desc(.data$n_patient_months), .data$seizure_type_combo)

type_patient_frequency <- events %>%
  dplyr::distinct(.data$patient_id, .data$seizure_type_combo) %>%
  dplyr::count(.data$seizure_type_combo, name = "n_patients") %>%
  dplyr::arrange(dplyr::desc(.data$n_patients), .data$seizure_type_combo)

top_types <- type_patient_frequency %>%
  dplyr::slice_head(n = max_sets) %>%
  dplyr::pull(.data$seizure_type_combo)

if (length(top_types) < 2) {
  warning("Fewer than 2 seizure types available for combination analysis. No outputs generated.")
  quit(save = "no", status = 0)
}

set_table <- tibble::tibble(
  seizure_type = top_types,
  set_key = paste0("s", seq_along(top_types)),
  y_axis_order = seq_along(top_types)
)

patient_members <- events %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(
    seizure_types = list(sort(unique(.data$seizure_type_combo))),
    .groups = "drop"
  )

upset_membership <- patient_members %>%
  dplyr::select(patient_id)

for (i in seq_len(nrow(set_table))) {
  set_type <- set_table$seizure_type[[i]]
  set_key <- set_table$set_key[[i]]
  upset_membership[[set_key]] <- purrr::map_lgl(
    patient_members$seizure_types,
    ~ set_type %in% .x
  )
}

upset_keys <- set_table$set_key

upset_intersections_all <- upset_membership %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(upset_keys))) %>%
  dplyr::summarise(
    n_patients = dplyr::n(),
    patient_ids = paste(sort(.data$patient_id), collapse = ";"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    n_sets = rowSums(as.data.frame(dplyr::across(dplyr::all_of(upset_keys)))),
    combination = purrr::pmap_chr(
      dplyr::across(dplyr::all_of(upset_keys)),
      function(...) {
        flags <- as.logical(c(...))
        types <- set_table$seizure_type[flags]
        if (length(types) == 0) {
          return("None")
        }
        paste(types, collapse = " + ")
      }
    )
  ) %>%
  dplyr::filter(.data$n_sets > 0) %>%
  dplyr::arrange(dplyr::desc(.data$n_patients), dplyr::desc(.data$n_sets), .data$combination) %>%
  dplyr::mutate(intersection_id = paste0("I", dplyr::row_number()))

upset_intersections_top <- upset_intersections_all %>%
  dplyr::slice_head(n = max_intersections_to_plot) %>%
  dplyr::arrange(dplyr::desc(.data$n_patients), dplyr::desc(.data$n_sets), .data$combination) %>%
  dplyr::mutate(
    plot_order = dplyr::row_number(),
    plot_x = factor(.data$plot_order, levels = .data$plot_order)
  )

n_upset_sets <- nrow(set_table)

upset_matrix <- upset_intersections_top %>%
  dplyr::select(plot_order, plot_x, intersection_id, dplyr::all_of(upset_keys)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(upset_keys),
    names_to = "set_key",
    values_to = "present"
  ) %>%
  dplyr::left_join(set_table, by = "set_key") %>%
  dplyr::mutate(
    y_idx = n_upset_sets - match(.data$seizure_type, top_types) + 1L,
    present = as.logical(.data$present)
  )

upset_segments <- upset_matrix %>%
  dplyr::filter(.data$present) %>%
  dplyr::group_by(.data$plot_order, .data$plot_x, .data$intersection_id) %>%
  dplyr::summarise(
    y_min = min(.data$y_idx),
    y_max = max(.data$y_idx),
    n_present = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(.data$n_present >= 2)

set_patient_totals <- set_table %>%
  dplyr::left_join(type_patient_frequency, by = c("seizure_type" = "seizure_type_combo")) %>%
  dplyr::mutate(
    n_patients = tidyr::replace_na(.data$n_patients, 0L),
    y_idx = n_upset_sets - match(.data$seizure_type, top_types) + 1L
  )

upset_bar_plot <- ggplot2::ggplot(
  upset_intersections_top,
  ggplot2::aes(x = .data$plot_x, y = .data$n_patients)
) +
  ggplot2::geom_col(fill = epilepsia_teal, width = 0.72) +
  ggplot2::geom_text(
    ggplot2::aes(label = .data$n_patients),
    vjust = -0.25,
    size = 4
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.14))
  ) +
  ggplot2::scale_x_discrete(
    limits = as.character(upset_intersections_top$plot_order),
    labels = NULL
  ) +
  ggplot2::labs(
    y = "Number of patients",
    x = NULL
  ) +
  ggplot2::theme_minimal(base_size = publication_base_size) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 13, color = "#222222"),
    axis.title.y = ggplot2::element_text(size = 14),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )

upset_matrix_plot <- ggplot2::ggplot(
  upset_matrix,
  ggplot2::aes(x = .data$plot_x, y = .data$y_idx)
) +
  ggplot2::geom_segment(
    data = upset_segments,
    ggplot2::aes(x = .data$plot_x, xend = .data$plot_x, y = .data$y_min, yend = .data$y_max),
    inherit.aes = FALSE,
    linewidth = 0.5,
    color = "#303030"
  ) +
  ggplot2::geom_point(color = "grey85", size = 2.2) +
  ggplot2::geom_point(
    data = upset_matrix %>% dplyr::filter(.data$present),
    color = epilepsia_teal,
    size = 2.4
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq_len(n_upset_sets),
    labels = rev(top_types),
    limits = c(0.5, n_upset_sets + 0.5),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_x_discrete(
    limits = as.character(upset_intersections_top$plot_order),
    labels = NULL
  ) +
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = publication_base_size) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 13, color = "#222222"),
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

set_total_plot <- ggplot2::ggplot(set_patient_totals) +
  ggplot2::geom_rect(
    ggplot2::aes(
      xmin = 0,
      xmax = .data$n_patients,
      ymin = .data$y_idx - 0.34,
      ymax = .data$y_idx + 0.34
    ),
    fill = epilepsia_teal
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = .data$n_patients, y = .data$y_idx, label = .data$n_patients),
    hjust = -0.25,
    size = 4
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq_len(n_upset_sets),
    labels = NULL,
    limits = c(0.5, n_upset_sets + 0.5),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.18))
  ) +
  ggplot2::labs(
    x = "Number of patients",
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = publication_base_size) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 13, color = "#222222"),
    axis.title.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

upset_plot <- patchwork::wrap_plots(
  upset_bar_plot,
  patchwork::plot_spacer(),
  upset_matrix_plot,
  set_total_plot,
  ncol = 2,
  widths = c(4.2, 1.35),
  heights = c(1.25, 3.1)
)

ggplot2::ggsave(
  filename = file.path(OUTPUT_FIG_DIR, "seizure_type_combinations.pdf"),
  plot = upset_plot,
  width = max(11, 0.33 * nrow(upset_intersections_top) + 4),
  height = max(6, 0.25 * n_upset_sets + 3)
)

readr::write_csv(
  type_patient_month_frequency,
  file.path(OUTPUT_TAB_DIR, "seizure_type_patient_month_frequency.csv")
)
readr::write_csv(
  type_patient_frequency,
  file.path(OUTPUT_TAB_DIR, "seizure_type_patient_frequency.csv")
)
readr::write_csv(
  set_table,
  file.path(OUTPUT_TAB_DIR, "seizure_type_combinations_sets.csv")
)
readr::write_csv(
  upset_membership,
  file.path(OUTPUT_TAB_DIR, "seizure_type_combinations_membership_patient.csv")
)
readr::write_csv(
  set_patient_totals %>%
    dplyr::select(set_key, seizure_type, y_idx, n_patients),
  file.path(OUTPUT_TAB_DIR, "seizure_type_combinations_set_patient_totals.csv")
)
readr::write_csv(
  upset_intersections_all %>%
    dplyr::select(
      intersection_id,
      n_patients,
      n_sets,
      combination,
      dplyr::all_of(upset_keys),
      patient_ids
    ),
  file.path(OUTPUT_TAB_DIR, "seizure_type_combinations_intersections_all.csv")
)
readr::write_csv(
  upset_intersections_top %>%
    dplyr::select(
      plot_order,
      intersection_id,
      n_patients,
      n_sets,
      combination,
      dplyr::all_of(upset_keys),
      patient_ids
    ),
  file.path(OUTPUT_TAB_DIR, "seizure_type_combinations_intersections_top.csv")
)
readr::write_csv(
  upset_matrix %>%
    dplyr::select(plot_order, intersection_id, set_key, seizure_type, y_idx, present) %>%
    dplyr::arrange(.data$plot_order, .data$y_idx),
  file.path(OUTPUT_TAB_DIR, "seizure_type_combinations_matrix.csv")
)
