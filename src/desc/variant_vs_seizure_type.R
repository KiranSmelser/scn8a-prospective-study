# Variant vs. seizure type analysis.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

source("src/desc/seizure_type_standardization.R")

OUTPUT_FIG_DIR <- "output/figs/seizure_patterns"
OUTPUT_TAB_DIR <- "output/tabs/seizure_patterns"
dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
png_output_path <- file.path(OUTPUT_FIG_DIR, "variant_vs_seizure_type.png")
if (file.exists(png_output_path)) {
  invisible(file.remove(png_output_path))
}

analysis_end <- analysis_end_date()

events <- readr::read_csv("data/events.csv", show_col_types = FALSE) %>%
  standardize_seizure_events(filter_to_seizure = TRUE) %>%
  dplyr::filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    !is.na(.data$patient_id),
    !is.na(.data$event_date),
    .data$event_date <= analysis_end,
    !.data$requires_review
  ) %>%
  dplyr::mutate(seizure_type_plot = .data$seizure_type_primary) %>%
  dplyr::filter(
    !is.na(.data$seizure_type_plot),
    !.data$seizure_type_plot %in% c(
      "Unknown/Unmapped",
      "Unspecified",
      "Device-detected event",
      "Mixed/Multiple"
    )
  )

if (nrow(events) == 0) {
  warning("No reviewed seizure events available. No outputs generated.")
  quit(save = "no", status = 0)
}

variants <- readr::read_csv("data/whatsapp_status.csv", show_col_types = FALSE) %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    status_timestamp_utc = suppressWarnings(
      lubridate::ymd_hms(stringr::str_replace_all(as.character(.data$run_timestamp), "_", " "), quiet = TRUE, tz = "UTC")
    ),
    variant = stringr::str_squish(as.character(.data$variant_p))
  ) %>%
  dplyr::mutate(
    variant = dplyr::na_if(.data$variant, ""),
    variant = dplyr::if_else(
      !is.na(.data$variant) & stringr::str_to_lower(.data$variant) %in% c("na", "n/a", "nan", "unknown"),
      NA_character_,
      .data$variant
    ),
    status_timestamp_utc = dplyr::coalesce(
      .data$status_timestamp_utc,
      as.POSIXct("1900-01-01 00:00:00", tz = "UTC")
    )
  ) %>%
  dplyr::filter(.data$patient_id %in% TARGET_PATIENT_IDS, !is.na(.data$patient_id)) %>%
  dplyr::arrange(.data$patient_id, dplyr::desc(.data$status_timestamp_utc)) %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(variant = dplyr::first(.data$variant), .groups = "drop")

events_with_variant <- events %>%
  dplyr::select("patient_id", "seizure_type_plot") %>%
  dplyr::left_join(variants, by = "patient_id") %>%
  dplyr::mutate(variant = dplyr::coalesce(.data$variant, "Unknown")) %>%
  dplyr::filter(.data$variant != "Unknown")

if (nrow(events_with_variant) == 0) {
  warning("No reviewed seizure events with known variants available. No outputs generated.")
  quit(save = "no", status = 0)
}

variant_patient_counts <- events_with_variant %>%
  dplyr::distinct(.data$patient_id, .data$variant) %>%
  dplyr::count(.data$variant, name = "n_unique_patients") %>%
  dplyr::mutate(variant_label = paste0(.data$variant, "\n(n=", .data$n_unique_patients, ")"))

base_counts <- events_with_variant %>%
  dplyr::count(.data$variant, .data$seizure_type_plot, name = "n_events")

type_order <- base_counts %>%
  dplyr::group_by(.data$seizure_type_plot) %>%
  dplyr::summarise(n_events_total = sum(.data$n_events), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(.data$n_events_total), .data$seizure_type_plot) %>%
  dplyr::pull(.data$seizure_type_plot)

variant_order <- base_counts %>%
  dplyr::group_by(.data$variant) %>%
  dplyr::summarise(n_events_variant = sum(.data$n_events), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(.data$n_events_variant), .data$variant) %>%
  dplyr::pull(.data$variant)

variant_label_lookup <- stats::setNames(variant_patient_counts$variant_label, variant_patient_counts$variant)
variant_levels <- rev(unname(variant_label_lookup[variant_order]))

heatmap_counts <- tidyr::expand_grid(
  variant = variant_order,
  seizure_type_plot = type_order
) %>%
  dplyr::left_join(base_counts, by = c("variant", "seizure_type_plot")) %>%
  tidyr::replace_na(list(n_events = 0L)) %>%
  dplyr::left_join(
    variant_patient_counts %>% dplyr::select("variant", "variant_label"),
    by = "variant"
  ) %>%
  dplyr::mutate(
    variant = factor(.data$variant_label, levels = variant_levels),
    seizure_type_plot = factor(.data$seizure_type_plot, levels = type_order)
  )

max_events <- max(heatmap_counts$n_events, na.rm = TRUE)
if (!is.finite(max_events)) {
  max_events <- 0
}

heatmap_counts <- heatmap_counts %>%
  dplyr::mutate(
    text_color = dplyr::if_else(
      max_events > 0 & .data$n_events >= 0.55 * max_events,
      "white",
      "black"
    )
  )

p <- ggplot2::ggplot(
  heatmap_counts,
  ggplot2::aes(x = .data$seizure_type_plot, y = .data$variant, fill = .data$n_events)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::geom_text(
    ggplot2::aes(label = .data$n_events, color = .data$text_color),
    size = 3
  ) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_fill_gradient(
    low = "#FDECEC",
    high = "#C62828",
    name = "# of seizures"
  ) +
  ggplot2::labs(
    title = "Variants vs Seizure Types",
    x = "Seizure type",
    y = "Variant"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold")
  )

plot_width <- max(8, 1.0 * dplyr::n_distinct(heatmap_counts$seizure_type_plot) + 3)
plot_height <- max(5, 0.35 * dplyr::n_distinct(heatmap_counts$variant) + 2.5)

ggplot2::ggsave(
  filename = file.path(OUTPUT_FIG_DIR, "variant_vs_seizure_type.pdf"),
  plot = p,
  width = plot_width,
  height = plot_height
)

readr::write_csv(
  heatmap_counts %>%
    dplyr::mutate(
      variant = as.character(.data$variant),
      seizure_type_plot = as.character(.data$seizure_type_plot)
    ) %>%
    dplyr::select("variant", "seizure_type_plot", "n_events") %>%
    dplyr::arrange(.data$variant, .data$seizure_type_plot),
  file.path(OUTPUT_TAB_DIR, "variant_vs_seizure_type_event_counts.csv")
)

variant_event_totals <- base_counts %>%
  dplyr::group_by(.data$variant) %>%
  dplyr::summarise(n_events_with_variant = sum(.data$n_events), .groups = "drop") %>%
  dplyr::left_join(
    variant_patient_counts %>% dplyr::select("variant", "n_unique_patients"),
    by = "variant"
  ) %>%
  dplyr::arrange(dplyr::desc(.data$n_events_with_variant), .data$variant)

readr::write_csv(
  variant_event_totals,
  file.path(OUTPUT_TAB_DIR, "variant_event_totals.csv")
)
