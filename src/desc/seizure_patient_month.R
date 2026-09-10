# Patient-month stacked seizure counts by type for all active patients.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(jsonlite)
})

source("src/desc/seizure_type_standardization.R")

OUTPUT_FIG_DIR <- "output/figs/seizure_patterns"
OUTPUT_TAB_DIR <- "output/tabs/seizure_patterns"
dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
png_output_path <- file.path(OUTPUT_FIG_DIR, "seizure_patient_month.png")
if (file.exists(png_output_path)) {
  invisible(file.remove(png_output_path))
}

analysis_end <- last_complete_analysis_month_end()
analysis_month_cap <- as.Date(lubridate::floor_date(analysis_end, "month"))

SEIZURE_TYPE_COLORS <- c(
  "Tonic-clonic" = "#5698a3",
  "Focal" = "#ffde76",
  "Tonic" = "#67771a",
  "Myoclonic" = "#0076c0",
  "Absence" = "#e37c1d",
  "Spasms" = "#7a5072",
  "Status epilepticus" = "#B22222",
  "Clonic" = "#3B6FB6",
  "Generalized motor (non tonic-clonic)" = "#8C564B",
  "Atonic" = "#BCBD22",
  "Other" = "#B0B0B0"
)

all_seizure_events <- read_events_corrected() %>%
  standardize_seizure_events(filter_to_seizure = TRUE) %>%
  dplyr::filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    !is.na(.data$patient_id),
    !is.na(.data$month),
    !is.na(.data$event_date),
    .data$event_date <= analysis_end
  )

whatsapp_names <- readr::read_csv("data/whatsapp_status.csv", show_col_types = FALSE) %>%
  dplyr::transmute(
    patient_id = as.character(.data$patient_id),
    patient_name = stringr::str_squish(stringr::str_trim(paste(.data$first_name, .data$last_name)))
  ) %>%
  dplyr::filter(.data$patient_id %in% TARGET_PATIENT_IDS, !is.na(.data$patient_id), .data$patient_name != "") %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(patient_name = dplyr::first(.data$patient_name), .groups = "drop")

excluded_patient_names <- c(
  "elyse tanzer",
  "mary poppins",
  "elliott conecker",
  "old henry jaki"
)
excluded_patient_ids <- whatsapp_names %>%
  dplyr::mutate(patient_name_key = stringr::str_to_lower(stringr::str_squish(.data$patient_name))) %>%
  dplyr::filter(.data$patient_name_key %in% excluded_patient_names) %>%
  dplyr::pull(.data$patient_id) %>%
  unique()

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
    )
  ) %>%
  dplyr::filter(.data$patient_id %in% TARGET_PATIENT_IDS, !is.na(.data$patient_id), !is.na(.data$variant)) %>%
  dplyr::arrange(.data$patient_id, dplyr::desc(.data$status_timestamp_utc)) %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(variant = dplyr::first(.data$variant), .groups = "drop")

events <- all_seizure_events %>%
  dplyr::filter(
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

read_active_months <- function() {
  if (!file.exists("output/tabs/app_usage.json")) {
    return(tibble::tibble())
  }

  usage_records <- jsonlite::fromJSON("output/tabs/app_usage.json", simplifyVector = FALSE)
  purrr::map_dfr(
    usage_records,
    function(rec) {
      patient <- rec$patient_id
      usage <- rec$app_usage_per_month
      if (is.null(patient) || is.null(usage) || length(usage) == 0) {
        return(tibble::tibble())
      }

      usage_vec <- unlist(usage, use.names = TRUE)
      if (is.null(usage_vec) || length(usage_vec) == 0) {
        return(tibble::tibble())
      }

      tibble::tibble(
        patient_id = patient,
        month = as.Date(paste0(names(usage_vec), "-01")),
        app_usage_days = as.numeric(usage_vec)
      )
    }
  ) %>%
    dplyr::mutate(month = as.Date(lubridate::floor_date(.data$month, "month"))) %>%
    dplyr::filter(
      .data$patient_id %in% TARGET_PATIENT_IDS,
      !is.na(.data$patient_id),
      !is.na(.data$month),
      .data$app_usage_days >= 1,
      .data$month <= analysis_month_cap
    ) %>%
    dplyr::distinct(.data$patient_id, .data$month)
}

active_months <- read_active_months()

if (nrow(active_months) == 0) {
  active_months <- events %>%
    dplyr::distinct(.data$patient_id, .data$month)
}

patients_with_seizures <- events %>%
  dplyr::distinct(.data$patient_id) %>%
  dplyr::pull(.data$patient_id)

active_months <- active_months %>%
  dplyr::filter(.data$patient_id %in% patients_with_seizures)

if (nrow(active_months) == 0) {
  warning("No active patient-months available. No outputs were generated.")
  quit(save = "no", status = 0)
}

global_start_month <- as.Date("2025-01-01")

global_end_month <- analysis_month_cap

active_months <- active_months %>%
  dplyr::filter(.data$month >= global_start_month)

event_counts_type <- events %>%
  dplyr::semi_join(active_months, by = c("patient_id", "month")) %>%
  dplyr::count(.data$patient_id, .data$month, .data$seizure_type_plot, name = "seizure_count")

monthly_totals <- active_months %>%
  dplyr::left_join(
    event_counts_type %>%
      dplyr::group_by(.data$patient_id, .data$month) %>%
      dplyr::summarise(total_seizures = sum(.data$seizure_count), .groups = "drop"),
    by = c("patient_id", "month")
  ) %>%
  tidyr::replace_na(list(total_seizures = 0L))

if (length(excluded_patient_ids) > 0) {
  monthly_totals <- monthly_totals %>%
    dplyr::filter(!.data$patient_id %in% excluded_patient_ids)

  event_counts_type <- event_counts_type %>%
    dplyr::filter(!.data$patient_id %in% excluded_patient_ids)
}

included_patients <- monthly_totals %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(total_seizures_all_active_months = sum(.data$total_seizures), .groups = "drop") %>%
  dplyr::filter(.data$total_seizures_all_active_months > 0) %>%
  dplyr::pull(.data$patient_id)

monthly_totals <- monthly_totals %>%
  dplyr::filter(.data$patient_id %in% included_patients)

event_counts_type <- event_counts_type %>%
  dplyr::filter(.data$patient_id %in% included_patients)

if (nrow(monthly_totals) == 0 || nrow(event_counts_type) == 0) {
  warning("No active patients with reported reviewed seizures in the plotting window. No outputs were generated.")
  quit(save = "no", status = 0)
}

patient_order <- monthly_totals %>%
  dplyr::group_by(.data$patient_id) %>%
  dplyr::summarise(total_seizures_all_active_months = sum(.data$total_seizures), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(.data$total_seizures_all_active_months), .data$patient_id) %>%
  dplyr::pull(.data$patient_id)

patient_label_lookup <- tibble::tibble(patient_id = patient_order) %>%
  dplyr::left_join(whatsapp_names, by = "patient_id") %>%
  dplyr::left_join(variants, by = "patient_id") %>%
  dplyr::mutate(
    patient_name = dplyr::coalesce(.data$patient_name, .data$patient_id),
    variant = dplyr::coalesce(.data$variant, "Unknown"),
    variant = dplyr::recode(
      .data$variant,
      "K1473K + Pro1428_Lys1473del [predicted inframe exon skipping]" =
        "c.4419+1A>G"
    ),
    patient_label = paste0(.data$patient_name, " (", .data$variant, ")")
  ) %>%
  dplyr::select(patient_id, patient_label)

patient_label_map <- stats::setNames(patient_label_lookup$patient_label, patient_label_lookup$patient_id)

type_order <- event_counts_type %>%
  dplyr::group_by(.data$seizure_type_plot) %>%
  dplyr::summarise(total_events = sum(.data$seizure_count), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(.data$total_events), .data$seizure_type_plot) %>%
  dplyr::pull(.data$seizure_type_plot)

monthly_totals <- monthly_totals %>%
  dplyr::mutate(patient_id = factor(.data$patient_id, levels = patient_order))

plot_counts <- event_counts_type %>%
  dplyr::mutate(
    patient_id = factor(.data$patient_id, levels = patient_order),
    seizure_type_plot = factor(.data$seizure_type_plot, levels = type_order)
  )

missing_types <- setdiff(type_order, names(SEIZURE_TYPE_COLORS))
if (length(missing_types) > 0) {
  extra_pal <- setNames(scales::hue_pal()(length(missing_types)), missing_types)
  type_colors <- c(SEIZURE_TYPE_COLORS, extra_pal)
} else {
  type_colors <- SEIZURE_TYPE_COLORS
}
type_colors <- type_colors[intersect(names(type_colors), type_order)]
legend_seed_data <- tibble::tibble(
  month = global_start_month,
  seizure_count = 0,
  seizure_type_plot = factor(type_order, levels = type_order)
)

n_patients <- length(patient_order)
ncol_facets <- if (n_patients <= 16) 4 else if (n_patients <= 30) 5 else 6
nrow_facets <- ceiling(n_patients / ncol_facets)

patient_row_map <- tibble::tibble(
  patient_id = factor(patient_order, levels = patient_order),
  panel_index = seq_along(patient_order),
  facet_row = ((.data$panel_index - 1L) %/% ncol_facets) + 1L
) %>%
  dplyr::select(patient_id, facet_row)

monthly_totals <- monthly_totals %>%
  dplyr::left_join(patient_row_map, by = "patient_id")

plot_counts <- plot_counts %>%
  dplyr::left_join(patient_row_map, by = "patient_id")

row_scale_reference <- monthly_totals %>%
  dplyr::group_by(.data$facet_row) %>%
  dplyr::summarise(row_max_monthly_seizures = max(.data$total_seizures, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    row_max_monthly_seizures = dplyr::if_else(.data$row_max_monthly_seizures <= 0, 1L, .data$row_max_monthly_seizures)
  )

build_row_plot <- function(row_id) {
  row_patients <- patient_row_map %>%
    dplyr::filter(.data$facet_row == row_id) %>%
    dplyr::pull(.data$patient_id) %>%
    as.character()

  row_totals <- monthly_totals %>%
    dplyr::filter(.data$facet_row == row_id) %>%
    dplyr::mutate(patient_id = factor(as.character(.data$patient_id), levels = row_patients))

  row_counts <- plot_counts %>%
    dplyr::filter(.data$facet_row == row_id) %>%
    dplyr::mutate(patient_id = factor(as.character(.data$patient_id), levels = row_patients))

  row_max <- row_scale_reference %>%
    dplyr::filter(.data$facet_row == row_id) %>%
    dplyr::pull(.data$row_max_monthly_seizures) %>%
    first()

  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = legend_seed_data,
      ggplot2::aes(x = .data$month, y = .data$seizure_count, fill = .data$seizure_type_plot),
      width = 25,
      alpha = 0,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +
    ggplot2::geom_col(
      data = row_totals,
      ggplot2::aes(x = .data$month, y = .data$total_seizures),
      fill = "#ECECEC",
      width = 25
    ) +
    ggplot2::geom_col(
      data = row_counts,
      ggplot2::aes(x = .data$month, y = .data$seizure_count, fill = .data$seizure_type_plot),
      width = 25
    ) +
    ggplot2::facet_wrap(
      ~ patient_id,
      ncol = ncol_facets,
      scales = "fixed",
      drop = FALSE,
      labeller = ggplot2::as_labeller(patient_label_map)
    ) +
    ggplot2::scale_x_date(
      date_labels = "%m/%y",
      date_breaks = "1 month",
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, row_max),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    ggplot2::scale_fill_manual(
      values = type_colors,
      breaks = type_order,
      limits = type_order,
      drop = FALSE
    ) +
    ggplot2::labs(
      x = "Month",
      y = "Seizure count",
      fill = "Seizure type"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "#F5F5F5", color = "#DDDDDD"),
      strip.text = ggplot2::element_text(size = 8),
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom",
      legend.key.height = ggplot2::unit(0.35, "cm"),
      legend.key.width = ggplot2::unit(0.9, "cm")
    ) +
    ggplot2::coord_cartesian(xlim = c(global_start_month, global_end_month))
}

row_ids <- sort(unique(patient_row_map$facet_row))
row_plots <- purrr::map2(
  row_ids,
  seq_along(row_ids),
  function(row_id, idx) {
    row_plot <- build_row_plot(row_id)
    is_bottom_row <- idx == length(row_ids)
    if (!is_bottom_row) {
      row_plot <- row_plot +
        ggplot2::theme(
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          legend.position = "none"
        )
    }
    row_plot
  }
)

p <- (
  patchwork::wrap_plots(
    row_plots,
    ncol = 1
  ) +
    patchwork::plot_annotation(
      title = "Patient-Month Seizure Frequency",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold")
      )
    )
)

plot_width <- max(14, 2.6 * ncol_facets)
plot_height <- max(9, 2.0 * nrow_facets + 1.5)

ggplot2::ggsave(
  filename = file.path(OUTPUT_FIG_DIR, "seizure_patient_month.pdf"),
  plot = p,
  width = plot_width,
  height = plot_height
)

readr::write_csv(
  plot_counts %>%
    dplyr::mutate(
      patient_id = as.character(.data$patient_id),
      seizure_type_plot = as.character(.data$seizure_type_plot)
    ) %>%
    dplyr::arrange(.data$patient_id, .data$month, .data$seizure_type_plot),
  file.path(OUTPUT_TAB_DIR, "seizure_patient_month_counts_by_type.csv")
)

readr::write_csv(
  monthly_totals %>%
    dplyr::mutate(patient_id = as.character(.data$patient_id)) %>%
    dplyr::arrange(.data$patient_id, .data$month),
  file.path(OUTPUT_TAB_DIR, "seizure_patient_month_total_counts.csv")
)

readr::write_csv(
  row_scale_reference,
  file.path(OUTPUT_TAB_DIR, "seizure_patient_month_row_scale_reference.csv")
)
