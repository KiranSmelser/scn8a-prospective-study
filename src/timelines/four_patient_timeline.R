# Create a four-panel timeline figure following the layout used for Figure 16
# in the Citizen SCN2A study.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# Panels are arranged row-wise in the order requested:
# Petak | Bejarano
# Kowalak | Jaki
four_panel_patients <- c(
  "5dcb3f4d559dc95b43af4872", # Oktavija Petak
  "68d58414b783c7baf7292504", # Peto Bejarano
  "67ead18b773a9ee7d1911a39", # Michal Kowalak
  "6954f548ce6a468b362fce7b"  # Henry Jaki
)

# Loading the timeline analysis with a patient override avoids regenerating
# every individual patient PDF while retaining the canonical plotting logic.
options(scn8a.timeline_patient_filter = four_panel_patients)
source(file.path("src", "timelines", "timelines.R"))
options(scn8a.timeline_patient_filter = NULL)

missing_patients <- setdiff(four_panel_patients, patients_with_data)
if (length(missing_patients) > 0) {
  stop(
    "Four-panel patient IDs missing from timeline data: ",
    paste(missing_patients, collapse = ", ")
  )
}

make_four_panel_plot <- function(patient_id) {
  panel_plot <- plot_patient_timeline(patient_id)
  if (is.null(panel_plot)) {
    stop("No timeline could be generated for patient: ", patient_id)
  }

  # Match the reduced-size marks and typography used in the Citizen SCN2A
  # four-panel figure, while preserving the SCN8A plot's encodings.
  for (layer_index in seq_along(panel_plot$layers)) {
    if (!is.null(panel_plot$layers[[layer_index]]$aes_params$size)) {
      panel_plot$layers[[layer_index]]$aes_params$size <-
        panel_plot$layers[[layer_index]]$aes_params$size * 0.48
    }
    if (!is.null(panel_plot$layers[[layer_index]]$aes_params$linewidth)) {
      panel_plot$layers[[layer_index]]$aes_params$linewidth <-
        panel_plot$layers[[layer_index]]$aes_params$linewidth * 0.48
    }
  }

  panel_plot +
    guides(color = "none", fill = "none", shape = "none") +
    theme(
      plot.title = element_text(size = 7, hjust = 0, face = "bold"),
      plot.subtitle = element_text(size = 5.5),
      axis.title.x = element_text(size = 6.5),
      axis.text.x = element_text(size = 5, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 5.5, color = "#222222"),
      axis.ticks = element_line(linewidth = 0.25),
      panel.grid.major = element_line(color = "#A0A0A0", linewidth = 0.2),
      panel.border = element_rect(linewidth = 0.35),
      plot.margin = margin(12, 18, 12, 18)
    )
}

four_panel_plots <- lapply(four_panel_patients, make_four_panel_plot)

four_panel_patchwork <- wrap_plots(four_panel_plots, ncol = 2, byrow = TRUE) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(26, 32, 26, 32)
    )
  ) &
  theme(legend.position = "none")

dir.create(file.path("output", "figs"), recursive = TRUE, showWarnings = FALSE)
figure_pdf <- file.path("output", "figs", "four_patient_timelines.pdf")

ggsave(
  filename = figure_pdf,
  plot = four_panel_patchwork,
  width = 13.33,
  height = 7.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

message(
  "Four-patient timeline figure written to: ",
  normalizePath(figure_pdf)
)
