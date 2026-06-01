# Seizure type proportions over time

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

source("src/desc/seizure_type_standardization.R")

OUTPUT_FIG_DIR <- "output/figs/seizure_patterns"
dir.create(OUTPUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
png_output_path <- file.path(OUTPUT_FIG_DIR, "sz_props.png")
if (file.exists(png_output_path)) {
  invisible(file.remove(png_output_path))
}

SEIZURE_TYPE_COLORS <- c(
  `Tonic-clonic` = "#5698a3",  
  `Focal`        = "#ffde76",  
  `Tonic`        = "#67771a",  
  `Myoclonic`    = "#0076c0",  
  `Absence`      = "#e37c1d",  
  `Spasms`       = "#7a5072",  
  `Other`        = "#B0B0B0"  
)

# Read events
events <- readr::read_csv("data/events.csv", show_col_types = FALSE) %>%
  standardize_seizure_events(filter_to_seizure = TRUE) %>%
  dplyr::filter(
    .data$patient_id %in% TARGET_PATIENT_IDS,
    !is.na(.data$month),
    !is.na(.data$patient_id)
  )

# Count each patient once per seizure group per month
monthly_counts <- events %>%
  dplyr::distinct(.data$patient_id, .data$month, .data$seizure_group) %>%
  dplyr::count(.data$month, .data$seizure_group, name = "n")

# Order seizure groups
group_levels <- c("Tonic-clonic", "Focal", "Other", "Tonic", "Myoclonic", "Absence", "Spasms")
monthly_counts <- monthly_counts %>%
  dplyr::mutate(seizure_group = factor(.data$seizure_group, levels = group_levels)) %>%
  dplyr::arrange(.data$month, .data$seizure_group)

# Compute total number of patients per month
patients_total <- monthly_counts %>%
  group_by(month) %>%
  summarise(total = sum(n))

# Compute proportion of each seizure type per month
props_df <- left_join(monthly_counts, patients_total, by = "month") %>%
  mutate(prop = n / total) %>%
  filter(month >= "2025-01-01" & month <= "2026-03-01")

# Plot
p_sz_props <- ggplot(props_df, aes(x = month, y = prop, fill = seizure_group)) +
  geom_area(stat = "smooth", method = "loess", position = "identity") +
  geom_line(stat = "smooth", method = "loess", formula = y ~ x, se = FALSE, linetype = 1,
            aes(color = seizure_group), show.legend = FALSE) +
  labs(x = "Month/Year", y = "Proportion", fill = "Seizure Type") +
  coord_cartesian(ylim = c(0,0.6)) +
  scale_x_date(date_labels = "%m/%Y", date_breaks = "1 month", expand = expansion(mult = c(0.01, 0.01))) +
  theme_classic() +
  theme(legend.justification = c(0.05, 1), legend.position = c(0.05, 1)) +
  scale_fill_manual(values = SEIZURE_TYPE_COLORS) +
  scale_color_manual(values = SEIZURE_TYPE_COLORS)

ggsave(
  filename = file.path(OUTPUT_FIG_DIR, "sz_props.pdf"),
  plot = p_sz_props,
  width = 10,
  height = 7
)
