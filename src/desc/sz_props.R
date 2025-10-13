# Seizure type proportions over time (area plot)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
})

SEIZURE_TYPE_COLORS <- c(
  `Tonic-clonic` = "#71D0F5FF",  
  `Focal`        = "#FD8CC1FF",  
  `Tonic`        = "#FED439FF",  
  `Myoclonic`    = "#FD7446FF",  
  `Absence`      = "#C80813FF",  
  `Spasms`       = "#197EC0FF",  
  `Other`        = "#B0B0B0"  
)

# Map raw seizure_type strings
map_seizure_type <- function(x) {
  s <- tolower(trimws(x))
  s <- gsub("_", " ", s)
  s <- gsub("[^a-z ]", " ", s)
  s <- stringr::str_squish(s)

  dplyr::case_when(
    stringr::str_detect(s, "spasm") ~ "Spasms",
    stringr::str_detect(s, "absence") ~ "Absence",
    stringr::str_detect(s, "myoc") ~ "Myoclonic",
    stringr::str_detect(s, "tonic ?clonic|tonicclonic|grand mal|grote aanval|gtcs") ~ "Tonic-clonic",
    stringr::str_detect(s, "\\bfocal\\b|impaired awareness|aware|hyperkinetic|gelastic|dacrystic") ~ "Focal",
    stringr::str_detect(s, "\\btonic\\b") ~ "Tonic",
    stringr::str_detect(s, "clonic") ~ "Tonic",
    TRUE ~ "Other"
  )
}

# Read events
events <- readr::read_csv("data/events.csv", show_col_types = FALSE) %>%
  dplyr::filter(tolower(.data$type) == "seizure") %>%
  dplyr::mutate(
    date_time = lubridate::ymd_hms(.data$date, quiet = TRUE, tz = "UTC"),
    month = as.Date(lubridate::floor_date(.data$date_time, "month")),
    seizure_type_clean = dplyr::if_else(is.na(.data$seizure_type) | .data$seizure_type == "",
                                        "Other", .data$seizure_type),
    seizure_group = map_seizure_type(.data$seizure_type_clean)
  ) %>%
  dplyr::filter(!is.na(.data$month), !is.na(.data$patient_id))

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
  filter(month >= "2025-01-01" & month <= "2025-10-01")

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

ggsave(filename = "output/figs/sz_props.png", plot = p_sz_props, width = 10, height = 7, dpi = 600)