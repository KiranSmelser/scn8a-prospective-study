# Shared analysis date limits.

ANALYSIS_CUTOFF_DATE <- as.Date("2026-05-31")

TARGET_PATIENT_IDS <- c(
  "67b6126de8efcd6464eeb8c4",
  "67b59b49e8efcd6464eeb595",
  "690a3a45b2c2dee67c5423d6",
  "68d58414b783c7baf7292504",
  "691b83cb7db39ba79465d9cf",
  "6808969d5413b5941d3d58d8",
  "67fd6a03773a9ee7d191cbea",
  "67e82887773a9ee7d190fb0a",
  "5dcb3f4d559dc95b43af4872",
  "67ead18b773a9ee7d1911a39",
  "67d181092099f5b1ffcbc0db",
  "67b63d09e8efcd6464eebcba",
  "6836626b1191a9613d19ac21",
  "67d0c5202099f5b1ffcbbbc9",
  "67a8e3f7890c36004c181796",
  "67b529afe8efcd6464eeb287",
  "67e5f605773a9ee7d190e931",
  "6954f548ce6a468b362fce7b",
  "67b60d2ee8efcd6464eeb882",
  "688de8e1e20e57ca0b5335eb",
  "681876911191a9613d188b77",
  "67b562cae8efcd6464eeb345",
  "693e0b36ce6a468b362ee35d"
)

analysis_end_date <- function(reference_date = Sys.Date()) {
  min(as.Date(reference_date), ANALYSIS_CUTOFF_DATE)
}

last_complete_analysis_month_end <- function(reference_date = Sys.Date()) {
  reference_month_start <- as.Date(format(as.Date(reference_date), "%Y-%m-01"))
  previous_month_end <- reference_month_start - 1
  min(previous_month_end, ANALYSIS_CUTOFF_DATE)
}
