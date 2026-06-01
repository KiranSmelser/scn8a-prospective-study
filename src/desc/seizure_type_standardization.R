# Helpers to standardize patient-reported seizure types.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

source("src/analysis_config.R")

normalize_seizure_text_for_match <- function(x) {
  normalized <- x %>%
    as.character() %>%
    tidyr::replace_na("") %>%
    stringr::str_squish()

  transliterated <- iconv(normalized, from = "", to = "ASCII//TRANSLIT")
  transliterated[is.na(transliterated)] <- normalized[is.na(transliterated)]

  transliterated %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("['`]", "") %>%
    stringr::str_replace_all("[^a-z0-9]+", " ") %>%
    stringr::str_squish()
}

get_manual_seizure_overrides <- function() {
  c(
    "generalized motor seizure" = "Tonic-clonic",
    "uogolniony napad toniczny z krzykiem" = "Tonic",
    "uog lniony napad toniczny z krzykiem" = "Tonic",
    "seizure where she became stiff and not breathing well tonic" = "Tonic",
    "im assuming its generalized absence seizure" = "Absence",
    "i a m assuming its generalized absence seizure" = "Absence",
    "focal impaired awareness myoclonic seizure" = "Focal",
    "focal to bilateral clonic seizure" = "Focal",
    "generalized tonic seizure" = "Tonic",
    "had 2 more of these" = "Focal",
    "had 2 of these" = "Focal",
    "had at least 4 of these followed by vomiting" = "Focal",
    "he had at least 10 focal seizures today followed by vomiting both ears are infected again" = "Focal",
    "absent with vomiting following not sure if theyvare galastic" = "Absence",
    "looks like dystonic posturing" = "Focal",
    "multiple 30 second focals his trileptal levels came back low i think his antibiotic took his levels down" = "Focal"
  )
}

lookup_manual_seizure_override <- function(match_text) {
  overrides <- get_manual_seizure_overrides()
  value <- unname(overrides[match_text])
  if (length(value) == 0 || is.na(value) || value == "") {
    if (!is.na(match_text) && stringr::str_detect(match_text, "uog.*napad toniczny z krzykiem")) {
      return("Tonic")
    }
    return(NA_character_)
  }
  as.character(value[[1]])
}

extract_standardized_seizure_components <- function(raw_type) {
  text <- normalize_seizure_text_for_match(raw_type)

  if (is.na(text) || text == "") {
    return(character())
  }

  types <- character()

  add_type <- function(condition, value) {
    if (isTRUE(condition)) {
      types <<- c(types, value)
    }
  }

  is_device_detected <- stringr::str_detect(text, "nightwatch|detected event|gedecteerde")
  is_status <- stringr::str_detect(text, "\\bstatus\\b|subcliniczl")
  is_spasms <- stringr::str_detect(text, "\\bspasm")
  is_absence <- stringr::str_detect(text, "\\babsence\\b|\\babsent\\b|petite crise")
  is_myoclonic <- stringr::str_detect(text, "\\bmyoc")
  is_tonic_clonic <- stringr::str_detect(
    text,
    "focal to bilateral (tonic )?clonic|tonic ?clonic|tonicclonic|grand mal|\\bgtc\\b|grote aanval|tonic clinic"
  )
  is_generalized_motor <- stringr::str_detect(
    text,
    "uogolniony napad toniczny|generalized motor seizure|generalized non t c|generalized tonic seizure|\\bgeneralized motor\\b"
  )
  is_focal <- stringr::str_detect(
    text,
    "\\bfocal\\b|\\bfocals\\b|gelastic|dacrystic|hyperkinetic|focal motora|focal arm|\\bfocal r\\b|\\bfocal l\\b|occipital onset"
  )
  is_tonic <- stringr::str_detect(text, "\\btonic\\b") && !is_tonic_clonic && !is_generalized_motor
  is_clonic <- stringr::str_detect(text, "\\bclonic\\b") && !is_tonic_clonic
  is_atonic <- stringr::str_detect(text, "\\batonic\\b|\\bdrop\\b")

  add_type(is_device_detected, "Device-detected event")
  add_type(is_status, "Status epilepticus")
  add_type(is_spasms, "Spasms")
  add_type(is_absence, "Absence")
  add_type(is_myoclonic, "Myoclonic")
  add_type(is_tonic_clonic, "Tonic-clonic")
  add_type(is_generalized_motor, "Generalized motor (non tonic-clonic)")
  add_type(is_focal, "Focal")
  add_type(is_tonic, "Tonic")
  add_type(is_clonic, "Clonic")
  add_type(is_atonic, "Atonic")

  sort(unique(types))
}

select_primary_seizure_type <- function(components, raw_type) {
  if (length(components) == 0) {
    raw_text <- stringr::str_squish(as.character(raw_type))
    if (is.na(raw_text) || raw_text == "") {
      return("Unspecified")
    }
    return("Unknown/Unmapped")
  }

  if (length(components) > 1) {
    return("Mixed/Multiple")
  }

  components[[1]]
}

map_components_to_analysis_group <- function(components) {
  if (length(components) == 0) {
    return("Other")
  }

  if ("Spasms" %in% components) {
    return("Spasms")
  }
  if ("Absence" %in% components) {
    return("Absence")
  }
  if ("Myoclonic" %in% components) {
    return("Myoclonic")
  }
  if ("Tonic-clonic" %in% components) {
    return("Tonic-clonic")
  }
  if ("Focal" %in% components) {
    return("Focal")
  }
  if (any(c("Tonic", "Clonic") %in% components)) {
    return("Tonic")
  }

  "Other"
}

build_seizure_review_reasons <- function(raw_type, match_text, mapping_count, components, manual_override_type) {
  if (!is.na(manual_override_type) && manual_override_type != "") {
    return(NA_character_)
  }

  reasons <- character()

  raw_text <- stringr::str_squish(as.character(raw_type))
  raw_text_lower <- stringr::str_to_lower(raw_text)
  raw_is_blank <- is.na(raw_text) || raw_text == ""

  if (mapping_count == 0 && !raw_is_blank) {
    reasons <- c(reasons, "unmapped_free_text")
  }
  if (mapping_count > 1) {
    reasons <- c(reasons, "multiple_seizure_types_mentioned")
  }
  if (stringr::str_detect(raw_text_lower, "\\?|unknown|not sure|difficult to say|assuming|looks like")) {
    reasons <- c(reasons, "uncertain_language")
  }
  if (stringr::str_detect(match_text, "nightwatch|detected event|gedecteerde")) {
    reasons <- c(reasons, "device_detected_event_needs_confirmation")
  }
  if (stringr::str_detect(match_text, "vomiting|not breathing|infected|antibiotic|scotoma|dystonic|eyes face")) {
    reasons <- c(reasons, "symptom_or_context_text")
  }
  if ("Generalized motor (non tonic-clonic)" %in% components) {
    reasons <- c(reasons, "generalized_non_specific_type")
  }

  reasons <- unique(reasons)
  if (length(reasons) == 0) {
    return(NA_character_)
  }

  paste(reasons, collapse = "; ")
}

standardize_seizure_events <- function(events_df, filter_to_seizure = TRUE) {
  events <- events_df

  if (!"type" %in% names(events)) {
    events$type <- NA_character_
  }
  if (!"seizure_type" %in% names(events)) {
    events$seizure_type <- NA_character_
  }
  if (!"date" %in% names(events)) {
    events$date <- NA_character_
  }

  events <- events %>%
    dplyr::mutate(
      type_clean = stringr::str_to_lower(stringr::str_squish(as.character(.data$type)))
    )

  if (isTRUE(filter_to_seizure)) {
    events <- events %>%
      dplyr::filter(.data$type_clean == "seizure")
  }

  events %>%
    dplyr::mutate(
      seizure_type_raw = stringr::str_squish(as.character(.data$seizure_type)),
      seizure_type_raw = dplyr::if_else(is.na(.data$seizure_type_raw), "", .data$seizure_type_raw),
      seizure_type_match_text = normalize_seizure_text_for_match(.data$seizure_type_raw),
      manual_override_type = purrr::map_chr(.data$seizure_type_match_text, lookup_manual_seizure_override),
      seizure_components = purrr::map2(
        .data$seizure_type_raw,
        .data$manual_override_type,
        function(raw_type, override_type) {
          if (!is.na(override_type) && override_type != "") {
            return(override_type)
          }
          extract_standardized_seizure_components(raw_type)
        }
      ),
      mapped_component_count = lengths(.data$seizure_components),
      seizure_type_standardized = purrr::map_chr(
        .data$seizure_components,
        ~ if (length(.x) == 0) NA_character_ else paste(.x, collapse = " + ")
      ),
      seizure_type_standardized = dplyr::coalesce(
        .data$seizure_type_standardized,
        dplyr::if_else(.data$seizure_type_raw == "", "Unspecified", "Unknown/Unmapped")
      ),
      seizure_type_primary = purrr::map2_chr(
        .data$seizure_components,
        .data$seizure_type_raw,
        select_primary_seizure_type
      ),
      seizure_group = purrr::map_chr(.data$seizure_components, map_components_to_analysis_group),
      review_reason = purrr::pmap_chr(
        list(
          .data$seizure_type_raw,
          .data$seizure_type_match_text,
          .data$mapped_component_count,
          .data$seizure_components,
          .data$manual_override_type
        ),
        build_seizure_review_reasons
      ),
      requires_review = !is.na(.data$review_reason),
      standardization_status = dplyr::case_when(
        !is.na(.data$manual_override_type) ~ "mapped_manual",
        .data$mapped_component_count == 0 ~ "unmapped",
        .data$mapped_component_count == 1 ~ "mapped_single",
        TRUE ~ "mapped_multiple"
      ),
      date_time = suppressWarnings(lubridate::ymd_hms(.data$date, quiet = TRUE, tz = "UTC")),
      event_date = as.Date(.data$date_time),
      month = as.Date(lubridate::floor_date(.data$date_time, "month"))
    ) %>%
    dplyr::filter(is.na(.data$event_date) | .data$event_date <= ANALYSIS_CUTOFF_DATE) %>%
    dplyr::select(-seizure_components)
}
