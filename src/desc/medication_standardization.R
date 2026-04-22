# Helpers to standardize non-rescue epilepsy medication names.

suppressPackageStartupMessages({
  library(tidyverse)
})

normalize_med_boolean_flag <- function(x, default = FALSE) {
  if (is.logical(x)) {
    return(tidyr::replace_na(x, default))
  }

  normalized <- stringr::str_to_lower(stringr::str_trim(as.character(x)))
  parsed <- dplyr::case_when(
    normalized %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    normalized %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
  tidyr::replace_na(parsed, default)
}

normalize_med_text_for_match <- function(x) {
  normalized <- x %>%
    as.character() %>%
    tidyr::replace_na("") %>%
    stringr::str_squish()

  transliterated <- iconv(normalized, from = "", to = "ASCII//TRANSLIT")
  transliterated[is.na(transliterated)] <- normalized[is.na(transliterated)]

  transliterated %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("&", " and ") %>%
    stringr::str_replace_all("[^a-z0-9]+", " ") %>%
    stringr::str_squish()
}

extract_standardized_med_components <- function(raw_name) {
  text <- normalize_med_text_for_match(raw_name)

  if (is.na(text) || text == "") {
    return(character())
  }

  meds <- character()

  add_med <- function(condition, value) {
    if (isTRUE(condition)) {
      meds <<- c(meds, value)
    }
  }

  add_med(stringr::str_detect(text, "\\b(cenobamate|xcopri|cenobamato)\\b"), "Cenobamate")
  add_med(stringr::str_detect(text, "\\b(brivaracetam|briviact)\\b"), "Brivaracetam")
  add_med(stringr::str_detect(text, "\\b(cannabidiol|epidiolex|epidyolex|cbd)\\b"), "Cannabidiol")
  add_med(stringr::str_detect(text, "\\b(carbamazepine|carmazepine|tegretol|timonil|cbz)\\b"), "Carbamazepine")
  add_med(stringr::str_detect(text, "\\b(clobazam|onfi|sympazan|clb)\\b"), "Clobazam")
  add_med(stringr::str_detect(text, "\\b(clonazepam|klonopin|czp)\\b"), "Clonazepam")
  add_med(stringr::str_detect(text, "\\b(eslicarbazepine|aptiom)\\b"), "Eslicarbazepine Acetate")
  add_med(stringr::str_detect(text, "\\b(ethosuximide|zarontin)\\b"), "Ethosuximide")
  add_med(stringr::str_detect(text, "\\b(felbamate|felbatol)\\b"), "Felbamate")
  add_med(stringr::str_detect(text, "\\bgabapentin\\b"), "Gabapentin")
  add_med(stringr::str_detect(text, "\\b(lacosamide|lacasomide|lacosamida|vimpat)\\b"), "Lacosamide")
  add_med(stringr::str_detect(text, "\\b(lamotrigine|lamotrigina|lamictal|neural)\\b"), "Lamotrigine")
  add_med(stringr::str_detect(text, "\\b(levetiracetam|keppra|lev)\\b"), "Levetiracetam")
  add_med(stringr::str_detect(text, "\\b(oxcarbazepine|trileptal|oxtellar|oxc)\\b"), "Oxcarbazepine")
  add_med(stringr::str_detect(text, "\\b(phenobarbital|phb)\\b"), "Phenobarbital")
  add_med(stringr::str_detect(text, "\\b(phenytoin|fenitoina|hidantal|pht)\\b"), "Phenytoin")
  add_med(stringr::str_detect(text, "\\b(relutrigine|praxis)\\b"), "Praxis")
  add_med(stringr::str_detect(text, "\\b(rufinamide|banzel)\\b"), "Rufinamide")
  add_med(stringr::str_detect(text, "\\b(topiramate|topamax|topiramato)\\b"), "Topiramate")
  add_med(stringr::str_detect(text, "\\b(valproic acid|valproate sodium|sodium valproate|divalproex|depakene|depakote|depacon|vpa|valporate)\\b"), "Valproic Acid")
  add_med(stringr::str_detect(text, "\\b(vigabatrin|sabril)\\b"), "Vigabatrin")
  add_med(stringr::str_detect(text, "\\b(zonisamide|zonegran)\\b"), "Zonisamide")

  sort(unique(meds))
}

build_med_review_reasons <- function(match_text, mapping_count, treatment_type_clean) {
  reasons <- character()

  if (mapping_count == 0) {
    reasons <- c(reasons, "unmapped_name")
  }
  if (mapping_count > 1) {
    reasons <- c(reasons, "multiple_medications_in_name")
  }
  if (!is.na(treatment_type_clean) && treatment_type_clean != "drug") {
    reasons <- c(reasons, "non_drug_treatment_type")
  }
  if (stringr::str_detect(match_text, "\\bdemo\\b")) {
    reasons <- c(reasons, "demo_entry")
  }
  if (stringr::str_detect(match_text, "\\b(trial|praxis)\\b|open\\s+label") &&
      !stringr::str_detect(match_text, "\\b(relutrigine|praxis)\\b")) {
    reasons <- c(reasons, "trial_or_placeholder_name")
  }
  if (stringr::str_detect(match_text, "\\b(keto|ketogenic|cetogenica|diet|exercise|exercices)\\b")) {
    reasons <- c(reasons, "diet_or_exercise_entry")
  }
  if (stringr::str_detect(match_text, "painpatrol|menthol|melatonin|sodium chloride|sodium bicarbonate|amoxicillin|vaxigrip|pepcid|levocarnitine|oracit|potassium citrate|promethazine|methylphenidate") &&
      !stringr::str_detect(match_text, "\\b(cannabidiol|epidiolex|epidyolex|cbd)\\b")) {
    reasons <- c(reasons, "possible_non_asm_or_support_med")
  }

  reasons <- unique(reasons)
  if (length(reasons) == 0) {
    return(NA_character_)
  }
  paste(reasons, collapse = "; ")
}

standardize_non_rescue_epilepsy_medications <- function(medications_df, include_non_drug = FALSE) {
  medications <- medications_df

  if (!"is_deleted" %in% names(medications)) {
    medications$is_deleted <- FALSE
  }
  if (!"is_rescue_med" %in% names(medications)) {
    medications$is_rescue_med <- FALSE
  }
  if (!"reason" %in% names(medications)) {
    medications$reason <- NA_character_
  }
  if (!"treatment_type" %in% names(medications)) {
    medications$treatment_type <- NA_character_
  }
  if (!"name" %in% names(medications)) {
    medications$name <- NA_character_
  }

  medications <- medications %>%
    dplyr::mutate(
      is_deleted = normalize_med_boolean_flag(.data$is_deleted),
      is_rescue_med = normalize_med_boolean_flag(.data$is_rescue_med),
      reason_clean = stringr::str_to_lower(stringr::str_squish(as.character(.data$reason))),
      treatment_type_clean = stringr::str_to_lower(stringr::str_squish(as.character(.data$treatment_type))),
      name = stringr::str_squish(as.character(.data$name)),
      match_text = normalize_med_text_for_match(.data$name),
      med_components = purrr::map(.data$name, extract_standardized_med_components),
      mapped_component_count = lengths(.data$med_components),
      standardized_components = purrr::map_chr(
        .data$med_components,
        ~ if (length(.x) == 0) NA_character_ else paste(.x, collapse = " + ")
      ),
      name_standardized = dplyr::coalesce(
        .data$standardized_components,
        dplyr::if_else(is.na(.data$name) | .data$name == "", "Unknown medication", .data$name)
      ),
      review_reason = purrr::pmap_chr(
        list(.data$match_text, .data$mapped_component_count, .data$treatment_type_clean),
        build_med_review_reasons
      ),
      requires_review = !is.na(.data$review_reason),
      standardization_status = dplyr::case_when(
        .data$mapped_component_count == 0 ~ "unmapped",
        .data$mapped_component_count == 1 ~ "mapped_single",
        TRUE ~ "mapped_multiple"
      )
    ) %>%
    dplyr::filter(
      !.data$is_deleted,
      !.data$is_rescue_med,
      .data$reason_clean == "epilepsy"
    )

  if (!include_non_drug) {
    medications <- medications %>%
      dplyr::filter(.data$treatment_type_clean == "drug")
  }

  medications %>%
    dplyr::select(-med_components, -match_text)
}
