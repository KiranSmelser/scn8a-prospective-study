# Run the full SCN8A prospective-study analysis pipeline.

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
trailing_args <- commandArgs(trailingOnly = TRUE)
continue_on_error <- "--continue-on-error" %in% trailing_args ||
  identical(Sys.getenv("RUN_ALL_CONTINUE_ON_ERROR"), "true")
dry_run <- "--dry-run" %in% trailing_args

script_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(script_arg) > 0) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("src/run_all.R", mustWork = FALSE)
}

find_project_root <- function(start_dir) {
  current <- normalizePath(start_dir, mustWork = TRUE)

  repeat {
    has_project_markers <- dir.exists(file.path(current, "src")) &&
      dir.exists(file.path(current, "data")) &&
      length(list.files(current, pattern = "\\.Rproj$", full.names = TRUE)) > 0

    if (has_project_markers) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find project root from: ", start_dir, call. = FALSE)
    }

    current <- parent
  }
}

project_root <- find_project_root(dirname(script_path))
setwd(project_root)

required_inputs <- c(
  "data/events.csv",
  "data/forms.csv",
  "data/form_answers.csv",
  "data/med_dosages.csv",
  "data/med_intakes.csv",
  "data/medications.csv",
  "data/patient_summary_metrics.json",
  "data/patients.csv",
  "data/prospective_development_milestones.csv",
  "data/prospective_surveys.csv",
  "data/registry.csv",
  "data/whatsapp_status.csv"
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop(
    "Required input file(s) are missing from data/:\n  ",
    paste(missing_inputs, collapse = "\n  "),
    call. = FALSE
  )
}

pipeline_steps <- data.frame(
  phase = c(
    rep("shared helpers", 4),
    rep("raw-data summaries", 14),
    rep("cohort build", 2),
    "patient-month panel",
    "cohort seizure summary",
    rep("clustering", 5),
    rep("prediction", 5),
    "development",
    rep("modeling", 2),
    rep("changepoints", 3),
    rep("periodicity", 4),
    rep("episodicity", 2),
    rep("timelines", 2)
  ),
  script = c(
    "src/analysis_config.R",
    "src/desc/medication_standardization.R",
    "src/desc/seizure_type_standardization.R",
    "src/data_corrections.R",
    "src/desc/app_usage.R",
    "src/desc/module_overlap_upset.R",
    "src/desc/weekly_survey.R",
    "src/desc/medication_standardization_report.R",
    "src/desc/medication_patterns.R",
    "src/desc/med_avgs.R",
    "src/desc/seizure_type_standardization_report.R",
    "src/desc/seizures_by_day.R",
    "src/desc/sz_avgs.R",
    "src/desc/sz_props.R",
    "src/desc/seizure_type_combinations.R",
    "src/desc/seizure_patient_month.R",
    "src/desc/seizure_burden_by_type.R",
    "src/desc/variant_vs_seizure_type.R",
    "src/seizures/build_cohort_seizure_dataset.R",
    "src/seizures/check_poisson_fit.R",
    "src/modeling/patient_month_panel.R",
    "src/desc/cohort_seizure_summary.R",
    "src/clustering/seizure_freq_features.R",
    "src/clustering/clustering.R",
    "src/clustering/summarize_clusters.R",
    "src/clustering/clusters_alluvial.R",
    "src/clustering/plot_clustering.R",
    "src/modeling/prediction/cluster_pred_sz.R",
    "src/modeling/prediction/cluster_pred_registry.R",
    "src/modeling/prediction/cluster_pred_combined.R",
    "src/modeling/prediction/combined_feature_table.R",
    "src/modeling/prediction/cluster_pred_table.R",
    "src/desc/dev_attainment.R",
    "src/modeling/modeling.R",
    "src/modeling/plot_model_coefficients.R",
    "src/changepoints/change_point_analysis.R",
    "src/changepoints/cluster_fluctuation.R",
    "src/changepoints/cluster_thresholds.R",
    "src/periodicity/patient_level_permutation.R",
    "src/periodicity/patient_level_periodicity.R",
    "src/periodicity/plot_patient_level_periodicity.R",
    "src/periodicity/plot_patient_level_permutation.R",
    "src/episodicity/patient_level_episodicity.R",
    "src/episodicity/plot_patient_level_episodicity.R",
    "src/timelines/timelines.R",
    "src/timelines/four_patient_timeline.R"
  ),
  stringsAsFactors = FALSE
)

all_src_scripts <- sort(normalizePath(
  list.files("src", pattern = "\\.R$", recursive = TRUE, full.names = TRUE),
  mustWork = TRUE
))
runner_script <- normalizePath("src/run_all.R", mustWork = TRUE)
manifest_scripts <- normalizePath(pipeline_steps$script, mustWork = FALSE)

missing_manifest_scripts <- pipeline_steps$script[!file.exists(pipeline_steps$script)]
if (length(missing_manifest_scripts) > 0) {
  stop(
    "Pipeline manifest references missing script(s):\n  ",
    paste(missing_manifest_scripts, collapse = "\n  "),
    call. = FALSE
  )
}

unlisted_scripts <- setdiff(all_src_scripts, c(manifest_scripts, runner_script))
if (length(unlisted_scripts) > 0) {
  stop(
    "R script(s) under src/ are not listed in pipeline_steps:\n  ",
    paste(unlisted_scripts, collapse = "\n  "),
    "\nAdd them to src/run_all.R in dependency order.",
    call. = FALSE
  )
}

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
rscript_path <- file.path(R.home("bin"), "Rscript")

as_r_string <- function(x) {
  paste0('"', gsub('(["\\\\])', '\\\\\\1', x), '"')
}

run_step <- function(index, total, phase, script) {
  label <- sprintf("[%02d/%02d] %s: %s", index, total, phase, script)
  message(timestamp(), " START ", label)

  elapsed <- system.time({
    r_expr <- paste0(
      "options(warn = -1); ",
      "source(",
      as_r_string(script),
      ", chdir = FALSE)"
    )
    status <- system2(rscript_path, args = c("--vanilla", "-e", shQuote(r_expr)))
  })

  if (!identical(status, 0L)) {
    stop("Step exited with status ", status, call. = FALSE)
  }

  message(
    timestamp(),
    " DONE  ",
    label,
    sprintf(" (%.1fs)", unname(elapsed[["elapsed"]]))
  )
}

message("Project root: ", project_root)
message("Running ", nrow(pipeline_steps), " R scripts.")
message("continue_on_error: ", continue_on_error)

if (dry_run) {
  message("dry_run: TRUE")
  for (i in seq_len(nrow(pipeline_steps))) {
    step <- pipeline_steps[i, ]
    message(sprintf("[%02d/%02d] %s: %s", i, nrow(pipeline_steps), step$phase, step$script))
  }
  quit(status = 0)
}

failures <- list()
for (i in seq_len(nrow(pipeline_steps))) {
  step <- pipeline_steps[i, ]
  tryCatch(
    run_step(i, nrow(pipeline_steps), step$phase, step$script),
    error = function(err) {
      message(timestamp(), " FAIL  ", step$script)
      message(conditionMessage(err))
      failures[[length(failures) + 1L]] <<- list(script = step$script, error = conditionMessage(err))

      if (!continue_on_error) {
        stop(err)
      }
    }
  )
}

if (length(failures) > 0) {
  failure_lines <- vapply(
    failures,
    function(x) paste0(x$script, ": ", x$error),
    character(1)
  )
  stop(
    "Pipeline completed with failure(s):\n  ",
    paste(failure_lines, collapse = "\n  "),
    call. = FALSE
  )
}

message(timestamp(), " Pipeline complete.")
