run_trial <- function(n_per_arm, scenario, p_interim, enroll_rate) {
  data <-  gen_data(n_per_arm, scenario, enroll_rate)
  int1_results_partial <- run_analysis(
    data,
    p_interim = p_interim,
    prefix = "interim_partial",
    partial_subject = TRUE
  )
  int1_results <- run_analysis(
    data,
    p_interim = p_interim,
    prefix = "interim",
    partial_subject = FALSE
  )
  int1_results_no_long <- run_analysis_no_long(
    data,
    p_interim = p_interim,
    prefix = "interim_no_long",
    time = max(data$time)
  )
  final_results <- run_analysis(data, p_interim = NULL, prefix = "final")
  final_results_no_long <- run_analysis_no_long(
    data,
    prefix = "final_no_long",
    time = max(data$time)
  )
  bind_cols(
    int1_results,
    int1_results_partial,
    int1_results_no_long,
    final_results,
    final_results_no_long
  )
}

run_trial_binary <- function(n_per_arm, scenario, p_interim, enroll_rate) {
  data <-  gen_data_binary(n_per_arm, scenario, enroll_rate)
  eois <- c(0.3, 0.35, 0.375, 0.40)
  int1_results_partial <- run_analysis_binary(
    data,
    p_interim = p_interim,
    prefix = "interim_partial",
    partial_subject = TRUE,
    eois = eois
  )
  int1_results <- run_analysis_binary(
    data,
    p_interim = p_interim,
    prefix = "interim",
    partial_subject = FALSE,
    eois = eois
  )
  int1_results_no_long <- run_analysis_no_long_binary(
    data,
    p_interim = p_interim,
    prefix = "interim_no_long",
    time = max(data$time),
    eois = eois
  )
  final_results <- run_analysis_binary(
    data,
    p_interim = NULL,
    prefix = "final",
    eois = eois
  )
  final_results_no_long <- run_analysis_no_long_binary(
    data,
    prefix = "final_no_long",
    time = max(data$time),
    eois = eois
  )
  bind_cols(
    int1_results,
    int1_results_partial,
    int1_results_no_long,
    final_results,
    final_results_no_long
  )
}
