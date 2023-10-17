library(targets)
library(tarchetypes)
targets::tar_option_set(
  packages = c(
    "targets",
    "tarchetypes",
    "dreamer",
    "dplyr",
    "tidyr",
    "rlang",
    "purrr",
    "ggplot2",
    "stringr",
    "tibble",
    "fs"
  ),
  storage = "worker",
  retrieval = "worker"
)

options(clustermq.scheduler = "sge")
options(clustermq.template = "sge.tmpl")
options(crayon.enabled = FALSE)

source("R/data_analysis.R")
source("R/data_generation.R")
source("R/eoi.R")
source("R/run_trial.R")
source("R/scenarios.R")
source("R/operating_characteristics.R")

list(
  tar_map_rep(
    sims,
    run_trial(
      n_per_arm = n_per_arm,
      scenario = scenario,
      p_interim = p_interim,
      enroll_rate = enroll_rate
    ),
    values = scenarios,
    reps = 100,
    batches = 10
  ),
  tar_map_rep(
    sims_binary,
    run_trial_binary(
      n_per_arm = n_per_arm,
      scenario = scenario,
      p_interim = p_interim,
      enroll_rate = enroll_rate
    ),
    values = scenarios_binary,
    reps = 10,
    batches = 100
  ),
  tar_target(
    ocs,
    get_ocs(
      sims,
      interim_eoi = "interim_pr_eoi_2.75_arm1_15_arm2_0_time_24_bma",
      final_eoi = "final_pr_eoi_2.75_arm1_15_arm2_0_time_24_bma",
      prth1_int = 0.3,
      prth2_int = 0.95,
      prth_final = 0.8
    )
  ),
  tar_target(
    ocs_binary,
    get_ocs(
      sims_binary,
      interim_eoi = "interim_pr_eoi_0.35_arm1_15_arm2_0_time_24_bma",
      final_eoi = "final_pr_eoi_0.35_arm1_15_arm2_0_time_24_bma",
      prth1_int = 0.3,
      prth2_int = 0.95,
      prth_final = 0.8
    )
  ),
  tar_target(
    ocs_partial,
    get_ocs(
      sims,
      interim_eoi = "interim_partial_pr_eoi_2.75_arm1_15_arm2_0_time_24_bma",
      final_eoi = "final_pr_eoi_2.75_arm1_15_arm2_0_time_24_bma",
      prth1_int = 0.3,
      prth2_int = 0.95,
      prth_final = 0.8
    )
  ),
  tar_target(
    ocs_partial_binary,
    get_ocs(
      sims_binary,
      interim_eoi = "interim_partial_pr_eoi_0.35_arm1_15_arm2_0_time_24_bma",
      final_eoi = "final_pr_eoi_0.35_arm1_15_arm2_0_time_24_bma",
      prth1_int = 0.3,
      prth2_int = 0.95,
      prth_final = 0.8
    )
  ),
  tar_target(
    ocs_no_long,
    get_ocs(
      sims,
      interim_eoi = "interim_no_long_pr_eoi_2.75_arm1_15_arm2_0_time_none_bma",
      final_eoi = "final_no_long_pr_eoi_2.75_arm1_15_arm2_0_time_none_bma",
      prth1_int = 0.3,
      prth2_int = 0.95,
      prth_final = 0.8
    )
  ),
  tar_target(
    ocs_no_long_binary,
    get_ocs(
      sims_binary,
      interim_eoi = "interim_no_long_pr_eoi_0.35_arm1_15_arm2_0_time_none_bma",
      final_eoi = "final_no_long_pr_eoi_0.35_arm1_15_arm2_0_time_none_bma",
      prth1_int = 0.3,
      prth2_int = 0.95,
      prth_final = 0.8
    )
  ),
  tar_target(
    ocs_plot,
    create_plot(ocs, ocs_partial, ocs_no_long, p_interim = .5, enroll_rate = 10)
  ),
  tar_target(
    ocs_plot_binary,
    create_plot(
      ocs_binary,
      ocs_partial_binary,
      ocs_no_long_binary,
      p_interim = .5,
      enroll_rate = 10
    )
  ),
  tar_target(
    weights_plot,
    plot_weights(
      sims = sims,
      n_per_arm = 50,
      p_interim = 0.5,
      enroll_rate = 10,
      milestone = "final"
    )
  ),
  tar_target(
    weights_plot_binary,
    plot_weights(
      sims = sims_binary,
      n_per_arm = 50,
      p_interim = 0.5,
      enroll_rate = 10,
      milestone = "final"
    )
  ),
  tar_target(
    mse_plot,
    plot_mse(
      mse = mse,
      n_per_arm = 50,
      p_interim = 0.5,
      enroll_rate = 10,
      milestone = "final",
      time = NULL,
      dose = NULL
    )
  ),
  tar_target(
    mse_plot_binary,
    plot_mse(
      mse = mse_binary,
      n_per_arm = 50,
      p_interim = 0.5,
      enroll_rate = 10,
      milestone = "final",
      time = NULL,
      dose = NULL
    )
  ),
  tar_target(
    mse_width_plot,
    plot_mse_and_widths(
      mse = mse,
      int_widths = int_widths,
      n_per_arm = 50,
      p_interim = 0.5,
      enroll_rate = 10,
      milestone = "final",
      dose = 15,
      time = 24
    )
  ),
  tar_target(
    mse_width_plot_binary,
    plot_mse_and_widths(
      mse = mse_binary,
      int_widths = int_widths_binary,
      n_per_arm = 50,
      p_interim = 0.5,
      enroll_rate = 10,
      milestone = "final",
      dose = 15,
      time = 24
    )
  ),
  tar_target(mse, calc_mse(sims)),
  tar_target(mse_binary, calc_mse_binary(sims_binary)),
  tar_target(int_widths, get_int_widths(sims)),
  tar_target(int_widths_binary, get_int_widths(sims_binary)),
  tar_target(scenarios_plot, plot_scenarios(sims)),
  tar_target(scenarios_plot_binary, plot_scenarios_binary(sims_binary)),
  tar_target(
    ocs_plot_file,
    save_gg(
      "figures/ocs.eps",
      ocs_plot,
      width = 7.5,
      height = 3.75
    )
  ),
  tar_target(
    ocs_plot_binary_file,
    save_gg(
      "figures/ocs_binary.eps",
      ocs_plot_binary,
      width = 7.5,
      height = 3.75
    )
  ),
  tar_target(
    weights_plot_file,
    save_gg(
      "figures/weights_plot.eps",
      weights_plot,
      width = 7.5,
      height = 3.75
    ),
    format = "file"
  ),
  tar_target(
    mse_plot_file,
    save_gg(
      "figures/mse_plot.eps",
      mse_plot,
      width = 7.5,
      height = 3.75
    ),
    format = "file"
  ),
  tar_target(
    mse_width_plot_file,
    save_gg(
      "figures/mse_width_plot.eps",
      mse_width_plot,
      width = 7.5,
      height = 3.75
    ),
    format = "file"
  ),
  tar_target(
    scenarios_plot_file,
    save_gg(
      "figures/scenarios_plot.eps",
      scenarios_plot,
      width = 7.5,
      height = 3.75
    ),
    format = "file"
  ),
  tar_target(
    scenarios_plot_binary_file,
    save_gg(
      "figures/scenarios_plot_binary.eps",
      scenarios_plot_binary,
      width = 7.5,
      height = 3.75
    ),
    format = "file"
  )
)
