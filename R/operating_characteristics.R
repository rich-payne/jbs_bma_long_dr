get_ocs <- function(
  sims,
  interim_eoi,
  final_eoi,
  prth1_int,
  prth2_int,
  prth_final
) {
  interim_eoi <- rlang::sym(interim_eoi)
  final_eoi <- rlang::sym(final_eoi)

  ocs <- select(
    sims,
    !!interim_eoi,
    !!final_eoi,
    all_of(c("scenario", "n_per_arm", "p_interim", "enroll_rate"))
  ) %>%
    group_by(across(all_of(
      c("scenario", "n_per_arm", "p_interim", "enroll_rate")
    ))) %>%
    mutate(
      interim = case_when(
        !!interim_eoi < !!prth1_int ~ "no efficacy",
        !!interim_eoi > !!prth2_int ~ "efficacy",
        TRUE ~ "sufficient efficacy"
      ),
      final = case_when(
        !!final_eoi > !!prth_final ~ "efficacy",
        TRUE ~ "no efficacy"
      ),
      stop = case_when(
        interim != "sufficient efficacy" ~ "interim",
        TRUE ~ "final"
      ),
      stop_outcome = case_when(
        stop == "interim" & interim == "no efficacy" ~ "futility",
        stop == "interim" & interim == "efficacy" ~ "success",
        stop == "final" & final == "no efficacy" ~ "futility",
        stop == "final" & final == "efficacy" ~ "success"
      )
    )
  return(ocs)
}

plot_weights <- function(sims, n_per_arm, p_interim, enroll_rate, milestone) {
  dat_plot <- sims %>%
    dplyr::filter(
      n_per_arm == !!n_per_arm,
      p_interim == !!p_interim,
      enroll_rate == !!enroll_rate
    ) %>%
    select(scenario, contains("w_post")) %>%
    pivot_longer(contains("w_post"), names_to = "model", values_to = "w_post") %>%
    mutate(
      milestone = str_match(model, "(.*)_w")[, 2],
      model = sub(".*mod_", "", model) %>%
        {sub("_long", "", .)}
    ) %>%
    dplyr::filter(milestone == !!milestone)
  scen_levels <- c(
    "null-null",
    "emax-itp",
    "emax-idp",
    "plateau-plateau",
    "plateau-decrease"
  )
  dat_plot_tile <- dat_plot %>%
    mutate(
      dose_response = sub("_.*", "", model),
      longitudinal = sub(".*_", "", model)
    ) %>%
    group_by(scenario, dose_response, longitudinal) %>%
    summarize(w_post = median(w_post)) %>%
    mutate(
      scenario = factor(
        scenario,
        levels = !!scen_levels
      ),
      longitudinal = factor(toupper(longitudinal), levels = c("ITP", "IDP"))
    )

  ggplot(dat_plot_tile, aes(longitudinal, dose_response, fill = w_post)) +
    geom_tile() +
    facet_wrap(vars(scenario)) +
    geom_text(aes(label = sprintf("%.2f", w_post))) +
    scale_fill_gradient(low = "white", high = "orange3") +
    ylab("Dose response") +
    xlab("Longitudinal") +
    labs(fill = "Median weight") +
    ggtitle("Posterior Weights") +
    theme(plot.title = element_text(hjust = .5))
}

calc_mse <- function(sims) {
  scenarios <- sims %>% expand(scenario)
  true_means <- pmap_df(
    scenarios,
    gen_data,
    n_per_arm = 1,
    enroll_rate = 1,
    sigma = 0,
    .id = "scenario"
  ) %>%
    mutate(
      scenario = sapply(
        as.numeric(scenario),
        function(xx) ((!!scenarios)$scenario)[xx]
      )
    ) %>%
    select(scenario, dose, time, true_mean = response)

  out <- sims %>%
    select(scenario, n_per_arm, p_interim, enroll_rate, contains("post_mean_dose")) %>%
    pivot_longer(contains("post_mean_dose"), values_to = "mean") %>%
    mutate(
      dose = as.numeric(str_match(name, "dose_(.*)_time")[, 2]),
      time = as.numeric(str_match(name, "time_(.*?)_")[, 2]),
      model = str_match(name, "_time_.*?_(?:mod_|)(.*)")[, 2],
      milestone = str_match(name, "(.*)_post")[, 2]
    ) %>%
    left_join(true_means, by = c("scenario", "dose", "time")) %>%
    mutate(
      difference = mean - true_mean
    )
  out_all <- out %>%
    group_by(scenario, n_per_arm, p_interim, enroll_rate, milestone, model) %>%
    summarize(
      mse = mean(difference ^ 2),
      bias = mean(difference)
    )
  out_times <- out %>%
    group_by(scenario, n_per_arm, p_interim, enroll_rate, milestone, model, time, dose) %>%
    summarize(
      mse = mean(difference ^ 2),
      bias = mean(difference)
    )
  dplyr::bind_rows(out_all, out_times)
}

calc_mse_binary <- function(sims) {
  scenarios <- sims %>% expand(scenario)
  true_means <- pmap_df(
    scenarios,
    gen_data_binary,
    n_per_arm = 1,
    enroll_rate = 1,
    .id = "scenario"
  ) %>%
    mutate(
      scenario = sapply(
        as.numeric(scenario),
        function(xx) ((!!scenarios)$scenario)[xx]
      )
    ) %>%
    select(scenario, dose, time, true_mean = mean)

  out <- sims %>%
    select(scenario, n_per_arm, p_interim, enroll_rate, contains("post_mean_dose")) %>%
    pivot_longer(contains("post_mean_dose"), values_to = "mean") %>%
    dplyr::filter(!grepl("time_none", name)) %>%
    mutate(
      dose = as.numeric(str_match(name, "dose_(.*)_time")[, 2]),
      time = as.numeric(str_match(name, "time_(.*?)_")[, 2]),
      model = str_match(name, "_time_.*?_(?:mod_|)(.*)")[, 2],
      milestone = str_match(name, "(.*)_post")[, 2]
    ) %>%
    left_join(true_means, by = c("scenario", "dose", "time")) %>%
    mutate(
      difference = mean - true_mean
    )
  out_all <- out %>%
    group_by(scenario, n_per_arm, p_interim, enroll_rate, milestone, model) %>%
    summarize(
      mse = mean(difference ^ 2),
      bias = mean(difference),
      .groups = "drop"
    )
  out_times <- out %>%
    group_by(
      scenario,
      n_per_arm,
      p_interim,
      enroll_rate,
      milestone,
      model,
      time,
      dose
    ) %>%
    summarize(
      mse = mean(difference ^ 2),
      bias = mean(difference),
      .groups = "drop"
    )
  dplyr::bind_rows(out_all, out_times)
}

plot_scenarios <- function(sims) {
  scenarios <- sims %>% expand(scenario)
  true_means <- pmap_df(
    scenarios,
    gen_data,
    n_per_arm = 1,
    enroll_rate = 1,
    sigma = 0,
    .id = "scenario"
  ) %>%
    mutate(
      scenario = sapply(
        as.numeric(scenario),
        function(xx) ((!!scenarios)$scenario)[xx]
      )
    ) %>%
    select(scenario, dose, time, true_mean = response) %>%
    mutate(
      dose = factor(dose, levels = c("15", "10", "5", "1", "0")),
      scenario = factor(
        scenario,
        c("null-null", "plateau-plateau", "plateau-decrease", "emax-itp", "emax-idp")
      )
    )

  ggplot(true_means, aes(time, true_mean, group = dose, color = dose)) +
    facet_wrap(~scenario) +
    geom_line() +
    geom_point() +
    labs(
      x = "Week",
      y = "Mean response",
      color = "Dose"
    ) +
    scale_x_continuous(breaks = unique(true_means$time))
}

plot_scenarios_binary <- function(sims) {
  scenarios <- sims %>% expand(scenario)
  true_means <- pmap_df(
    scenarios,
    gen_data_binary,
    n_per_arm = 1,
    enroll_rate = 1,
    .id = "scenario"
  ) %>%
    mutate(
      scenario = sapply(
        as.numeric(scenario),
        function(xx) ((!!scenarios)$scenario)[xx]
      )
    ) %>%
    select(scenario, dose, time, true_mean = mean) %>%
    mutate(
      dose = factor(dose, levels = c("15", "10", "5", "1", "0")),
      scenario = factor(scenario, levels = c("null-null", "plateau-plateau", "plateau-decrease", "emax-itp", "emax-idp"))
    )

  ggplot(true_means, aes(time, true_mean, group = dose, color = dose)) +
    facet_wrap(~scenario) +
    geom_line() +
    geom_point() +
    labs(
      x = "Week",
      y = "Mean response",
      color = "Dose"
    ) +
    scale_x_continuous(breaks = unique(true_means$time)) +
    scale_y_continuous(breaks = seq(0, 1, .1), limits = c(0, 1))
}

get_int_widths <- function(sims) {
  out <- sims %>%
    select(
      matches(".*_post_.*dose"),
      scenario,
      n_per_arm,
      p_interim,
      enroll_rate,
      tar_batch,
      tar_rep
    ) %>%
    pivot_longer(
      cols = matches(".*_post_.*dose")
    ) %>%
    mutate(
      dose = as.numeric(str_match(name, "dose_(.*)_time")[, 2]),
      time = as.numeric(str_match(name, "time_(.*?)_")[, 2]),
      quant = str_match(name, "post_(.*?)_")[, 2],
      model = str_match(name, "_time_.*?_(?:mod_|)(.*)")[, 2],
      milestone = str_match(name, "(.*)_post")[, 2]
    )
  out2 <- out %>%
    pivot_wider(
      id_cols = c(
        "model",
        "dose",
        "time",
        "scenario",
        "n_per_arm",
        "p_interim",
        "enroll_rate",
        "milestone",
        "tar_batch",
        "tar_rep"
      ),
      names_from = "quant",
      values_from = "value"
    )

  out3 <- out2 %>%
    mutate(width = `97.50%` - `2.50%`) %>%
    group_by(scenario, n_per_arm, p_interim, enroll_rate, milestone, model, time, dose) %>%
    summarize(median_width = median(width)) %>%
    ungroup() %>%
    mutate(
      scen_dose_model = str_match(scenario, "(.*?)-")[, 2],
      scen_long_model = str_match(scenario, "-(.*)")[, 2]
    ) %>%
    mutate(
      correct_dose_model = map2_lgl(scen_dose_model, model, grepl),
      correct_long_model = map2_lgl(scen_long_model, model, grepl),
    )
  return(out3)
}

save_gg <- function(filename, plot, ...) {
  path <- path_wd(filename)
  ggsave(path, plot, ...)
  return(path)
}


plot_mse <- function(mse, n_per_arm, p_interim, enroll_rate, milestone, time, dose) {
  scenario_levels <- c(
    "null-null",
    "emax-itp",
    "emax-idp",
    "plateau-plateau",
    "plateau-decrease"
  )
  dat_plot <- mse %>%
    dplyr::filter(
      n_per_arm == !!n_per_arm,
      p_interim == !!p_interim,
      enroll_rate == !!enroll_rate,
      milestone == !!milestone
    ) %>%
    mutate(
      dose_response = sub("_.*", "", model),
      longitudinal = sub(".*_", "", model),
      scenario = factor(
        scenario,
        levels = scenario_levels
      )
    )
  if (!is.null(time)) {
    dat_plot <- dplyr::filter(dat_plot, time == !!time)
  } else {
    dat_plot <- dplyr::filter(dat_plot, is.na(time))
  }
  if (!is.null(dose)) {
    dat_plot <- dplyr::filter(dat_plot, dose == !!dose)
  } else {
    dat_plot <- dplyr::filter(dat_plot, is.na(dose))
  }
  ggplot(dat_plot, aes(longitudinal, dose_response, fill = mse)) +
    geom_tile() +
    facet_wrap(vars(scenario)) +
    geom_text(aes(label = sprintf("%.3f", mse))) +
    scale_fill_gradient(low = "orange3", high = "white") +
    ylab("Dose response") +
    xlab("Longitudinal") +
    labs(fill = "MSE") +
    ggtitle("Posterior Mean-Squared-Error") +
    theme(plot.title = element_text(hjust = .5))
}

plot_mse_and_widths <- function(
  mse,
  int_widths,
  n_per_arm,
  p_interim,
  enroll_rate,
  milestone,
  dose,
  time
) {
  scenario_levels <- c(
    "null-null",
    "emax-itp",
    "emax-idp",
    "plateau-plateau",
    "plateau-decrease"
  )
  dat_mse <- mse %>%
    dplyr::filter(
      n_per_arm == !!n_per_arm,
      p_interim == !!p_interim,
      enroll_rate == !!enroll_rate,
      milestone == !!milestone,
      time == !!time,
      dose == !!dose
    )
  dat_widths <- int_widths %>%
    dplyr::filter(
      n_per_arm == !!n_per_arm,
      p_interim == !!p_interim,
      enroll_rate == !!enroll_rate,
      milestone == !!milestone,
      time == !!time,
      dose == !!dose
    )
  dat <- dplyr::full_join(
    dat_mse,
    dat_widths,
    by = c(
      "scenario",
      "n_per_arm",
      "p_interim",
      "enroll_rate",
      "milestone",
      "model",
      "time",
      "dose"
    )
  ) %>%
    mutate(
      width = median_width,
      Model = case_when(
        model == "bma" ~ "BMA",
        scenario == "emax-itp" & model == "emax_long_itp" ~ "True",
        scenario == "emax-idp" & model == "emax_long_idp" ~ "True",
        TRUE ~ "Other"
      ),
      Model = factor(Model, levels = c("BMA", "True", "Other")),
      scenario = factor(
        scenario,
        levels = scenario_levels
      )
    )
  ggplot(dat, aes(width, mse, color = Model, shape = Model)) +
    geom_point(size = 4) +
    facet_wrap(vars(scenario)) +
    scale_shape_manual(values = c(17, 15, 1)) +
    labs(
      y = "MSE",
      x = "Interval Width"
    )
}

summarize_stop_outcome <- function(x, stop_outcome) {
  out <- x %>%
    dplyr::mutate(n_scenario = n()) %>%
    group_by(stop, .add = TRUE) %>%
    summarize(
      prob = sum(stop_outcome == !!stop_outcome) / n_scenario[1],
      .groups = "drop"
    ) %>%
    expand_stops(c("interim", "final")) %>%
    dplyr::mutate(stop_outcome = !!stop_outcome)
  return(out)
}

expand_stops <- function(x, milestone_names) {
  grps <- groups(x)
  x %>%
    ungroup() %>%
    mutate(stop = factor(stop, levels = milestone_names)) %>%
    complete(!!!grps, stop, fill = list(prob = 0)) %>%
    mutate(stop = as.character(stop)) %>%
    group_by(!!!grps)
}

create_plot <- function(ocs, ocs_partial, ocs_no_long, p_interim, enroll_rate, remove_null = TRUE) {
  data_plot <- bind_rows(
    summarize_stop_outcome(ocs, "success") %>%
      add_column(analysis = "only completers"),
    summarize_stop_outcome(ocs_partial, "success") %>%
      add_column(analysis = "all enrolled"),
    summarize_stop_outcome(ocs_no_long, "success") %>%
      add_column(analysis = "BMA-Mod")
  ) %>%
    dplyr::filter(
      p_interim == !!p_interim,
      enroll_rate == !!enroll_rate
    )
  if (remove_null) {
    data_plot <- data_plot %>% dplyr::filter(scenario != "null-null")
  }
  data_plot <- data_plot %>%
    mutate(
      analysis = factor(
        analysis,
        levels = c("BMA-Mod", "only completers", "all enrolled")
      ),
      stop = factor(stop, levels = c("interim", "final"))
    ) %>%
    rename(`n per arm` = n_per_arm) %>%
    pivot_wider(
      names_from = "stop",
      values_from = "prob"
    ) %>%
    mutate(
      total = interim + final
    ) %>%
    pivot_longer(
      interim:total,
      names_to = "stop",
      values_to = "prob"
    ) %>%
    dplyr::filter(stop != "final")

  p <- ggplot(data_plot, aes(stop, prob, group = analysis, fill = analysis)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(
      `n per arm` ~ scenario,
      labeller = labeller(`n per arm` = label_both, scenario = label_value)
    ) +
    scale_fill_grey(start = .7, end = .1) +
    labs(
      fill = "Analysis",
      x = "Milestone",
      y = "Pr(Stop for Success)"
    )
  return(p)
}
