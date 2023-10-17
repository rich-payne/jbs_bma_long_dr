gen_data <- function(
  n_per_arm,
  scenario,
  enroll_rate,
  times = c(0, 4, 8, 12, 24),
  doses = c(0, 1, 5, 10, 15),
  sigma = 1.5
) {
  t_max <- max(times)
  a <- .25
  n_doses <- length(doses)
  n_times <- length(times)
  if (scenario == "null-null") {
    dat <- sim_patients(
      times = times,
      doses = doses,
      intercept = a,
      response_doses = rep(0, n_doses),
      response_times = rep(0, n_times),
      n_per_arm = n_per_arm,
      sigma = sigma
    )
  } else if (scenario == "plateau-plateau") {
    dat <- sim_patients(
      times = times,
      doses = doses,
      intercept = a,
      response_doses = c(1, 2.5, 3.9, 4, 4),
      response_times = c(0, .5, .9, 1, 1),
      n_per_arm = n_per_arm,
      sigma = sigma
    )
  } else if (scenario == "plateau-decrease") {
    dat <- sim_patients(
      times = times,
      doses = doses,
      intercept = a,
      response_doses = c(1, 2.5, 3.9, 4, 4) / .90,
      response_times = c(0, .5, .9, 1, .90),
      n_per_arm = n_per_arm,
      sigma = sigma
    )
  } else if (scenario == "emax-itp") {
    dat <- dreamer_data_emax(
      n_cohorts = rep(n_per_arm, n_doses),
      dose = doses,
      b1 = 1,
      b2 = 4.1,
      b3 = 1,
      b4 = 2,
      sigma = sigma,
      times = times,
      longitudinal = "itp",
      a = .25,
      c1 = .25,
      t_max = t_max
    )
  } else if (scenario == "emax-idp") {
    dat <- dreamer_data_emax(
      n_cohorts = rep(n_per_arm, n_doses),
      dose = doses,
      b1 = 1,
      b2 = 5.15,
      b3 = 1,
      b4 = 2,
      sigma = sigma,
      times = times,
      longitudinal = "idp",
      a = .25,
      c1 = .25,
      c2 = -.25,
      d1 = 10,
      d2 = 15,
      gam = .2,
      t_max = t_max
    )
  }
  dat <- dat %>%
    mutate(arm = as.character(dose)) %>%
    enroll_patients(doses, enroll_rate)
  return(dat)
}

logit <- function(x) log(x / (1 - x))
ilogit <- function(x) 1 / (1 + exp(- x))

gen_data_binary <- function(
  n_per_arm,
  scenario,
  enroll_rate,
  times = c(0, 4, 8, 12, 24),
  doses = c(0, 1, 5, 10, 15)
) {
  t_max <- max(times)
  a <- logit(.10)
  n_doses <- length(doses)
  n_times <- length(times)
  if (scenario == "null-null") {
    dat <- sim_patients(
      times = times,
      doses = doses,
      intercept = a,
      response_doses = rep(0, n_doses),
      response_times = rep(0, n_times),
      n_per_arm = n_per_arm,
      sigma = 0
    ) %>%
      make_binary()
  } else if (scenario == "plateau-plateau") {
    dat <- sim_patients(
      times = times,
      doses = doses,
      intercept = a,
      response_doses = c(1, 2.5, 3.9, 4, 4) / 1.5,
      response_times = c(0, .5, .9, 1, 1),
      n_per_arm = n_per_arm,
      sigma = 0
    ) %>%
      make_binary()
  } else if (scenario == "plateau-decrease") {
    dat <- sim_patients(
      times = times,
      doses = doses,
      intercept = a,
      response_doses = c(1, 2.5, 3.9, 4, 4) / (.90 * 1.5),
      response_times = c(0, .5, .9, 1, .90),
      n_per_arm = n_per_arm,
      sigma = 0
    ) %>%
      make_binary()
  } else if (scenario == "emax-itp") {
    dat <- dreamer_data_emax(
      n_cohorts = rep(n_per_arm, n_doses),
      dose = doses,
      b1 = 1,
      b2 = 4.08 / 1.5,
      b3 = 1,
      b4 = 2,
      sigma = 0,
      times = times,
      longitudinal = "itp",
      a = a,
      c1 = .25,
      t_max = t_max
    ) %>%
      mutate(mean = response) %>%
      make_binary()
  } else if (scenario == "emax-idp") {
    dat <- dreamer_data_emax(
      n_cohorts = rep(n_per_arm, n_doses),
      dose = doses,
      b1 = 1,
      b2 = 4.09 / (1.5 * 0.8),
      b3 = 1,
      b4 = 2,
      sigma = 0,
      times = times,
      longitudinal = "idp",
      a = a,
      c1 = .25,
      c2 = -.25,
      d1 = 10,
      d2 = 15,
      gam = .2,
      t_max = t_max
    ) %>%
      mutate(mean = response) %>%
      make_binary()
  }
  dat <- dat %>%
    mutate(arm = as.character(dose)) %>%
    enroll_patients(doses, enroll_rate)
  return(dat)
}

sim_patients <- function(
  times,
  doses,
  intercept,
  response_doses,
  response_times,
  n_per_arm,
  sigma
) {
  check_lengths(times, response_times)
  check_lengths(doses, response_doses)
  dat <- expand.grid(time = times, dose = doses) %>%
    mutate(
      time_ind = sapply(time, function(xx) which(xx == !!times)),
      dose_ind = sapply(dose, function(xx) which(xx == !!doses)),
      mean = !!intercept +
        (!!response_times)[time_ind] * (!!response_doses)[dose_ind]
    ) %>%
    select(-time_ind, -dose_ind) %>%
    tidyr::expand_grid(subject_in_arm = 1:n_per_arm) %>%
    mutate(
      subject = paste0(dose, "_", subject_in_arm),
      response = rnorm(n(), mean, !!sigma),
      arm = as.character(dose)
    ) %>%
    dplyr::select(-subject_in_arm)
  return(dat)
}

make_binary <- function(x) {
  x %>%
    mutate(
      mean = ilogit(mean),
      response = rbinom(n(), 1, mean)
    )
}

check_lengths <- function(x, y) {
  if (length(x) != length(y)) {
    stop(
      "times/response_times or doses/response_doses lengths differ.",
      call. = FALSE
    )
  }
}

enroll_patients <- function(dat, doses, enroll_rate) {
  dat %>%
    simpatico::randomize(block_ratios = get_block_ratios(doses)) %>%
    simpatico::add_chronology(enroll_rate = enroll_rate, time_var = "time")
}

get_block_ratios <- function(doses) {
  out <- rep(1, length(doses))
  names(out) <- as.character(doses)
  return(out)
}

summarize_data <- function(x) {
  x %>%
    mutate(n_enrolled = length(unique(subject))) %>%
    group_by(dose, time, n_enrolled) %>%
    summarize(
      obs_mean = mean(response),
      obs_sd = sd(response),
      obs_se = obs_sd / sqrt(n()),
      .groups = "keep"
    ) %>%
    pivot_wider(
      names_from = c("dose", "time"),
      values_from = c("obs_mean", "obs_sd", "obs_se"),
      names_glue = "{.value}_dose_{dose}_time_{time}"
    ) %>%
    ungroup()
}
