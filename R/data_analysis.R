run_analysis <- function(data, p_interim = NULL, prefix = "", partial_subjects) {
  if (!is.null(p_interim)) {
    n_interim <- round(p_interim * length(unique(data$subject)))
    data <- data %>%
      select_subjects_by_n(
        n = n_interim,
        partial_subjects = partial_subjects,
        visit_time_var = "time"
      )
  }
  mcmc <- run_mcmc(data)
  mcmc_sumry <- summarize_mcmc(mcmc)
  post_sumry <- summarize_post(mcmc)
  weight_sumry <- summarize_weights(mcmc)
  eoi_probs <- get_eois(mcmc, doses = mcmc$doses)
  data_summary <- summarize_data(data)
  out <- bind_cols(mcmc_sumry, post_sumry, weight_sumry, eoi_probs, data_summary)
  colnames(out) <- paste0(format_prefix(prefix), colnames(out))
  return(out)
}

run_analysis_binary <- function(
  data,
  p_interim = NULL,
  prefix = "",
  partial_subjects,
  eois
) {
  if (!is.null(p_interim)) {
    n_interim <- round(p_interim * length(unique(data$subject)))
    data <- data %>%
      select_subjects_by_n(
        n = n_interim,
        partial_subjects = partial_subjects,
        visit_time_var = "time"
      )
  }
  mcmc <- run_mcmc_binary(data)
  mcmc_sumry <- summarize_mcmc(mcmc)
  post_sumry <- summarize_post(mcmc)
  weight_sumry <- summarize_weights(mcmc)
  eoi_probs <- get_eois(mcmc, doses = mcmc$doses, eois)
  data_summary <- summarize_data(data)
  out <- bind_cols(mcmc_sumry, post_sumry, weight_sumry, eoi_probs, data_summary)
  colnames(out) <- paste0(format_prefix(prefix), colnames(out))
  return(out)
}


run_analysis_no_long <- function(data, p_interim = NULL, prefix = "", time) {
  if (!is.null(p_interim)) {
    n_interim <- round(p_interim * length(unique(data$subject)))
    data <- data %>%
      select_subjects_by_n(
        n = n_interim,
        partial_subjects = FALSE,
        visit_time_var = "time"
      )
  }
  data <- data %>% dplyr::filter(time == !!time)
  mcmc <- run_mcmc_no_long(data)
  mcmc_sumry <- summarize_mcmc(mcmc)
  post_sumry <- summarize_post(mcmc)
  weight_sumry <- summarize_weights(mcmc)
  eoi_probs <- get_eois(mcmc, doses = mcmc$doses)
  data_summary <- summarize_data(data)
  out <- bind_cols(mcmc_sumry, post_sumry, weight_sumry, eoi_probs, data_summary)
  colnames(out) <- paste0(format_prefix(prefix), colnames(out))
  return(out)
}

run_analysis_no_long_binary <- function(
  data,
  p_interim = NULL,
  prefix = "",
  time,
  eois
) {
  if (!is.null(p_interim)) {
    n_interim <- round(p_interim * length(unique(data$subject)))
    data <- data %>%
      simpatico::select_subjects_by_n(
        n = n_interim,
        partial_subjects = FALSE,
        visit_time_var = "time"
      )
  }
  data <- data %>% dplyr::filter(time == !!time)
  mcmc <- run_mcmc_no_long_binary(data)
  mcmc_sumry <- summarize_mcmc(mcmc)
  post_sumry <- summarize_post(mcmc)
  weight_sumry <- summarize_weights(mcmc)
  eoi_probs <- get_eois(mcmc, doses = mcmc$doses, eois)
  data_summary <- summarize_data(data)
  out <- bind_cols(mcmc_sumry, post_sumry, weight_sumry, eoi_probs, data_summary)
  colnames(out) <- paste0(format_prefix(prefix), colnames(out))
  return(out)
}

format_prefix <- function(x) {
  if(x != "") {
    x <- paste0(x, "_")
  }
  return(x)
}

run_mcmc <- function(data) {
  models <- get_models(t_max = max(data$time))
  mcmc <- rlang::exec(
    dreamer_mcmc,
    data = data,
    n_chains = 4,
    n_iter = 5000,
    silent = TRUE,
    convergence_warn = FALSE,
    !!!models
  )
}

run_mcmc_binary <- function(data) {
  models <- get_models_binary(t_max = max(data$time))
  mcmc <- rlang::exec(
    dreamer_mcmc,
    data = data,
    n_chains = 4,
    n_iter = 5000,
    silent = TRUE,
    convergence_warn = FALSE,
    !!!models
  )
}

run_mcmc_no_long <- function(data) {
  models <- get_models(t_max = max(data$time), longitudinal = FALSE)
  mcmc <- rlang::exec(
    dreamer_mcmc,
    data = data,
    n_chains = 4,
    n_iter = 5000,
    silent = TRUE,
    convergence_warn = FALSE,
    !!!models
  )
}

run_mcmc_no_long_binary <- function(data) {
  models <- get_models_binary(t_max = max(data$time), longitudinal = FALSE)
  mcmc <- rlang::exec(
    dreamer_mcmc,
    data = data,
    n_chains = 4,
    n_iter = 5000,
    silent = TRUE,
    convergence_warn = FALSE,
    !!!models
  )
}

summarize_weights <- function(x) {
  out <- x$w_post
  names(out) <- paste0("w_post_", names(out))
  out <- data.frame(t(out))
  return(out)
}

summarize_mcmc <- function(x) {
  ind <- mcmc_indices(x)
  map2(x[ind], names(x[ind]), summarize_mcmc_impl) %>%
    bind_cols()
}

summarize_post <- function(x) {
  mod_index <- sapply(x, inherits, "dreamer_mcmc") %>%
    {names(.)[.]} %>%
    c("bma") %>%
    {data.frame(model_name = .)}
  pmap(mod_index, summarize_post_impl, x = x) %>%
    bind_cols()
}

summarize_post_impl <- function(x, model_name) {
  out <- select_model(x, model_name) %>%
    {dreamer::posterior(.)$stats}
  if (rlang::has_name(out, "time")) {
    out <- out %>%
      pivot_wider(
        names_from = c("dose", "time"),
        values_from = c("mean", "2.50%", "97.50%"),
        names_glue = paste0("post_{.value}_dose_{dose}_time_{time}_", model_name)
      )
  } else {
    out <- out %>%
      pivot_wider(
        names_from = "dose",
        values_from = c("mean", "2.50%", "97.50%"),
        names_glue = paste0("post_{.value}_dose_{dose}_time_none_", model_name)
      )
  }
  return(out)
}

select_model <- function(x, model_name) {
  if (model_name == "bma") {
    return(x)
  }
  x[[model_name]]
}

mcmc_indices <- function(x) {
  vapply(
    x,
    function(xx) inherits(xx, "dreamer_mcmc"),
    logical(1)
  ) %>%
    which()
}

# x is a single jags model object
summarize_mcmc_impl <- function(x, model_name) {
  gman <- summary(x) %>%
    select(
      param,
      post_mean = mean,
      post_median = `50%`,
      post_lb = `2.5%`,
      post_ub = `97.5%`,
      gelman_upper
    ) %>%
    pivot_wider(
      names_from = "param",
      values_from = c("gelman_upper", "post_mean", "post_lb", "post_ub", "post_median"),
      names_prefix = paste0(model_name, "_")
    )
}

get_models <- function(t_max, longitudinal = TRUE) {
  mods <- list(
    mod_quad = expr(model_quad(
      mu_b1 = 0,
      sigma_b1 = 10,
      mu_b2 = 0,
      sigma_b2 = 10,
      mu_b3 = 0,
      sigma_b3 = 10,
      shape = 1,
      rate = .01
    )),
    mod_logquad = expr(model_logquad(
      mu_b1 = 0,
      sigma_b1 = 10,
      mu_b2 = 0,
      sigma_b2 = 10,
      mu_b3 = 0,
      sigma_b3 = 10,
      shape = 1,
      rate = .01
    )),
    mod_emax = expr(model_emax(
      mu_b1 = 0,
      sigma_b1 = 10,
      mu_b2 = 0,
      sigma_b2 = 10,
      mu_b3 = log(8),
      sigma_b3 = 1,
      mu_b4 = 1,
      sigma_b4 = 3,
      shape = 1,
      rate = .01
    )),
    mod_exp = expr(model_exp(
      mu_b1 = 0,
      sigma_b1 = 10,
      mu_b2 = 0,
      sigma_b2 = 10,
      mu_b3 = -1,
      sigma_b3 = 5,
      shape = 1,
      rate = .01
    ))
  )
  mods_long <- list(mods_long = NULL)
  if (longitudinal) {
    mods_long <- list(
      long_itp = expr(
        model_longitudinal_itp(
          mu_a = 0, sigma_a = 10, a_c1 = 0, b_c1 = 5, t_max = !!t_max
        )
      ),
      long_idp = expr(
        model_longitudinal_idp(
          mu_a = 0, sigma_a = 10, a_c1 = 0, b_c1 = 5, a_c2 = -5, b_c2 = 0, t_max = !!t_max
        )
      )
    )
  }
  mod_grid <- expand_grid(mods, mods_long)
  w_prior <- 1 / nrow(mod_grid)
  out <- purrr::pmap(
    mod_grid,
    function(mods, mods_long, w_prior) {
      call_modify(mods, longitudinal = mods_long, w_prior = w_prior) %>%
        eval()
    },
    w_prior = w_prior
  )
  nmes <- names(mods_long)
  if (is.null(nmes)) {
    nmes <- "none"
  }
  names(out) <- paste0(names(out), "_", nmes)
  return(out)
}

get_models_binary <- function(t_max, longitudinal = TRUE) {
  mods <- list(
    mod_quad = expr(model_quad_binary(
      mu_b1 = 0,
      sigma_b1 = 2,
      mu_b2 = 0,
      sigma_b2 = 2,
      mu_b3 = 0,
      sigma_b3 = 2,
      link = "logit"
    )),
    mod_logquad = expr(model_logquad_binary(
      mu_b1 = 0,
      sigma_b1 = 2,
      mu_b2 = 0,
      sigma_b2 = 2,
      mu_b3 = 0,
      sigma_b3 = 2,
      link = "logit"
    )),
    mod_emax = expr(model_emax_binary(
      mu_b1 = 0,
      sigma_b1 = 2,
      mu_b2 = 0,
      sigma_b2 = 2,
      mu_b3 = log(8),
      sigma_b3 = 2,
      mu_b4 = 1,
      sigma_b4 = 2,
      link = "logit"
    )),
    mod_exp = expr(model_exp_binary(
      mu_b1 = 0,
      sigma_b1 = 2,
      mu_b2 = 0,
      sigma_b2 = 2,
      mu_b3 = -1,
      sigma_b3 = 2,
      link = "logit"
    ))
  )
  mods_long <- list(mods_long = NULL)
  if (longitudinal) {
    mods_long <- list(
      long_itp = expr(
        model_longitudinal_itp(
          mu_a = 0, sigma_a = 2, a_c1 = 0, b_c1 = 5, t_max = !!t_max
        )
      ),
      long_idp = expr(
        model_longitudinal_idp(
          mu_a = 0, sigma_a = 2, a_c1 = 0, b_c1 = 5, a_c2 = -5, b_c2 = 0, t_max = !!t_max
        )
      )
    )
  }
  mod_grid <- expand_grid(mods, mods_long)
  w_prior <- 1 / nrow(mod_grid)
  out <- purrr::pmap(
    mod_grid,
    function(mods, mods_long, w_prior) {
      call_modify(mods, longitudinal = mods_long, w_prior = w_prior) %>%
        eval()
    },
    w_prior = w_prior
  )
  nmes <- names(mods_long)
  if (is.null(nmes)) {
    nmes <- "none"
  }
  names(out) <- paste0(names(out), "_", nmes)
  return(out)
}

select_subjects_by_n <- function(
    data,
    n,
    partial_subjects = TRUE,
    visit_time_var = "visit_time",
    trial_time_var = "trial_time",
    subject_var = "subject",
    enroll_time_var = "enroll_time"
) {
  visit_time_var <- sym(visit_time_var)
  trial_time_var <- sym(trial_time_var)
  subject_var <- sym(subject_var)
  enroll_time_var <- sym(enroll_time_var)
  finish_time <- get_finish_time(data, subject_var, trial_time_var, n)
  out <- data %>%
    group_by(!!subject_var) %>%
    mutate(partial_data = !all(!!trial_time_var <= !!finish_time)) %>%
    ungroup() %>%
    filter(!!trial_time_var <= !!finish_time) %>%
    remove_partials(partial_subjects)
  return(out)
}

remove_partials <- function(data, partial_subjects) {
  if (!partial_subjects) {
    data <- data %>% dplyr::filter(!.data$partial_data)
  }
  data <- select(data, -.data$partial_data)
  return(data)
}

get_finish_time <- function (data, subject_var, trial_time_var, n) {
  data %>%
    group_by(!!subject_var) %>%
    mutate(last_subject_time = max(!!trial_time_var)) %>%
    ungroup() %>%
    arrange(.data$last_subject_time) %>% group_by(!!subject_var) %>%
    nest() %>%
    ungroup() %>%
    slice(!!n) %>%
    unnest("data") %>%
    slice(1) %>%
    pull(.data$last_subject_time)
}
