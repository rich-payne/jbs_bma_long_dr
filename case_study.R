library(dreamer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

n_per_arm <- 20
doses <- c(0, 1, 5, 10, 15)
times <- c(0, 4, 8, 12, 16, 20, 24)
t_max <- max(times)

set.seed(2221)

data <- expand_grid(subject = 1:n_per_arm, dose = doses) %>%
  mutate(subject = 1:n()) %>%
  expand_grid(time = times) %>%
  mutate(
    dose_index = vapply(dose, function(xx) which(xx == !!doses), integer(1)),
    time_index = vapply(time, function(xx) which(xx == !!times), integer(1)),
    response_doses = (c(1, 2.5, 3.5, 4, 4) / .90)[dose_index],
    response_times = c(0, .5, .9, 1, 0.95, 0.925, .90)[time_index],
    mean = 0.25 + response_doses * response_times,
    response = mean + rnorm(n(), 0, sd = 1.5)
  )

# Table 1
data_sumry <- distinct(data, dose, time, mean) %>%
  pivot_wider(names_from = "time", values_from = "mean", names_prefix = "time_")

mcmc <- dreamer_mcmc(
  data,
  mod_quad_itp = model_quad(
    mu_b1 = 0,
    sigma_b1 = 10,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = 0,
    sigma_b3 = 10,
    shape = 1,
    rate = .01,
    longitudinal = model_longitudinal_itp(
      mu_a = 0, sigma_a = 10, a_c1 = 0, b_c1 = 5, t_max = t_max
    ),
    w_prior = 1 / 8
  ),
  mod_logquad_itp = model_logquad(
    mu_b1 = 0,
    sigma_b1 = 10,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = 0,
    sigma_b3 = 10,
    shape = 1,
    rate = .01,
    longitudinal = model_longitudinal_itp(
      mu_a = 0, sigma_a = 10, a_c1 = 0, b_c1 = 5, t_max = t_max
    ),
    w_prior = 1 / 8
  ),
  mod_emax_itp = model_emax(
    mu_b1 = 0,
    sigma_b1 = 10,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = log(8),
    sigma_b3 = 1,
    mu_b4 = 1,
    sigma_b4 = 3,
    shape = 1,
    rate = .01,
    longitudinal = model_longitudinal_itp(
      mu_a = 0, sigma_a = 10, a_c1 = 0, b_c1 = 5, t_max = t_max
    ),
    w_prior = 1 / 8
  ),
  mod_exp_itp = model_exp(
    mu_b1 = 0,
    sigma_b1 = 10,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = -1,
    sigma_b3 = 5,
    shape = 1,
    rate = .01,
    longitudinal = model_longitudinal_itp(
      mu_a = 0, sigma_a = 10, a_c1 = 0, b_c1 = 5, t_max = t_max
    ),
    w_prior = 1 / 8
  ),
  mod_quad_idp = model_quad(
    mu_b1 = 0,
    sigma_b1 = 10,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = 0,
    sigma_b3 = 10,
    shape = 1,
    rate = .01,
    longitudinal = model_longitudinal_idp(
      mu_a = 0, sigma_a = 10, a_c1 = 0, b_c1 = 5, a_c2 = -5, b_c2 = 0, t_max = t_max
    ),
    w_prior = 1 / 8
  ),
  mod_logquad_idp = model_logquad(
    mu_b1 = 0,
    sigma_b1 = 10,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = 0,
    sigma_b3 = 10,
    shape = 1,
    rate = .01,
    longitudinal = model_longitudinal_idp(
      mu_a = 0, sigma_a = 10, a_c1 = 0, b_c1 = 5, a_c2 = -5, b_c2 = 0, t_max = t_max
    ),
    w_prior = 1 / 8
  ),
  mod_emax_idp = model_emax(
    mu_b1 = 0,
    sigma_b1 = 10,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = log(8),
    sigma_b3 = 1,
    mu_b4 = 1,
    sigma_b4 = 3,
    shape = 1,
    rate = .01,
    longitudinal = model_longitudinal_idp(
      mu_a = 0, sigma_a = 10, a_c1 = 0, b_c1 = 5, a_c2 = -5, b_c2 = 0, t_max = t_max
    ),
    w_prior = 1 / 8
  ),
  mod_exp_idp = model_exp(
    mu_b1 = 0,
    sigma_b1 = 10,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = -1,
    sigma_b3 = 5,
    shape = 1,
    rate = .01,
    longitudinal = model_longitudinal_idp(
      mu_a = 0, sigma_a = 10, a_c1 = 0, b_c1 = 5, a_c2 = -5, b_c2 = 0, t_max = t_max
    ),
    w_prior = 1 / 8
  )
)

mcmc$w_post

# Figure 4
p1 <- plot(mcmc$mod_quad_itp, data = data, times = 24) +
  ggtitle(sprintf(
    "Quadratic, ITP, Weight = %.0f%%",
    100 * mcmc$w_post["mod_quad_itp"]
  ))
p2 <- plot(mcmc$mod_logquad_itp, data = data, times = 24) +
  ggtitle(sprintf("Log-Quadratic, ITP, Weight = %.0f%%", 100 * mcmc$w_post["mod_logquad_itp"]))
p3 <- plot(mcmc$mod_emax_itp, data = data, times = 24) +
  ggtitle(sprintf("EMAX, ITP, Weight = %.0f%%", 100 * mcmc$w_post["mod_emax_itp"]))
p4 <- plot(mcmc$mod_exp_itp, data = data, times = 24) +
  ggtitle(sprintf("Exponential, ITP, Weight = %.0f%%", 100 * mcmc$w_post["mod_exp_itp"]))
p5 <- plot(mcmc$mod_quad_idp, data = data, times = 24) +
  ggtitle(sprintf("Quadratic, IDP, Weight = %.0f%%", 100 * mcmc$w_post["mod_quad_idp"]))
p6 <- plot(mcmc$mod_logquad_idp, data = data, times = 24) +
  ggtitle(sprintf("Log-Quadratic IDP, Weight = %.0f%%", 100 * mcmc$w_post["mod_logquad_idp"]))
p7 <- plot(mcmc$mod_emax_idp, data = data, times = 24) +
  ggtitle(sprintf("EMAX, IDP, Weight = %.0f%%", 100 * mcmc$w_post["mod_emax_idp"]))
p8 <- plot(mcmc$mod_exp_idp, data = data, times = 24) +
  ggtitle(sprintf("Exponential, IDP, Weight = %.0f%%", 100 * mcmc$w_post["mod_exp_idp"]))
p9 <- plot(mcmc, data = data, times = 24) +
  ggtitle("Bayesian Model Averaging")
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9

# Figure 5
p <- plot(mcmc, data = data)
p + facet_wrap(~dose, labeller = label_both) +
  ggtitle("") +
  scale_color_manual(values = rep("black", 5)) +
  theme(legend.position = "none")
