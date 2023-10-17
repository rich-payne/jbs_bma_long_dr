scenarios <- expand.grid(
  n_per_arm = 50,
  scenario = c(
    "null-null", "plateau-plateau", "plateau-decrease", "emax-itp", "emax-idp"
  ),
  p_interim = c(.5, .75),
  enroll_rate = c(10, 30),
  stringsAsFactors = FALSE
)

scenarios_binary <- expand.grid(
  n_per_arm = 50,
  scenario = c(
    "null-null", "plateau-plateau", "plateau-decrease", "emax-itp", "emax-idp"
  ),
  p_interim = .5,
  enroll_rate = 10,
  stringsAsFactors = FALSE
)
