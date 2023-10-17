library(targets)
library(dplyr)

tar_load(int_widths)
tmp <- int_widths %>%
  dplyr::filter(
    scenario == "emax-itp",
    p_interim == .5,
    enroll_rate == 10,
    milestone == "final",
    time == max(time, na.rm = TRUE),
    dose == 15
  )

bma <- tmp %>%
  dplyr::filter(model == "bma", n_per_arm == 50) %>%
  pull(median_width)
truth <- tmp %>%
  dplyr::filter(model == "emax_long_itp", n_per_arm == 50) %>%
  pull(median_width)
1 - truth / bma

tar_load(mse)
tmp <- mse %>%
  ungroup() %>%
  dplyr::filter(
    scenario == "emax-itp",
    p_interim == .5,
    enroll_rate == 10,
    milestone == "final",
    time == max(time, na.rm = TRUE),
    dose == 15
  )

bma <- tmp %>%
  dplyr::filter(model == "bma", n_per_arm == 50) %>%
  pull(mse)
truth <- tmp %>%
  dplyr::filter(model == "emax_long_itp", n_per_arm == 50) %>%
  pull(mse)
1 - truth / bma
