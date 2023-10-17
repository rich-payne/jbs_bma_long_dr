library(dplyr)
library(tidyr)
library(gt)
library(dreamer)

source("R/data_generation.R")

data <- bind_rows(
  gen_data(n_per_arm = 1, scenario = "null-null", enroll_rate = 1, sigma = 0) %>%
    mutate(scenario = "Null-Null"),
  gen_data(n_per_arm = 1, scenario = "plateau-plateau", enroll_rate = 1, sigma = 0) %>%
    mutate(scenario = "Plateau-Plateau"),
  gen_data(n_per_arm = 1, scenario = "plateau-decrease", enroll_rate = 1, sigma = 0) %>%
    mutate(scenario = "Plateau-Decrease"),
  gen_data(n_per_arm = 1, scenario = "emax-itp", enroll_rate = 1, sigma = 0) %>%
    mutate(scenario = "EMAX-ITP"),
  gen_data(n_per_arm = 1, scenario = "emax-idp", enroll_rate = 1, sigma = 0) %>%
    mutate(scenario = "EMAX-IDP"),
) %>%
  select(Scenario = scenario, Week = time, Dose = dose, Mean = response)

data_table <- data %>%
  pivot_wider(names_from = "Week", values_from = "Mean", names_prefix = "Week ") %>%
  arrange(
    match(
      Scenario,
      c(
        "Null-Null", "Plateau-Plateau",
        "Plateau-Decrease", "EMAX-ITP", "EMAX-IDP"
      )
    ),
    Dose
  ) %>%
  filter(!(Scenario %in% c("EMAX-ITP", "EMAX-IDP")))

# table 5
tab <- gt(data_table, rowname_col = "Scenario") %>%
  tab_stubhead("Scenario") %>%
  fmt_number(starts_with("Week"))
tab
# as_latex(tab)
