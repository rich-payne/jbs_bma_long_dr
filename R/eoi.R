get_eois_impl <- function(x, doses, model, eois) {
  eoi_grid <- expand.grid(
    eoi = eois,
    dose = doses,
    reference_dose = 0
  )
  out <- pmap_df(eoi_grid, pr_eoi, x = x)
  if (rlang::has_name(out, "time")) {
    out <- out %>%
      pivot_wider(
        values_from = "prob",
        names_from = c("eoi", "time", "dose", "reference_dose"),
        names_glue = paste0("pr_eoi_{eoi}_arm1_{dose}_arm2_{reference_dose}_time_{time}_", model)
      )
  } else {
    out <- out %>%
      pivot_wider(
        values_from = "prob",
        names_from = c("eoi", "dose", "reference_dose"),
        names_glue = paste0("pr_eoi_{eoi}_arm1_{dose}_arm2_{reference_dose}_time_none_", model)
      )
  }
}

get_eois <- function(x, doses, eois = c(2.5, 2.75)) {
  bma_eois <- get_eois_impl(x, doses, model = "bma", eois = eois)
  ind <- mcmc_indices(x)
  out <- map2(
    x[ind],
    names(x[ind]),
    get_eois_impl,
    doses = doses,
    eois = eois
  ) %>%
    bind_cols(bma_eois)
  return(out)
}
