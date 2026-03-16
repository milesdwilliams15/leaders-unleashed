# -------------------------------------------------------------------------
# Functions with extra goodies for {peacesciencer}
# -------------------------------------------------------------------------

library(peacesciencer)



# x_spell() ---------------------------------------------------------------

x_spell <- function(x) {
  nx <- 1 - x
  for(i in 1:length(x)) {
    if(i == 1) {
      nx[i] <- nx[i]
    } else {
      nx[i] <- (1 + nx[i-1]) * (1 - x[i]) +
        x[i] * (1 - x[i-1]) * (1 + nx[i - 1])
      if(x[i] == 0 & x[i-1] == 1) nx[i] <- 1
    }
  }
  ifelse(nx > 0, nx - 1, nx)
}
# add_icd_mics() ----------------------------------------------------------

add_icd_mics <- function(
    data, level
) {
  # check for CoW system codes:
  yes <- any(colnames(data) %in% c("ccode", "ccode1"))
  if(!yes) {
    stop("You can only use this function with CoW system codes.")
  }
  
  # whittle MIE data to MIC level:
  dt <- suppressMessages(
    read_csv(
      "https://raw.githubusercontent.com/milesdwilliams15/foreign-figures/refs/heads/main/_data/mie-1.0.csv"
    )
  )
  
  dt |>
    group_by(
      ccode1, ccode2, styear, micnum
    ) |>
    summarize(
      init = max((eventnum == 1), na.rm = T),
      across(fatalmin1:fatalmax2, sum),
      hostlev = max(hostlev),
      .groups = "drop"
    ) |> 
    group_by(micnum) |>
    mutate(
      miconset = (min(styear) == styear) + 0,
      hostlev = max(hostlev)
    ) -> dt
  
  # only keep those of certain hostlev:
  if(!missing(level)) {
    dt |>
      filter(hostlev >= level) -> dt
  }
  
  # some more aggregating:
  dt |>
    group_by(ccode1, ccode2, styear, micnum) |>
    summarize(
      micongoing = 1,
      miconset = max(miconset),
      micongoing_init1 = max(init),
      miconset_init1 = max(init) * miconset,
      across(fatalmin1:fatalmax2, sum),
      .groups = "drop"
    ) -> dt
  
  # the data isn't mirrored:
  bind_rows(
    dt,
    dt |>
      rename(
        ccode2 = ccode1,
        ccode1 = ccode2,
        micongoing_init2 = micongoing_init1,
        miconset_init2 = miconset_init1,
        fatalmin2 = fatalmin1,
        fatalmin1 = fatalmin2,
        fatalmax2 = fatalmax1,
        fatalmax1 = fatalmax2
      )
  ) |> 
    rename(year = styear) |>
    mutate(
      across(everything(), ~ replace_na(.x, 0))
    ) |> 
    select(
      ccode1, 
      ccode2, 
      year, 
      everything()
    ) -> dt
  
  # if the data is dyadic:
  if(any(colnames(data) %in% c("ccode1"))) {
    data |>
      left_join(
        dt, by = c("ccode1", "ccode2", "year")
      ) -> dt
    
    # fill in missings:
    dt |>
      mutate(
        across(
          starts_with("mic") | starts_with("fatal"),
          ~ replace_na(.x, 0)
        )
      ) -> dt
    
    # compute peace spell
    dt |>
      mutate(id = paste0(ccode1, "-", ccode2)) |>
      stevemisc::ps_spells(
        micongoing,
        year,
        id,
        ongoing = FALSE
      ) |>
      rename(
        micspell = spell
      ) |>
      select(-id) -> dt
    
    # preserve attributes
    attributes(dt)$ps_data_type <- "dyad_year"
    attributes(dt)$ps_system <- "cow"
  } else {
    dt |>
      rename(ccode = ccode1) |>
      group_by(ccode, year) |>
      summarize(
        micongoing = 1,
        miconset = max(miconset),
        micongoing_init = max(micongoing_init1),
        miconset_init = max(miconset_init1),
        fatalmin = sum(fatalmin1),
        fatalmax = sum(fatalmax1),
        fatalmin_total = sum(fatalmin1 + fatalmin2),
        fatalmax_total = sum(fatalmax1 + fatalmax2),
        .groups = "drop"
      ) -> dt
    
    # merge:
    data |>
      left_join(
        dt,
        by = c("ccode", "year")
      ) -> dt
    
    # fill missings:
    dt |>
      mutate(
        across(
          starts_with("mic") | starts_with("fatal"),
          ~ replace_na(.x, 0)
        )
      ) -> dt
    
    # compute peace spell
    dt |>
      stevemisc::ps_spells(
        micongoing,
        year,
        ccode,
        ongoing = FALSE
      ) |>
      rename(
        micspell = spell
      ) -> dt
    
    ## preserve attributes
    attributes(dt)$ps_data_type <- "state_year"
    attributes(dt)$ps_system <- "cow"
  }
  
  # return
  dt
}


# add_opportunity() ---------------------------------------------------------

get_opportunity_data <- function(data) {
  create_dyadyears(subset_years = 1816:2014) |>
    add_cow_majors() |>
    add_contiguity() |>
    add_capital_distance() -> ddy
  
  logit <- function(x) 1 / (1 + exp(-x))
  
  ddy |>
    mutate(
      contig = ifelse(conttype >= 1, 1, 0),
      majdyad = pmax(cowmaj1, cowmaj2)
    ) |>
    transmute(
      ccode1, ccode2, year,
      contdyads = contig,
      prd = ifelse(contig == 1 | majdyad == 1, 1, 0),
      bcd = logit(
        4.801 + 4.50*contig - 1.051*log(capdist) + 2.901*majdyad)
    ) -> ddy
  
  ddy |>
    group_by(year) |>
    summarize(
      dyads = n(),
      contdyads = sum(contdyads),
      prd = sum(prd),
      bcd = sum(bcd)
    )
}


# glm_robust() ------------------------------------------------------------

glm_robust <- function(formula, family, data, se_type = "HC0", clusters = NULL) {
  ## the fit
  if(missing(family)) {
    fit <- glm(formula, family = binomial, data)
  } else {
    fit <- glm(formula, family, data)
  }
  
  ## the cov matrix
  if(is.null(clusters)) {
    vcv <- sandwich::vcovCL(fit, type = se_type)
  } else {
    vcv <- sandwich::vcovCL(fit, type = se_type, cluster = data[, clusters])
  }
  
  ## prep the summary
  ctest <- lmtest::coeftest(fit, vcv) |> broom::tidy()
  cis   <- lmtest::coefci(fit, vcov. = vcv) 
  ctest <- ctest |>
    mutate(lower = cis[, 1], upper = cis[, 2])
  
  ## add some extra attributes
  attributes(ctest)$fit <- fit
  attributes(ctest)$vcov <- vcv
  
  ## return
  ctest
}


# gam_robust() ------------------------------------------------------------

gam_robust <- function(formula, family, data, se_type = "HC0", clusters = NULL) {
  ## the fit
  if(missing(family)) {
    fit <- gam(formula, family = binomial, data)
  } else {
    fit <- gam(formula, family, data)
  }
  
  ## the cov matrix
  if(is.null(clusters)) {
    vcv <- sandwich::vcovCL(fit, type = se_type)
  } else {
    vcv <- sandwich::vcovCL(fit, type = se_type, cluster = data[, clusters])
  }
  
  ## prep the summary
  ctest <- lmtest::coeftest(fit, vcv) |> broom::tidy()
  cis   <- lmtest::coefci(fit, vcov. = vcv) 
  ctest <- ctest |>
    mutate(lower = cis[, 1], upper = cis[, 2])
  
  ## add some extra attributes
  attributes(ctest)$fit <- fit
  attributes(ctest)$vcov <- vcv
  
  ## return
  ctest
}




# coef_plot() -------------------------------------------------------------

coef_plot <- function(tidy_fit, coef_map) {
  if(missing(coef_map)) {
    coef_map <- unique(tidy_fit$term)
    names(coef_map) <- coef_map
  } else {
    names(coef_map) <- unique(tidy_fit$term)
  }
  if(with(tidy_fit, !exists("model"))) {
    tidy_fit$model <- NA
  }
  tidy_fit |>
    mutate(
      ord = n():1
    ) -> tidy_fit
  tidy_fit |>
    ggplot() +
    aes(estimate, reorder(term, ord), xmin = lower, xmax = upper, color = model) +
    geom_pointrange(
      size = .1,
      position = position_dodge(
        ifelse(any(is.na(tidy_fit$model)), 0, -.5)
      )
    ) +
    geom_vline(
      xintercept = 0,
      lty = 2
    ) +
    scale_y_discrete(
      labels = coef_map
    ) +
    labs(
      x = "Estimate",
      y = NULL,
      color = NULL
    ) +
    theme(
      plot.title.position = "plot",
      legend.position = ifelse(
        any(is.na(tidy_fit$model)),
        "", "right"
      )
    )
}


# sim_pred() --------------------------------------------------------------

sim_pred <- function(model, newdata, its = 1000) {
  if(missing(newdata)) stop("Missing 'newdata'.")
  fit <- attributes(model)$fit
  vcv <- attributes(model)$vcov
  L   <- function(x) 1 / (1 + exp(-x))
  sim_preds <- replicate(
    n = its,
    expr = {
      nfit <- fit
      nfit$coefficients <- MASS::mvrnorm(
        n = length(fit$coefficients),
        mu = fit$coefficients,
        Sigma = vcv
      ) |> diag()
      ndata <- newdata
      ndata$p <- predict(
        nfit, newdata = newdata
      ) |> L()
      list(ndata)
    }
  )
  out <- dplyr::bind_rows(sim_preds)
  
  ## return
  out
}