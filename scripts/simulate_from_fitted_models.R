library(tidyverse)
library(here)
library(janitor)
library(blantyreESBL)
library(patchwork)
library(kableExtra)
library(tidytext)
library(tidybayes)


mod_sum_df <-
  read_rds(here("data_processed/gp_mod_sum_df.rda"))

mod_draws <-
  read_rds(here("data_processed/gp_drawslistout.rda"))

outcome_vars <-
  read_rds(here("data_processed/gp_outcome_vars.rda"))


simulate_exposure <- function(outcome_var, exposure_var, mod_draws, outcome_vars,
                              runin = 7, t_max = 14, b2_exp = 10, cov_exp = 7) {
  # outcome_var = string, outcome AMR subclass
  # exposure_var = name of covariate from model
  # exposure_string = human readable string of exposure for plots
  # mod_draws = list of draws frommodel fits
  # outcome_vars = list of outcome vars saved from model fit
  if (length(outcome_var) > 1 | length(exposure_var) > 1) {
    stop("more than one outcome or exposure passed to simulate_exposure")
  }

  if (!exposure_var %in% c("ceftriaxone", "cotrimoxazole", "ciprofloxacin", "amoxicillin")) {
    stop("exposure var not found")
  }
  if (!outcome_var %in% outcome_vars) {
    stop("outcome var not found")
  }

  current_outcome_pos <- which(outcome_vars == outcome_var)

  d <- mod_draws[[current_outcome_pos]]

  res_post <-
    d[, , dimnames(d)[[3]][grepl("beta|tau", dimnames(d)[[3]])]] |>
    as_tibble(rownames = "itr") |>
    pivot_longer(-itr) |>
    mutate(chain = str_extract(name, "[0-9]")) |>
    mutate(var = gsub("[1-4]\\.", "", name)) |>
    select(-name) |>
    pivot_wider(id_cols = c(itr, chain), names_from = var, values_from = value) |>
    janitor::clean_names()

  generate_coef <- function(cov, tau) {
    coef <-
      if (cov == 0) {
        return(0)
      } else if (cov > 0) {
        return(1)
      } else if (cov < 0) {
        return(exp(cov / tau))
      }
  }
  # runin = 7
  # t_max = 100
  # b1_exp <- 7
  # b2_exp <- 10

  if (exposure_var == "ceftriaxone") {
    cov_beta_1 <- c(rep(1, cov_exp), seq(0, -(t_max - cov_exp) + 1, -1))
  } else {
    cov_beta_1 <- 0
  }

  if (exposure_var == "cotrimoxazole") {
    cov_beta_3 <- c(rep(1, cov_exp), seq(0, -(t_max - cov_exp) + 1, -1))
  } else {
    cov_beta_3 <- 0
  }
  if (exposure_var == "ciprofloxacin") {
    cov_beta_4 <- c(rep(1, cov_exp), seq(0, -(t_max - cov_exp) + 1, -1))
  } else {
    cov_beta_4 <- 0
  }
  if (exposure_var == "amoxicillin") {
    cov_beta_5 <- c(rep(1, cov_exp), seq(0, -(t_max - cov_exp) + 1, -1))
  } else {
    cov_beta_5 <- 0
  }

  sim_df <-
    bind_rows(
      tibble(
        t = seq(-runin, 0, 1),
        cov_beta_1 = 0,
        cov_beta_2 = 0,
        cov_beta_3 = 0,
        cov_beta_4 = 0,
        cov_beta_5 = 0
      ),
      tibble(
        t = 1:t_max
      ) |>
        mutate(
          cov_beta_1 = cov_beta_1,
          cov_beta_2 = c(rep(1, b2_exp), seq(0, -(t_max - b2_exp) + 1, -1)),
          cov_beta_3 = cov_beta_3,
          cov_beta_4 = cov_beta_4,
          cov_beta_5 = cov_beta_5
        )
    ) |>
    cross_join(res_post) |>
    rowwise() |>
    mutate(
      pr = plogis(
        beta_0 +
          beta_1 * generate_coef(cov_beta_1, tau) +
          beta_2 * generate_coef(cov_beta_2, tau) +
          beta_3 * generate_coef(cov_beta_3, tau) +
          beta_4 * generate_coef(cov_beta_4, tau) +
          beta_5 * generate_coef(cov_beta_5, tau)
      ),
      exposure = exposure_var,
      outcome = outcome_var
    )

  return(sim_df)
}

simulate_exposure_multi <- function(outcome_var_to_sim, exposure_vars, mod_draws, outcome_vars,
                                    runin = 7, t_max = 14, b2_exp = 10, cov_exp = 7) {
  i <- 1
  listout <- list()
  for (outcome_var in outcome_var_to_sim) {
    for (exposure_var in exposure_vars) {
      listout[[i]] <- simulate_exposure(
        exposure_var = exposure_var,
        outcome_var = outcome_var,
        mod_draws = mod_draws,
        outcome_vars = outcome_vars,
        runin = runin,
        t_max = t_max,
        b2_exp = b2_exp,
        cov_exp = 7
      )
      i <- i + 1
    }
  }

  return(bind_rows(listout))
}

df <-
  simulate_exposure_multi(
    exposure_var = c("ceftriaxone", "amoxicillin", "cotrimoxazole", "ciprofloxacin"),
    outcome_var = c("aminoglycoside", "cephalosporin", "clindamycin_erythromycin_streptogramin b", "quinolone", "carbapenem"),
    mod_draws = mod_draws,
    outcome_vars = outcome_vars
  )


df |>
  group_by(t, exposure, outcome) |>
  summarise(
    med = median(pr),
    lci = quantile(pr, 0.025),
    uci = quantile(pr, 0.975)
  ) |>
  ggplot(aes(t, med, ymin = lci, ymax = uci, color = exposure, fill = exposure)) +
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.3) +
  xlim(c(-7, 20)) +
  facet_grid(exposure ~ outcome)


df |>
  filter(
    exposure %in% c("amoxicillin", "ceftriaxone", "cotrimoxazole"),
    outcome %in% c("cephalosporin", "aminoglycoside", "clindamycin_erythromycin_streptogramin b")
  ) |>
  group_by(t, exposure, outcome) |>
  mutate(
    exposure = factor(exposure, levels = c("ceftriaxone", "amoxicillin", "cotrimoxazole")),
    outcome = factor(outcome, c("cephalosporin", "aminoglycoside", "clindamycin_erythromycin_streptogramin b"))
  ) |>
  summarise(
    med = median(pr),
    lci = quantile(pr, 0.025),
    uci = quantile(pr, 0.975)
  ) |>
  ggplot(aes(t, med,
    # ymin = lci, ymax = uci,
    linetype = exposure
  )) +
  geom_line() +
  # geom_ribbon(color = NA, alpha = 0.3) +
  xlim(c(-4, 15)) +
  facet_grid(~outcome) +
  theme_bw() +
  labs(y = "Proportion", linetype = "Exposure")

ggsave(here("eccmid_plot1.pdf"), width = 10, height = 3)


ggsave(here("eccmid_plot1.png"), width = 10, height = 3)

# taxonomy


tax_mod_sum_df <-
  read_rds(here("data_processed/gp_taxonomy_mod_sum_df.rda"))

tax_mod_draws <-
  read_rds(here("data_processed/gp_taxonomy_drawslistout.rda"))

tax_outcome_vars <-
  read_rds(here("data_processed/gp_taxonomy_outcome_vars.rda"))



tax_simulate_exposure <- function(outcome_var, exposure_var, mod_draws, outcome_vars,
                                  runin = 7, t_max = 14, b2_exp = 10, cov_exp = 7) {
  # outcome_var = string, outcome AMR subclass
  # exposure_var = name of covariate from model
  # exposure_string = human readable string of exposure for plots
  # mod_draws = list of draws frommodel fits
  # outcome_vars = list of outcome vars saved from model fit
  if (length(outcome_var) > 1 | length(exposure_var) > 1) {
    stop("more than one outcome or exposure passed to simulate_exposure")
  }

  if (!exposure_var %in% c("ceftriaxone", "cotrimoxazole", "ciprofloxacin", "amoxicillin")) {
    stop("exposure var not found")
  }
  if (!outcome_var %in% outcome_vars) {
    stop("outcome var not found")
  }

  current_outcome_pos <- which(outcome_vars == outcome_var)

  d <- mod_draws[[current_outcome_pos]]

  res_post <-
    d[, , dimnames(d)[[3]][grepl("beta|tau|phi", dimnames(d)[[3]])]] |>
    as_tibble(rownames = "itr") |>
    pivot_longer(-itr) |>
    mutate(chain = str_extract(name, "[0-9]")) |>
    mutate(var = gsub("[1-4]\\.", "", name)) |>
    select(-name) |>
    pivot_wider(id_cols = c(itr, chain), names_from = var, values_from = value) |>
    janitor::clean_names()

  generate_coef <- function(cov, tau) {
    coef <-
      if (cov == 0) {
        return(0)
      } else if (cov > 0) {
        return(1)
      } else if (cov < 0) {
        return(exp(cov / tau))
      }
  }
  # runin = 7
  # t_max = 100
  # b1_exp <- 7
  # b2_exp <- 10

  if (exposure_var == "ceftriaxone") {
    cov_beta_1 <- c(rep(1, cov_exp), seq(0, -(t_max - cov_exp) + 1, -1))
  } else {
    cov_beta_1 <- 0
  }

  if (exposure_var == "cotrimoxazole") {
    cov_beta_3 <- c(rep(1, cov_exp), seq(0, -(t_max - cov_exp) + 1, -1))
  } else {
    cov_beta_3 <- 0
  }
  if (exposure_var == "ciprofloxacin") {
    cov_beta_4 <- c(rep(1, cov_exp), seq(0, -(t_max - cov_exp) + 1, -1))
  } else {
    cov_beta_4 <- 0
  }
  if (exposure_var == "amoxicillin") {
    cov_beta_5 <- c(rep(1, cov_exp), seq(0, -(t_max - cov_exp) + 1, -1))
  } else {
    cov_beta_5 <- 0
  }

  sim_df <-
    bind_rows(
      tibble(
        t = seq(-runin, 0, 1),
        cov_beta_1 = 0,
        cov_beta_2 = 0,
        cov_beta_3 = 0,
        cov_beta_4 = 0,
        cov_beta_5 = 0
      ),
      tibble(
        t = 1:t_max
      ) |>
        mutate(
          cov_beta_1 = cov_beta_1,
          cov_beta_2 = c(rep(1, b2_exp), seq(0, -(t_max - b2_exp) + 1, -1)),
          cov_beta_3 = cov_beta_3,
          cov_beta_4 = cov_beta_4,
          cov_beta_5 = cov_beta_5
        )
    ) |>
    cross_join(res_post) |>
    rowwise() |>
    mutate(
      reads =
        exp(
          beta_0 +
            beta_1 * generate_coef(cov_beta_1, tau) +
            beta_2 * generate_coef(cov_beta_2, tau) +
            beta_3 * generate_coef(cov_beta_3, tau) +
            beta_4 * generate_coef(cov_beta_4, tau) +
            beta_5 * generate_coef(cov_beta_5, tau)
          # log(2.5e7)
        ),
      exposure = exposure_var,
      outcome = outcome_var
    )

  return(sim_df)
}


tax_simulate_exposure_multi <- function(outcome_var_to_sim, exposure_vars, mod_draws, outcome_vars,
                                        runin = 7, t_max = 14, b2_exp = 10, cov_exp = 7) {
  i <- 1
  listout <- list()
  for (outcome_var in outcome_var_to_sim) {
    for (exposure_var in exposure_vars) {
      listout[[i]] <- tax_simulate_exposure(
        exposure_var = exposure_var,
        outcome_var = outcome_var,
        mod_draws = mod_draws,
        outcome_vars = outcome_vars,
        runin = runin,
        t_max = t_max,
        b2_exp = b2_exp,
        cov_exp = cov_exp
      )
      i <- i + 1
    }
  }

  return(bind_rows(listout))
}


df2 <-
  tax_simulate_exposure_multi(
    exposure_var = c("ceftriaxone", "amoxicillin", "cotrimoxazole", "ciprofloxacin"),
    outcome_var = c("phylum_Proteobacteria", "order_Enterobacterales", "genus_Escherichia"),
    mod_draws = tax_mod_draws,
    outcome_vars = tax_outcome_vars
  )

# df2 |>
#   group_by(t, exposure, outcome) |>
#   summarise(
#     med = median(reads),
#     lci = quantile(reads, 0.025),
#     uci = quantile(reads, 0.975)
#   ) |>
#   ggplot(aes(t, med, ymin = lci, ymax = uci, color = exposure, fill = exposure)) +
#   geom_line() +
#   geom_ribbon(color = NA, alpha = 0.3) +
#   xlim(c(-7, 20)) +
#   facet_grid(exposure ~ outcome)

df2 |>
  filter(exposure %in% c("amoxicillin", "ceftriaxone", "cotrimoxazole")) |>
  group_by(t, exposure, outcome) |>
  mutate(exposure = factor(exposure, levels = c("ceftriaxone", "amoxicillin", "cotrimoxazole"))) |>
  summarise(
    med = median(reads),
    lci = quantile(reads, 0.025),
    uci = quantile(reads, 0.975)
  ) |>
  ggplot(aes(t, med,
    # ymin = lci, ymax = uci,
    linetype = exposure
  )) +
  geom_line() +
  # geom_ribbon(color = NA, alpha = 0.3) +
  xlim(c(-4, 15)) +
  facet_grid(~outcome, scales = "free") +
  theme_bw() +
  labs(y = "Proportion", linetype = "Exposure")
  

ggsave(here("eccmid_plot2.pdf"), width = 10, height = 3)
ggsave(here("eccmid_plot2.png"), width = 10, height = 3)
