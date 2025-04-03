
library(tidyverse)
library(phyloseq)
library(here)
library(blantyreESBL)
library(microViz)
library(janitor)
library(cmdstanr)
library(bayesplot)
library(patchwork)
library(pheatmap)
library(ggpubr)
library(viridis)
library(kableExtra)
library(tidytext)
library(scales)
library(ggsci)
library(ggplotify)
library(blantyreSepsis)

samples_to_drop <- read_lines(here("data_processed/samples_to_drop.txt"))

df <-
  read_tsv(
    here(
      "data_raw/ECO_AMRfinder_all_metagenomes_default_17_11_23/AMR_binary_matrix.tsv"
    )
  ) |>
  janitor::clean_names()

df <-
  df |>
  pivot_longer(-c(class, subclass, gene_symbol)) |>
  mutate(
    name =
      str_extract(
        name,
        "sample_[0-9]{1,3}_(.*)_(:?uncorrected_amr|amr)",
        group = 1
      )
  ) |>
  mutate(name = gsub("_uncorrected", "", name)) |>
  pivot_wider(
    id_cols = c(class, subclass, gene_symbol),
  )

df2 <-
  df |>
  mutate(gene_symbol = make_clean_names(gene_symbol)) |>
  select(-c(class, subclass)) |>
  pivot_longer(-gene_symbol) |>
  filter(!grepl("^nc", name)) |>
  pivot_wider(id_cols = name, names_from = gene_symbol, values_from = value) |>
  mutate(name = gsub("_.*$", "", name)) |>
  left_join(
    btESBL_stoolESBL |>
      select(pid, lab_id, arm, visit, ESBL, enroll_date, data_date) |>
      mutate(lab_id = tolower(lab_id)),
    by = c("name" = "lab_id")
  )

df2 <-
  df2 |>
  filter(!name %in% tolower(samples_to_drop))

df <-
  df |>
  select(!matches(paste(tolower(samples_to_drop), collapse = "|")))
included_participants <-
  df2 |>
  pull(pid) |>
  unique()



d <- read_rds(here("data_processed/diversity_mod_draws.rda"))

# simulations

sim_diversity <- function(mod_draws, exposure_var = "ceftriaxone",
                          runin = 7, t_max = 14, b2_exp = 10, cov_exp = 7,
                          t_scale_factor = 55.25) {
  # mod_draws = list of draws frommodel fits

  d <- mod_draws

  res_post <-
    d[, , dimnames(d)[[3]][grepl("beta|tau|sigma", dimnames(d)[[3]])]] |>
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
        return(exp(cov / (tau * t_scale_factor)))
      }
  }
  # runin = 7
  # t_max = 100
  # b1_exp <- 7
  # b2_exp <- 10

  if (exposure_var == "ceftriaxone") {
    cov_beta_1 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
  } else {
    cov_beta_1 <- 0
  }
  if (exposure_var == "cotrimoxazole") {
    cov_beta_3 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
  } else {
    cov_beta_3 <- 0
  }
  if (exposure_var == "ciprofloxacin") {
    cov_beta_4 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
  } else {
    cov_beta_4 <- 0
  }
  if (exposure_var == "amoxicillin") {
    cov_beta_5 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
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
          cov_beta_2 = c(rep(1, b2_exp), seq(-1, -(t_max - b2_exp), -1)),
          cov_beta_3 = cov_beta_3,
          cov_beta_4 = cov_beta_4,
          cov_beta_5 = cov_beta_5
        )
    )

  if (b2_exp == 0) {
    sim_df$cov_beta_2 <- 0
  }

  sim_df <-
    sim_df |>
    cross_join(res_post) |>
    rowwise() |>
    mutate(
      pr =
        beta_0 +
          beta_1 * generate_coef(cov_beta_1, tau) +
          beta_2 * generate_coef(cov_beta_2, tau) +
          beta_3 * generate_coef(cov_beta_3, tau) +
          beta_4 * generate_coef(cov_beta_4, tau) +
          beta_5 * generate_coef(cov_beta_5, tau),
      exposure = exposure_var,
    )

  return(sim_df)
}

sim_df <-
  bind_rows(
    sim_diversity(d, exposure = "ceftriaxone", t_max = 100),
    sim_diversity(d, exposure = "ciprofloxacin", t_max = 100),
    sim_diversity(d, exposure = "cotrimoxazole", t_max = 100),
    sim_diversity(d, exposure = "amoxcicillin", t_max = 100),
    sim_diversity(d, exposure = "no antibiotic", t_max = 100),
    sim_diversity(d, exposure = "community", b2_exp = 0, t_max = 100)
  )

p_1 <-
  sim_df |>
  group_by(t, exposure) |>
  summarise(
    diversity = median(pr),
    lci = quantile(pr, 0.025),
    uci = quantile(pr, 0.975)
  ) |>
  mutate(
    exposure = case_when(
      grepl("community", exposure) ~ "Community\nNo Antibiotic",
      grepl("antibiotic", exposure) ~ "Hospitalised\nNo Antibiotic",
      grepl("amox", exposure) ~ "Hospitalised\nAmoxicillin",
      grepl("cotrim", exposure) ~ "Hospitalised\nCotrimoxazole",
      grepl("cipro", exposure) ~ "Hospitalised\nCiprofloxacin",
      grepl("ceft", exposure) ~ "Hospitalised\nCeftriaxone"
    ),
    exposure =
      factor(exposure,
        levels = c(
          "Community\nNo Antibiotic",
          "Hospitalised\nNo Antibiotic",
          "Hospitalised\nAmoxicillin",
          "Hospitalised\nCotrimoxazole",
          "Hospitalised\nCiprofloxacin",
          "Hospitalised\nCeftriaxone"
        )
      )
  ) |>
  ggplot(aes(t, diversity, ymin = lci, ymax = uci)) +
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.3) +
  facet_wrap(~exposure, nrow = 1) +
  xlim(c(-1, 50)) +
  theme_bw() +
  labs(x = "Time (days)", y = "Mean Shannon Diversity", title = "Shannon Diversity vs Time", subtitle = "Effect of antmicrobial varies by agent")


#### ------------------


tax_mod_sum_df <-
  read_rds(here("data_processed/gp_taxonomy_mod_sum_df.rda"))

tax_mod_draws <-
  read_rds(here("data_processed/gp_taxonomy_drawslistout.rda"))

tax_outcome_vars <-
  read_rds(here("data_processed/gp_taxonomy_outcome_vars.rda"))



tax_simulate_exposure <- function(outcome_var, exposure_var, mod_draws, outcome_vars,
                                  runin = 7, t_max = 14, b2_exp = 10, cov_exp = 7,
                                  t_scale_factor = 55.25) {
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
        return(exp(cov / (tau * t_scale_factor)))
      }
  }
  # runin = 7
  # t_max = 100
  # b1_exp <- 7
  # b2_exp <- 10

  if (exposure_var == "ceftriaxone") {
    cov_beta_1 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
  } else {
    cov_beta_1 <- 0
  }

  if (b2_exp == 0) {
    cov_beta_2 <- 0
  } else {
    cov_beta_2 <- c(rep(1, b2_exp), seq(-1, -(t_max - b2_exp), -1))
  }

  if (exposure_var == "cotrimoxazole") {
    cov_beta_3 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
  } else {
    cov_beta_3 <- 0
  }
  if (exposure_var == "ciprofloxacin") {
    cov_beta_4 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
  } else {
    cov_beta_4 <- 0
  }
  if (exposure_var == "amoxicillin") {
    cov_beta_5 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
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
          cov_beta_2 = cov_beta_2,
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
                                        runin = 7, t_max = 30, b2_exp = 10, cov_exp = 7,
                                        t_scale_factor = 55.25) {
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
        cov_exp = cov_exp,
        t_scale_factor = t_scale_factor
      )
      i <- i + 1
    }
  }

  return(bind_rows(listout))
}


df_taxsim <-
  tax_simulate_exposure_multi(
    exposure_var = c("ceftriaxone", "amoxicillin", "cotrimoxazole", "ciprofloxacin"),
    outcome_var = c("phylum_Proteobacteria", "order_Enterobacterales", "genus_Escherichia"),
    mod_draws = tax_mod_draws,
    outcome_vars = tax_outcome_vars,
  )

p_2a <-
  df_taxsim |>
  filter(exposure %in% c("ceftriaxone", "cotrimoxazole")) |>
  group_by(t, exposure, outcome) |>
  mutate(exposure = factor(exposure, levels = c("ceftriaxone", "amoxicillin", "cotrimoxazole"))) |>
  summarise(
    med = median(reads),
    lci = quantile(reads, 0.25),
    uci = quantile(reads, 0.75)
  ) |>
  mutate(outcome = str_to_title(gsub("_", "\\\n", outcome))) |>
  ggplot(aes(t, med,
    ymin = lci, ymax = uci,
    linetype = exposure, color = exposure, fill = exposure
  )) +
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.3) +
  xlim(c(-4, 30)) +
  facet_wrap(~outcome, scales = "free") +
  theme_bw() +
  labs(
    y = "Mean proportion of reads", linetype = "Exposure", color = "Exposure", fill =
      "Exposure",
      title = "Microbiome composition vs time", subtitle = "Ceftriaxone vs Co-trimoxazole"
  ) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = viridis_pal(option = "A")(6)[c(2, 4)]) +
  scale_fill_manual(values = viridis_pal(option = "A")(6)[c(2, 4)])

p_2b <-
  df_taxsim |>
  filter(exposure %in% c("ciprofloxacin", "amoxicillin")) |>
  group_by(t, exposure, outcome) |>
  mutate(exposure = factor(exposure, levels = c("ceftriaxone", "ciprofloxacin", "amoxicillin", "cotrimoxazole"))) |>
  summarise(
    med = median(reads),
    lci = quantile(reads, 0.25),
    uci = quantile(reads, 0.75)
  ) |>
  mutate(outcome = str_to_title(gsub("_", "\\\n", outcome))) |>
  ggplot(aes(t, med,
    ymin = lci, ymax = uci,
    linetype = exposure, color = exposure, fill = exposure
  )) +
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.3) +
  xlim(c(-4, 30)) +
  facet_wrap(~outcome, scales = "free") +
  theme_bw() +
  labs(
    y = "Mean proportion of reads", linetype = "Exposure", color = "Exposure", fill =
      "Exposure", 
      title = "Microbiome Composition vs Time", subtitle = "Ciprofloxacin vs Amoxicillin"
  ) +
  theme(legend.position = "bottom") +
  scale_color_manual(values = viridis_pal(option = "B")(6)[c(3, 5)]) +
  scale_fill_manual(values = viridis_pal(option = "B")(6)[c(3, 5)])

## ---------------------------------------


mod_sum_df <-
  read_rds(here("data_processed/gp_mod_sum_df.rda"))

mod_draws <-
  read_rds(here("data_processed/gp_drawslistout.rda"))

outcome_vars <-
  read_rds(here("data_processed/gp_outcome_vars.rda"))


mod_sum_df <-
  read_rds(here("data_processed/gp_mod_sum_df.rda"))

fact_levels <-
  mod_sum_df |>
  filter(parameter == "beta[1]") |>
  mutate(outcome_var = str_to_title(gsub("_", " ", outcome_var))) |>
  arrange(desc(m)) |>
  pull(outcome_var)

# simulations



simulate_exposure <- function(outcome_var, exposure_var, mod_draws, outcome_vars,
                              runin = 7, t_max = 14, b2_exp = 10, cov_exp = 7,
                              t_scale_factor = 55.25) {
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
        return(exp(cov / (tau * t_scale_factor)))
      }
  }
  # runin = 7
  # t_max = 100
  # b1_exp <- 7
  # b2_exp <- 10

  if (exposure_var == "ceftriaxone") {
    cov_beta_1 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
  } else {
    cov_beta_1 <- 0
  }

  if (exposure_var == "cotrimoxazole") {
    cov_beta_3 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
  } else {
    cov_beta_3 <- 0
  }
  if (exposure_var == "ciprofloxacin") {
    cov_beta_4 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
  } else {
    cov_beta_4 <- 0
  }
  if (exposure_var == "amoxicillin") {
    cov_beta_5 <- c(rep(1, cov_exp), seq(-1, -(t_max - cov_exp), -1))
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
          cov_beta_2 = c(rep(1, b2_exp), seq(-1, -(t_max - b2_exp), -1)),
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
                                    runin = 7, t_max = 30, b2_exp = 10, cov_exp = 7, t_scale_factor = 55.25) {
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
        cov_exp = 7,
        t_scale_factor = t_scale_factor
      )
      i <- i + 1
    }
  }

  return(bind_rows(listout))
}


df_ressim <-
  simulate_exposure_multi(
    exposure_var = c("ceftriaxone", "amoxicillin", "cotrimoxazole", "ciprofloxacin"),
    outcome_var = c("aminoglycoside", "cephalosporin", "clindamycin_erythromycin_streptogramin b", "quinolone", "azithromycin_erythromycin_spiramycin_telithromycin"),
    mod_draws = mod_draws,
    outcome_vars = outcome_vars
  )

p_3a <-
  df_ressim |>
  filter(
    exposure %in% c("ceftriaxone", "cotrimoxazole"),
    outcome %in% c("cephalosporin", "aminoglycoside", "clindamycin_erythromycin_streptogramin b")
  ) |>
  group_by(t, exposure, outcome) |>
  mutate(
    outcome = str_to_title(gsub("_", " ", outcome)),
    exposure = str_to_title(exposure),
    exposure = factor(exposure, levels = c("Ceftriaxone", "Amoxicillin", "Cotrimoxazole")),
    outcome = factor(outcome, c("Cephalosporin", "Aminoglycoside", "Clindamycin Erythromycin Streptogramin B"))
  ) |>
  summarise(
    med = median(pr),
    lci = quantile(pr, 0.025),
    uci = quantile(pr, 0.975)
  ) |>
  ggplot(aes(t, med,
    ymin = lci, ymax = uci,
    color = exposure, fill = exposure
  )) +
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.3) +
  xlim(c(-4, 30)) +
  facet_grid(~outcome, labeller = labeller(outcome = label_wrap_gen(25))) +
  theme_bw() +
  labs(y = "Prevalence", linetype = "Exposure") +
  scale_color_manual(values = viridis_pal(option = "A")(6)[c(2, 4)]) +
  scale_fill_manual(values = viridis_pal(option = "A")(6)[c(2, 4)]) +
  theme(legend.position = "bottom") +
  labs(color = "Exposure", fill = "Exposure",
    title = "AMR gene prevalence vs time", subtitle = "Ceftrixaone vs Co-trimoxazole")


p_3b <-
  df_ressim |>
  filter(
    exposure %in% c("ciprofloxacin", "amoxicillin"),
    outcome %in% c("cephalosporin", "aminoglycoside", "clindamycin_erythromycin_streptogramin b")
  ) |>
  group_by(t, exposure, outcome) |>
  mutate(
    outcome = str_to_title(gsub("_", " ", outcome)),
    exposure = str_to_title(exposure),
    exposure = factor(exposure, levels = c("Ciprofloxacin", "Amoxicillin")),
    outcome = factor(outcome, c("Cephalosporin", "Aminoglycoside", "Clindamycin Erythromycin Streptogramin B"))
  ) |>
  summarise(
    med = median(pr),
    lci = quantile(pr, 0.025),
    uci = quantile(pr, 0.975)
  ) |>
  ggplot(aes(t, med,
    ymin = lci, ymax = uci,
    color = exposure, fill = exposure
  )) +
  geom_line() +
  geom_ribbon(color = NA, alpha = 0.3) +
  xlim(c(-4, 30)) +
  facet_grid(~outcome, labeller = labeller(outcome = label_wrap_gen(25))) +
  theme_bw() +
  labs(y = "Prevalence", linetype = "Exposure") +
  scale_color_manual(values = viridis_pal(option = "B")(6)[c(3, 5)]) +
  scale_fill_manual(values = viridis_pal(option = "B")(6)[c(3, 5)]) +
  theme(legend.position = "bottom") +
  labs(color = "Exposure", fill = "Exposure",
    title = "AMR gene prevalence vs time", subtitle = "Ciprofloxacin vs Amoxicillin")



# -----------
p_1 / (p_2a + p_2b) + (p_3a + p_3b) + plot_annotation(tag_levels = "A")

ggsave(here("eccmid_plot.pdf"), width = 11, height =10)

ggsave(here("eccmid_plot.svg"), width = 11, height =10)
