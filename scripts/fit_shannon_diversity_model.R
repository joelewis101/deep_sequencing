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

biom <- import_biom(here("data_raw/vivien_taxonomy_reports/table.biom"))

samples_to_drop <- read_lines(here("data_processed/samples_to_drop.txt"))

# add metadata
biom <-
  ps_join(biom,
    blantyreESBL::btESBL_stoolESBL |>
      select(lab_id, pid, arm, visit, ESBL) |>
      mutate(arm = as.character(arm)) |>
      left_join(
        btESBL_participants |>
          transmute(
            pid = pid,
            abx = recieved_prehosp_ab,
            visit = 0
          ),
        by = join_by(pid, visit)
      ),
    match_sample_names = "lab_id"
  )

# mung headings
biom@tax_table@.Data <-
  substring(biom@tax_table@.Data, 4)

colnames(biom@tax_table@.Data) <-
  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# filter homo sapients / streptophyta


taxa_to_filter <-
  c(
    biom@tax_table@.Data |>
      as.data.frame() |>
      filter(grepl("sapiens", Species)) |>
      rownames(),
    biom@tax_table@.Data |>
      as.data.frame() |>
      filter(if_any(everything(), ~ grepl("Streptophyta", .x))) |>
      rownames()
  )

biom <- prune_taxa(!taxa_names(biom) %in% taxa_to_filter, biom)

# remove those with low reads
biom <-
  prune_samples(
    !sample_names(biom) %in% samples_to_drop,
    biom
  )


# alpha diversity


exps <-
  btESBL_exposures |>
  select(-died) |>
  pivot_longer(-c(pid, assess_type)) |>
  filter(value == 1) |>
  group_by(pid, name) |>
  summarise(
    start = min(assess_type),
    end = max(assess_type)
  )

df_mod <-
  estimate_richness(biom, measures = "Shannon") |>
  as_tibble(rownames = "lab_id") |>
  left_join(
    btESBL_stoolESBL |>
      select(lab_id, pid, arm, visit, t),
    by = join_by(lab_id)
  ) |>
  cross_join(
    tibble(exps = unique(exps$name))
  ) |>
  # then rolling join by pid, ab
  left_join(exps,
    by = join_by(pid, exps == name, closest(t > start))
  ) |>
  mutate(exps = paste0("exp_", exps)) |>
  rowwise() |>
  mutate(
    time_since_exp =
      case_when(
        is.na(end) ~ NA,
        t > start & t <= end ~ 0,
        t > end ~ t - end,
        TRUE ~ -999
      )
  ) |>
  select(-c(start, end)) |>
  pivot_wider(
    id_cols = !matches("exps"),
    names_from = exps,
    values_from = time_since_exp
  ) |>
  mutate(across(matches("exp"), ~ if_else(is.na(.x), -1, .x))) |>
  arrange(pid, t)


mod <- cmdstan_model(here("lm.stan"))


stan_data_list <-
  list(
    n = nrow(df_mod),
    n_participants = length(unique(df_mod$pid)),
    n_covariates = 5,
    t = scale(df_mod$t)[, 1],
    N = df_mod |>
      group_by(pid) |>
      summarise(n = n()) |>
      pull(n),
    t_e = matrix(
      c(
        df_mod$exp_cefo / sd(df_mod$t),
        df_mod$exp_hosp / sd(df_mod$t),
        df_mod$exp_cotri / sd(df_mod$t),
        df_mod$exp_cipro / sd(df_mod$t),
        # df_mod$exp_tb,
        df_mod$exp_amoxy / sd(df_mod$t)
      ),
      ncol = 5
    ),
    y = df_mod$Shannon
  )


fit <- mod$sample(
  data = stan_data_list,
  chains = 4,
  parallel_chains = 4,
  # iter_warmup = 1000,
  # iter_sampling = 1000,
  refresh = 100,
  # adapt_delta = 0.99
)


d <- fit$draws()


df_sum_out <- mcmc_intervals_data(d, prob_outer = 0.95)


write_rds(df_sum_out, here("data_processed/diversity_mod_sum_df.rda"), compress = "gz")

write_rds(fit, here("data_processed/diversity_mod_fit.rda"), compress = "gz")

write_rds(d, here("data_processed/diversity_mod_draws.rda"), compress = "gz")

mod_diagnostics_df <-
fit$diagnostic_summary() |>
  as.data.frame()


write_rds(mod_diagnostics_df,
  here("data_processed/diversity_mod_diagnostics_df.rda"),
  compress = "gz"
)
