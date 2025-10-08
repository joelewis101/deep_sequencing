# This script will clean up the data and fit the AMR model and save the outputs 
# diagnostics are plotted in model-amr-diagnostics.qmd
# and resistome.qmd plots the overall analysis 
# for gaussian process model, model-ame-diagnostics-GP.qmd plots diagnostics
# and buildin-correlated-random-eff-mod.qmd plots results



library(tidyverse)
library(here)
library(janitor)
library(blantyreESBL)
library(patchwork)
library(kableExtra)
library(cmdstanr)
library(bayesplot)


df <-
  read_tsv(
    here(
      "data_raw/ECO_AMRfinder_all_metagenomes_default_17_11_23/AMR_binary_matrix.tsv"
    )
  ) |>
  janitor::clean_names()

samples_to_drop <- read_lines(here("data_processed/samples_to_drop.txt"))

df <-
  df |>
  pivot_longer(-c(class, subclass, gene_symbol)) |>
  mutate(
    name = str_extract(
      name,
      "sample_[0-9]{1,3}_(.*)_(:?uncorrected_amr|amr)",
      group = 1
    )
  ) |>
  mutate(
    subclass = if_else(
      grepl("cfi|cbl|crx|cep|cfx", gene_symbol),
      paste0(subclass, "_b"),
      subclass
    )
  ) |>
  mutate(name = gsub("_uncorrected", "", name)) |>
  pivot_wider(
    id_cols = c(class, subclass, gene_symbol),
  ) |>
  filter(class == "BETA-LACTAM")

df2 <-
  df |>
  mutate(
    subclass = tolower(subclass),
    subclass = gsub("-|/", "_", subclass)
  ) |>
  select(-c(class, gene_symbol)) |>
  pivot_longer(-subclass) |>
  filter(!grepl("^nc", name)) |>
  mutate(name = gsub("_.*$", "", name)) |>
  filter(value == 1) |>
  select(-value) |>
  unique() |>
  pivot_wider(
    names_from = subclass,
    values_from = subclass,
    values_fn = length,
    values_fill = 0
  ) |>
  left_join(
    btESBL_stoolESBL |>
      mutate(lab_id = tolower(lab_id),
      assess_type = as.numeric(difftime(data_date, enroll_date, unit = "days"))) |>
      select(pid, lab_id, assess_type),
    by = c("name" = "lab_id")
  ) 

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
  # add in abs to the df2 df in prep for join
  df2 |>
  cross_join(
    tibble(exps = unique(exps$name))
  ) |>
  # then rolling join by pid, ab
  left_join(exps,
    by = join_by(pid, exps == name, closest(assess_type > start))
  ) |>
  mutate(exps = paste0("exp_", exps)) |>
  rowwise() |>
  mutate(
    time_since_exp =
      case_when(
        is.na(end) ~ NA,
        assess_type > start & assess_type <= end ~ 0,
        assess_type > end ~ assess_type - end,
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
  arrange(pid, assess_type) |>
  filter(!name %in% tolower(samples_to_drop)) |>
  left_join(
    btESBL_stoolESBL |>
      transmute(
        name = tolower(lab_id),
        sample_type = case_when(
          is.na(sample_type) ~ 1,
          !is.na(sample_type) & sample_type == "stool" ~ 1,
          TRUE ~ 0
        )
      ),
    by = join_by(name)
  ) |>
  relocate(c(pid, assess_type, sample_type))

df_mod |>
  select(pid, matches(("exp_"))) |>
  pivot_longer(-pid) |>
  filter(value != -1) |>
  select(pid, name) |>
  unique() |>
  count(name) |>
  arrange(desc(n)) |>
  as.data.frame()


# fit the multilevel Gaussian process model and save

# fit for all


outcome_vars <-
  names(df_mod)[!grepl("exp|name|pid|assess_type|sample_type", names(df_mod))]

fitlistout <- list()
drawslistout <- list()
summarylistout <- list()

for (i in seq_len(length(outcome_vars))) {
  print(paste0("Generating model data for outcome var ", outcome_vars[i]))
  print(paste0("Variable ", i, " of ", length(outcome_vars)))


  stan_data_list <-
    list(
      n = nrow(df_mod),
      n_participants = length(unique(df_mod$pid)),
      n_covariates = 5,
      t = scale(df_mod$assess_type)[, 1],
      N = df_mod |>
        group_by(pid) |>
        summarise(n = n()) |>
        pull(n),
      t_e = matrix(
        c(
          df_mod$exp_cefo / sd(df_mod$assess_type),
          df_mod$exp_hosp / sd(df_mod$assess_type),
          df_mod$exp_cotri / sd(df_mod$assess_type),
          df_mod$exp_cipro / sd(df_mod$assess_type),
          # df_mod$exp_tb,
          df_mod$exp_amoxy / sd(df_mod$assess_type)
        ),
        ncol = 5
      ),
      stool = df_mod$sample_type,
      y = df_mod[[outcome_vars[i]]]
    )

  print("Extracting draws and saving to list")

  fit <- mod2$sample(
    data = stan_data_list,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 100,
    adapt_delta = 0.99
  )

  fitlistout[[i]] <- fit

  d <- fit$draws()

  drawslistout[[i]] <- d

  df_sum_out <- mcmc_intervals_data(d, prob_outer = 0.95)
  df_sum_out$outcome_var <- outcome_vars[i]
  summarylistout[[i]] <- df_sum_out
}

mod_sum_df <- bind_rows(summarylistout)

write_rds(mod_sum_df, 
  here("data_processed/resfindermodel_summary-df_bacteroidesbls-stratifed.rda"), compress = "gz")

write_rds(fitlistout,
  here("data_processed/resfindermodel_fitted-model-list_bacteroidesbls-stratified.rda"), compress = "gz")

write_rds(drawslistout,
  here("data_processed/resfindermodel_posterior-draws-list_bacteroidesbls-stratified.rda"), compress = "gz")


diagnosticsumlistout <- list()

for (i in seq_len(length(outcome_vars))) {

  dfout <- bind_rows(fitlistout[[i]]$diagnostic_summary())
  dfout$outcome_var <- outcome_vars[i]

  diagnosticsumlistout[[i]] <- dfout

}

mod_diagnostics_df <-
  bind_rows(diagnosticsumlistout)


write_rds(mod_diagnostics_df, 
  here("data_processed/resfindermodel_diagnostics_bacteroidesbls-stratified.rda"), compress = "gz")

write_rds(outcome_vars, here("data_processed/resfindermodel_outcomes-var_bacteroidesbls-stratified.rda"))


