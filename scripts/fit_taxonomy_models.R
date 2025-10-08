library(tidyverse)
library(phyloseq)
library(here)
library(blantyreESBL)
library(microViz)
library(janitor)
library(cmdstanr)
library(bayesplot)

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

biom_bact <- subset_taxa(biom, Kingdom == "Bacteria")

otu_abs <-
  biom_bact@tax_table@.Data |>
  as_tibble(rownames = "tax_id") |>
  left_join(
    biom_bact@otu_table@.Data |>
      as_tibble(rownames = "tax_id"),
    by = join_by(tax_id)
  ) |>
  pivot_longer(-c(tax_id, Kingdom, Phylum, Class, Order, Family, Genus, Species)) |>
  left_join(
    btESBL_stoolESBL |>
      select(lab_id, pid, arm, visit, t),
    by = join_by(name == lab_id)
  ) |>
  filter(!is.na(visit)) |>
  filter(!name %in% toupper(samples_to_drop))


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


rank_bacterial_phyla <-
  otu_abs |>
  group_by(Phylum, name) |>
  summarise(value = sum(value)) |>
  group_by(Phylum) |>
  summarise(median_abundance = median(value)) |>
  arrange(-median_abundance)


rank_proteobacteria_orders <-
  otu_abs |>
  filter(Phylum == "Proteobacteria") |>
  group_by(Order, name) |>
  summarise(value = sum(value)) |>
  group_by(Order) |>
  summarise(median_abundance = median(value)) |>
  arrange(-median_abundance)

rank_enterobacterales_genus <-
  otu_abs |>
  filter(Order == "Enterobacterales") |>
  group_by(Genus, name) |>
  summarise(value = sum(value)) |>
  group_by(Genus) |>
  summarise(median_relative_abundance = median(value)) |>
  arrange(-median_relative_abundance)

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
  left_join(
    otu_abs |>
      group_by(name, pid, arm, visit, t) |>
      summarise(depth = sum(value)),
    otu_abs |>
      group_by(name, Phylum) |>
      summarise(reads = sum(value)) |>
      filter(Phylum != "", Phylum %in% rank_bacterial_phyla$Phylum[1:3]) |>
      mutate(Phylum = paste0("phylum_", Phylum)) |>
      pivot_wider(id_cols = name, names_from = Phylum, values_from = reads)
  ) |>
  left_join(
    otu_abs |>
      filter(Phylum == "Proteobacteria") |>
      group_by(name, Order) |>
      summarise(reads = sum(value)) |>
      filter(Order != "", Order %in% rank_proteobacteria_orders$Order[1:10]) |>
      mutate(Order = paste0("order_", Order)) |>
      pivot_wider(id_cols = name, names_from = Order, values_from = reads)
  ) |>
  left_join(
    otu_abs |>
      filter(Order == "Enterobacterales") |>
      group_by(name, Genus) |>
      summarise(reads = sum(value)) |>
      filter(Genus != "", Genus %in% rank_enterobacterales_genus$Genus[1:11]) |>
      mutate(Genus = paste0("genus_", Genus)) |>
      pivot_wider(id_cols = name, names_from = Genus, values_from = reads)
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
  arrange(pid, t) |>
  left_join(
    btESBL_stoolESBL |>
      transmute(
        name = lab_id,
        sample_type = case_when(
          is.na(sample_type) ~ 1,
          !is.na(sample_type) & sample_type == "stool" ~ 1,
          TRUE ~ 0
        )
      ),
    by = join_by(name)
  ) |>
  relocate(c(pid, sample_type))

outcome_vars <- names(df_mod)[grepl("phylum|order|genus", names(df_mod))]

fitlistout <- list()
drawslistout <- list()
summarylistout <- list()


mod <- cmdstan_model(here("negbin_gp_nopreds_stoolvswab.stan"))

problem_models <- outcome_vars[grepl("Firm|Hypho|Desulf|Burkhol|Serrat", outcome_vars)]

for (i in seq_len(length(outcome_vars))) {
  print(paste0("Generating model data for outcome var ", outcome_vars[i]))
  print(paste0("Variable ", i, " of ", length(outcome_vars)))

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
      y = df_mod[[outcome_vars[i]]],
      d = log(df_mod$depth),
      stool = df_mod$sample_type
    )

  print("Extracting draws and saving to list")


  fit <- mod$sample(
    data = stan_data_list,
    chains = 4,
    parallel_chains = 4,
    # iter_warmup = 1000,
    # iter_sampling = 1000,
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

write_rds(mod_sum_df, here("data_processed/gp_taxonomy_mod_sum_df.rda"), compress = "gz")

write_rds(fitlistout, here("data_processed/gp_taxonomy_fitlistout.rda"), compress = "gz")

write_rds(drawslistout, here("data_processed/gp_taxonomy_drawslistout.rda"), compress = "gz")

write_rds(outcome_vars, here("data_processed/gp_taxonomy_outcome_vars.rda"))

diagnosticsumlistout <- list()

for (i in seq_len(length(outcome_vars))) {
  dfout <- bind_rows(fitlistout[[i]]$diagnostic_summary())
  dfout$outcome_var <- outcome_vars[i]

  diagnosticsumlistout[[i]] <- dfout
}

mod_diagnostics_df <-
  bind_rows(diagnosticsumlistout)


write_rds(mod_diagnostics_df,
  here("data_processed/gp_taxonomy_mod_diagnostics_df.rda"),
  compress = "gz"
)

write_rds(outcome_vars, here("data_processed/gp_taxonomy_outcome_vars.rda"))
