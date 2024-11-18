
library(tidyverse)
library(phyloseq)
library(here)
library(blantyreESBL)
library(microViz)
library(janitor)
library(cmdstanr)
library(bayesplot)

biom <- import_biom(here("data_raw/vivien_taxonomy_reports/table.biom"))


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

distbiom <- distance(otu_table(biom), method = "bray")

pcoa_biom <-
ordinate( otu_table(biom), "PCoA", "bray")

plot_ordination(biom, pcoa_biom, color = "ESBL") +
  facet_grid(arm ~ visit)

# alpha diversity
estimate_richness(biom, measures = "Shannon") |>
  as_tibble(rownames = "lab_id") |>
  left_join(
    btESBL_stoolESBL |>
      select(lab_id, pid, arm, visit) |>
      # mutate(arm = as.character(arm)) |>
      left_join(
        btESBL_participants |>
          transmute(
            pid = pid,
            abx = recieved_prehosp_ab,
            visit = 0
          ),
        by = join_by(pid, visit)
      ),
    by = join_by("lab_id")
  ) |>
  ggplot(aes(visit, Shannon, group = interaction(visit, arm), color = as.factor(arm))) +
  geom_boxplot()

# extract taxa, get to relative abundance and plot

biom_bact <- subset_taxa(biom, Kingdom == "Bacteria")

biom_bact_r <- transform_sample_counts(biom_bact, function(x) x / sum(x))




# biom_bact_r <- tax_glom(biom_bact_r, taxrank = "Phylum")

otu_r_df <-
  biom_bact_r@tax_table@.Data |>
  as_tibble(rownames = "tax_id") |>
  left_join(
    biom_bact_r@otu_table@.Data |>
      as_tibble(rownames = "tax_id"),
    by = join_by(tax_id)
  ) |>
  pivot_longer(-c(tax_id, Kingdom, Phylum, Class, Order, Family, Genus, Species)) |>
  left_join(
    btESBL_stoolESBL |>
      select(lab_id, pid, arm, visit) |>
      # mutate(arm = as.character(arm)) |>
      left_join(
        btESBL_participants |>
          transmute(
            pid = pid,
            admission_abx = recieved_prehosp_ab,
            visit = 0
          ),
        by = join_by(pid, visit)
      ),
    by = join_by(name == lab_id)
  ) |>
  filter(!is.na(visit))

rank_bacterial_phyla <-
  otu_r_df |>
  group_by(Phylum, name) |>
  summarise(value = sum(value)) |>
  group_by(Phylum) |>
  summarise(median_relative_abundance = median(value)) |>
  arrange(-median_relative_abundance)

otu_r_df |>
  group_by(name, pid, arm, visit, Phylum) |>
  summarise(value = sum(value)) |>
  filter(Phylum %in% rank_bacterial_phyla$Phylum[1:10]) |>
  mutate(Phylum = factor(Phylum, levels = rank_bacterial_phyla$Phylum)) |>
  ggplot(aes(Phylum, value, fill = Phylum)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(alpha = 0.6, shape = 3) +
  facet_grid(arm ~ visit) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

rank_proteobacteria_orders <-
  otu_r_df |>
  filter(Phylum == "Proteobacteria") |>
  group_by(Order, name) |>
  summarise(value = sum(value)) |>
  group_by(Order) |>
  summarise(median_relative_abundance = median(value)) |>
  arrange(-median_relative_abundance)

otu_r_df |>
  filter(Phylum == "Proteobacteria") |>
  group_by(name, pid, arm, visit, Order) |>
  summarise(value = sum(value)) |>
  filter(Order %in% rank_proteobacteria_orders$Order[1:10]) |>
  mutate(Order = factor(Order, levels = rank_proteobacteria_orders$Order)) |>
  ggplot(aes(Order, value)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(alpha = 0.6, shape = 3) +
  facet_grid(arm ~ visit) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

rank_enterobacterales_genus <-
  otu_r_df |>
  filter(Order == "Enterobacterales") |>
  group_by(Genus, name) |>
  summarise(value = sum(value)) |>
  group_by(Genus) |>
  summarise(median_relative_abundance = median(value)) |>
  arrange(-median_relative_abundance)

otu_r_df |>
  filter(Order == "Enterobacterales") |>
  group_by(name, pid, arm, visit, Genus) |>
  summarise(value = sum(value)) |>
  filter(Genus %in% rank_enterobacterales_genus$Genus[1:11]) |>
  filter(Genus != "") |>
  mutate(Genus = factor(Genus, levels = rank_enterobacterales_genus$Genus)) |>
  ggplot(aes(Genus, value)) +
  geom_boxplot(outliers = FALSE) +
  geom_point(alpha = 0.6, shape = 3) +
  facet_grid(arm ~ visit) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 0.3))

# correlate phyla with AMR genes

# load resfinder stuff


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


df |>
  mutate(gene_symbol = make_clean_names(gene_symbol)) |>
  pivot_longer(-c(gene_symbol, class, subclass)) |>
  filter(!grepl("^nc", name)) |>
  mutate(
    name = toupper(gsub("_.*$", "", name)),
    subclass = tolower(gsub("/", "-", subclass))
  ) |>
  group_by(name, subclass) |>
  summarise(amr_present = if_else(any(value == 1), "Present", "Absent")) |>
  pivot_wider(id_cols = name, names_from = subclass, values_from = amr_present) |>
  left_join(
    otu_r_df |>
      group_by(name, pid, arm, visit, Phylum) |>
      summarise(value = sum(value)) |>
      filter(Phylum %in% rank_bacterial_phyla$Phylum[1:3]) |>
      mutate(Phylum = factor(Phylum, levels = rank_bacterial_phyla$Phylum)) |>
      ungroup() |>
      select(name, Phylum, value),
    by = join_by(name)
  ) |>
  pivot_longer(-c(name, Phylum, value), names_to = "amr_gene", values_to = "present") |>
  ggplot(aes(amr_gene, value, fill = present)) + 
  geom_boxplot(outliers = FALSE) +
  # geom_point(alpha = 0.6, shape = 3) +
  facet_wrap(~ Phylum, nrow = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
coord_flip()


df |>
  mutate(gene_symbol = make_clean_names(gene_symbol)) |>
  pivot_longer(-c(gene_symbol, class, subclass)) |>
  filter(!grepl("^nc", name)) |>
  mutate(
    name = toupper(gsub("_.*$", "", name)),
    subclass = tolower(gsub("/", "-", subclass))
  ) |>
  group_by(name, class) |>
  summarise(amr_present = if_else(any(value == 1), "Present", "Absent")) |>
  pivot_wider(id_cols = name, names_from = class, values_from = amr_present) |>
  left_join(
    otu_r_df |>
      group_by(name, pid, arm, visit, Phylum) |>
      summarise(value = sum(value)) |>
      filter(Phylum %in% rank_bacterial_phyla$Phylum[1:3]) |>
      mutate(Phylum = factor(Phylum, levels = rank_bacterial_phyla$Phylum)) |>
      ungroup() |>
      select(name, Phylum, value),
    by = join_by(name)
  ) |>
  pivot_longer(-c(name, Phylum, value), names_to = "amr_gene", values_to = "present") |>
  ggplot(aes(amr_gene, value, fill = present)) +
  geom_boxplot(outliers = FALSE, notch = FALSE) +
  # geom_point(alpha = 0.6, shape = 3) +
  facet_wrap(~Phylum, nrow = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

df |>
  mutate(gene_symbol = make_clean_names(gene_symbol)) |>
  pivot_longer(-c(gene_symbol, class, subclass)) |>
  filter(!grepl("^nc", name)) |>
  mutate(
    name = toupper(gsub("_.*$", "", name)),
    subclass = tolower(gsub("/", "-", subclass))
  ) |>
  group_by(name, class) |>
  summarise(amr_present = if_else(any(value == 1), "Present", "Absent")) |>
  pivot_wider(id_cols = name, names_from = class, values_from = amr_present) |>
  left_join(
    otu_r_df |>
      group_by(name, pid, arm, visit, Phylum) |>
      summarise(value = sum(value)) |>
      filter(Phylum %in% rank_bacterial_phyla$Phylum[1:3]) |>
      mutate(Phylum = factor(Phylum, levels = rank_bacterial_phyla$Phylum)) |>
      ungroup() |>
      select(name, Phylum, value),
    by = join_by(name)
  ) |>
  pivot_longer(-c(name, Phylum, value), names_to = "amr_gene", values_to = "present") |>
  group_by(amr_gene, Phylum, present) |>
  summarise(
    med = median(value),
    lci = quantile(value, 0.25),
    uci = quantile(value, 0.75)
  ) |>
  ggplot(aes(amr_gene, med, color = present, ymin = lci, ymax = uci)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(width = 0, position = position_dodge(width = 0.2)) +
  facet_wrap(~Phylum, nrow = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

# These all suck - heatmaps ?


df |>
  mutate(gene_symbol = make_clean_names(gene_symbol)) |>
  pivot_longer(-c(gene_symbol, class, subclass)) |>
  filter(!grepl("^nc", name)) |>
  mutate(
    name = toupper(gsub("_.*$", "", name)),
    subclass = tolower(gsub("/", "-", subclass))
  ) |>
  group_by(name, class) |>
  summarise(amr_present = if_else(any(value == 1), 1, 0)) |>
  pivot_wider(id_cols = name, names_from = class, values_from = amr_present) |>
  left_join(
    otu_r_df |>
      group_by(name, pid, arm, visit, Phylum) |>
      summarise(value = sum(value)) |>
      filter(Phylum %in% rank_bacterial_phyla$Phylum[1:3]) |>
      mutate(Phylum = factor(Phylum, levels = rank_bacterial_phyla$Phylum)) |>
      ungroup() |>
      pivot_wider(id_cols = name, names_from = Phylum, values_from = value),
    by = join_by(name)
  ) |>
  as.data.frame() -> cormat

rownames(cormat) <- cormat$name

cormat <- cor(select(cormat, -name), method = "spearman")

pheatmap(
  cormat[
    rownames(cormat) %in% c("Bacteroidetes", "Firmicutes", "Proteobacteria"),
    !dimnames(cormat)[[1]] %in% c("Bacteroidetes", "Firmicutes", "Proteobacteria")
  ]
)

df |>
  mutate(gene_symbol = make_clean_names(gene_symbol)) |>
  pivot_longer(-c(gene_symbol, class, subclass)) |>
  filter(!grepl("^nc", name)) |>
  mutate(
    name = toupper(gsub("_.*$", "", name)),
    subclass = tolower(gsub("/", "-", subclass))
  ) |>
  group_by(name, class) |>
  summarise(amr_present = if_else(any(value == 1), 1, 0)) |>
  pivot_wider(id_cols = name, names_from = class, values_from = amr_present) |>
  left_join(
otu_r_df |>
  filter(Phylum == "Proteobacteria") |>
  group_by(name, pid, arm, visit, Order) |>
  summarise(value = sum(value)) |>
  filter(Order %in% rank_proteobacteria_orders$Order[1:10]) |>
  mutate(Order = factor(Order, levels = rank_proteobacteria_orders$Order)) |>
      ungroup() |>
      pivot_wider(id_cols = name, names_from = Order, values_from = value),
    by = join_by(name)
  ) |>
  as.data.frame() -> cormat_order



rownames(cormat_order) <- cormat_order$name

cormat_order <- cor(select(cormat_order, -name), method = "spearman")

pheatmap(
  cormat_order[
    rownames(cormat_order) %in% rank_proteobacteria_orders$Order,
    !dimnames(cormat_order)[[1]] %in% rank_proteobacteria_orders$Order
  ]
)


df |>
  mutate(gene_symbol = make_clean_names(gene_symbol)) |>
  pivot_longer(-c(gene_symbol, class, subclass)) |>
  filter(!grepl("^nc", name)) |>
  mutate(
    name = toupper(gsub("_.*$", "", name)),
    subclass = tolower(gsub("/", "-", subclass))
  ) |>
  group_by(name, class) |>
  summarise(amr_present = if_else(any(value == 1), 1, 0)) |>
  pivot_wider(id_cols = name, names_from = class, values_from = amr_present) |>
  left_join(
otu_r_df |>
  filter(Order == "Enterobacterales") |>
  group_by(name, pid, arm, visit, Genus) |>
  summarise(value = sum(value)) |>
  filter(Genus %in% rank_enterobacterales_genus$Genus[1:11]) |>
  filter(Genus != "") |>
      ungroup() |>
      pivot_wider(id_cols = name, names_from = Genus, values_from = value),
    by = join_by(name)
  ) |>
  as.data.frame() -> cormat_genus



rownames(cormat_genus) <- cormat_genus$name

cormat_genus <- cor(select(cormat_genus, -name), method = "spearman")

pheatmap(
  cormat_genus[
    rownames(cormat_genus) %in% rank_enterobacterales_genus$Genus,
    !dimnames(cormat_genus)[[1]] %in% rank_enterobacterales_genus$Genus
  ]
)

# maybe leave this for assemblies ?
# next - model absolute abundance of reads vs time and antibiotic exposure

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
  filter(!is.na(visit))


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
left_join(
    otu_abs |>
      group_by(name, pid, arm, visit, t) |>
      summarise(depth = sum(value)),
    otu_abs |>
      filter(Phylum == "Proteobacteria") |>
      group_by(name, pid, arm, visit, t) |>
      summarise(reads = sum(value))
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

df_mod |>
  select(pid, matches(("exp_"))) |>
  pivot_longer(-pid) |>
  filter(value != -1) |>
  select(pid, name) |>
  unique() |>
  count(name) |>
  arrange(desc(n)) |>
  as.data.frame()

mod <- cmdstan_model(here("negbin-randintercept.stan"))

stan_data_list <-
  list(
    n = nrow(df_mod),
    n_participants = length(unique(df_mod$pid)),
    n_covariates = 5,
    pid = df_mod |>
      group_by(pid) |>
      mutate(int_id = cur_group_id()) |>
      pull(int_id),
    t_e = matrix(
      c(
        df_mod$exp_cefo,
        df_mod$exp_hosp,
        df_mod$exp_cotri,
        df_mod$exp_cipro,
        # df_mod$exp_tb,
        df_mod$exp_amoxy
      ),
      ncol = 5
    ),
    y = df_mod$reads,
    d = log(df_mod$depth)
  )

fit <- mod$sample(
  data = stan_data_list,
  chains = 4,
  parallel_chains = 4
  # iter_warmup = 1500,
  # iter_sampling = 3000
)



d <- fit$draws()

summary(d)

mcmc_trace(d,  regex_pars = "alpha|beta|sigma|tau")

mcmc_intervals(d, regex_pars = "alpha|beta|sigma", prob_outer = 0.95)

mcmc_intervals_data(d, regex_pars = "alpha|beta|sigma", prob_outer = 0.95)

mcmc_intervals(d, regex_pars = "gamma")

mcmc_intervals(d, regex_pars = "tau")


tibble(
y = stan_data_list$y,
summarise_draws(d[,,which(grepl("y_", dimnames(d)[[3]]))], joemed, joequant, joe50c)
) |> janitor::clean_names() |>
ggplot(aes(y, joemed, ymin = q5, ymax = q95)) +
geom_point() +
geom_ribbon(alpha = 0.2, color = NA) +
geom_ribbon(aes(ymin = x25_percent, ymax = x75_percent), alpha = 0.5, color = NA) + 
geom_abline(intercept = 0, slope = 1) +
coord_cartesian(ylim = c(0,1e+08))

# gausiian process


mean_assess_type <- mean(df_mod$assess_type)
sd_assess_type <- mean(df_mod$assess_type)


stan_data_list <-
  list(
    n = nrow(df_mod),
    n_participants = length(unique(df_mod$pid)),
    n_covariates = 5,
    N = df_mod |>
      group_by(pid) |>
      summarise(n = n()) |>
      pull(n),
    t = scale(df_mod$assess_type)[,1],
    t_e = matrix(
      c(
        df_mod$exp_cefo / sd(df_mod$assess_type),
        df_mod$exp_hosp / sd(df_mod$assess_type),
        df_mod$exp_cotri / sd(df_mod$assess_type),
        df_mod$exp_cipro /  sd(df_mod$assess_type),
        # df_mod$exp_tb,
        df_mod$exp_amoxy /   sd(df_mod$assess_type)
      ),
      ncol = 5
    ),
    y = df_mod$erythromycin_spiramycin_telithromycin
  )

stan_data_list <-
  list(
    n = nrow(df_mod),
    n_participants = length(unique(df_mod$pid)),
    n_covariates = 5,
    N = df_mod |>
      group_by(pid) |>
      summarise(n = n()) |>
      pull(n),
   t = scale(df_mod$t)[,1],
    t_e = matrix(
      c(
        df_mod$exp_cefo / sd(df_mod$t),
        df_mod$exp_hosp / sd(df_mod$t),
        df_mod$exp_cotri / sd(df_mod$t),
        df_mod$exp_cipro / sd(df_mod$t),
        # df_mod$exp_tb,
        df_mod$exp_amoxy /sd(df_mod$t)
      ),
      ncol = 5
    ),
    y = df_mod$reads,
    d = log(df_mod$depth)
  )
mod2 <- cmdstan_model(here("negbin_gp.stan"))


fit <- mod2$sample(
  data = stan_data_list,
  chains = 4,
  parallel_chains = 4
  # iter_warmup = 1500,
  # iter_sampling = 3000
)

d <- fit$draws()

d2 <- d[,,which(!grepl("y_|^p$", dimnames(d)[[3]]))]

summary(d2)



mcmc_intervals(d2, regex_pars = "alpha|beta|sigma|tau|length|phi", prob_outer = 0.95)

mcmc_trace(d2, regex_pars = "alpha|beta|sigma|tau|length|phi", prob_outer = 0.95)


mcmc_pairs(d, regex_pars = "alpha|beta|sigma|tau|length|phi", prob_outer = 0.95)

## h post pred check

joemed <- function(x) {
 median(x, na.rm = TRUE)
}

joequant <- function(x) {
 quantile2(x, na.rm = TRUE)
}


joe50c <- function(x) {
 quantile(x, c(0.25, 0.75), na.rm = TRUE)
}

tibble(
y = stan_data_list$y,
summarise_draws(d[,,which(grepl("y_", dimnames(d)[[3]]))], joemed, joequant, joe50c)
) |> janitor::clean_names() |>
ggplot(aes(y, joemed, ymin = q5, ymax = q95)) +
geom_point() +
geom_ribbon(alpha = 0.2, color = NA) +
geom_ribbon(aes(ymin = x25_percent, ymax = x75_percent), alpha = 0.5, color = NA) + 
geom_abline(intercept = 0, slope = 1) +
coord_cartesian(ylim = c(0,1e+08))


