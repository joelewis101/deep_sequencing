library(tidyverse)
library(phyloseq)
library(here)
library(blantyreESBL)
library(microViz)
library(janitor)
library(cmdstanr)
library(bayesplot)

library(readxl)

biom <- import_biom(here("data_raw/vivien_taxonomy_reports/table.biom"))

samplist <-
  read_xls(here("data_raw/sample_management/submission/2021-10-15_CGR-LIMS-submission.xls")) |>
  janitor::clean_names() |>
  mutate(
    name =
      if_else(grepl("NC", sample_name),
        sample_name,
        gsub("_.+", "", sample_name)
      ),
    row = gsub("[0-9]", "", plate_position_well_a1_b1_c1_etc_or_a1_a2_a3_etc),
    col = gsub("[A-Z]", "", plate_position_well_a1_b1_c1_etc_or_a1_a2_a3_etc)
  )

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


reads <-
  biom@tax_table@.Data |>
  as_tibble(rownames = "tax_id") |>
  left_join(
    biom@otu_table@.Data |>
      as_tibble(rownames = "tax_id"),
    by = join_by(tax_id)
  ) |>
  pivot_longer(-c(tax_id, Kingdom, Phylum, Class, Order, Family, Genus, Species)) |>
  group_by(name) |>
  summarise(reads = sum(value)) |>
  mutate(type = if_else(grepl("NC", name), "NC", "patient sample"))

reads <-
  reads |>
  left_join(
    select(
      samplist,
      c(name, plate_name_number, row, col, concentration_ng_ul)
    ),
    by = join_by(name == name)
  ) |>
  mutate(conc = as.numeric(concentration_ng_ul))

reads |>
  ggplot(aes(conc, reads)) +
  geom_point() +
  facet_wrap(~ type, scales = "free") +
  labs(y = "Nanodrop conc.")

reads |>
  ggplot(aes(type, reads)) +
  geom_jitter() 

reads |> filter(reads < 3e5, conc > 10) |> as.data.frame()

reads |> filter(grepl("NC", name)) |> as.data.frame()

reads |>
  filter(plate_name_number == 3) |>
  mutate(col = as.numeric(col)) |>
  ggplot(aes(col, row, fill = reads, label = name)) +
  geom_tile() +
  geom_text(size = 2) +
  scale_fill_viridis_c()
