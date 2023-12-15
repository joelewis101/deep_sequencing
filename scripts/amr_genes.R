library(tidyverse)
library(here)
library(janitor)
library(pheatmap)
library(ggplotify)
library(blantyreESBL)
library(FactoMineR)
library(factoextra)

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
  ungroup() |>
  select(-c(class, subclass, gene_symbol)) |>
  as.data.frame()

rownames(df2) <- make_clean_names(df$gene_symbol)

pl <-
  as.ggplot(pheatmap(df2))

ggsave("~/tmp/plot.pdf", pl, width = 40, height = 40)

# PCA on resistome, color by arm / abx exposure


df_pca <-
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

df_pca |>
  filter(is.na(pid))

# ok what we got

df_pca |>
  ggplot(aes(visit)) +
  geom_bar() +
  facet_wrap(~arm, ncol = 1)


p <-
  PCA(df_pca |>
    select(-c(
      "name",
      "pid",
      "arm",
      "visit",
      "ESBL",
      "enroll_date",
      "data_date"
    )))
#


armz <- df_pca$arm

p_ind <- bind_cols(
  p$ind$coord,
  arm = df_pca$arm,
  visit= df_pca$visit
)


p_ind |>
  ggplot(aes(Dim.1, Dim.2, color = as.character(arm))) +
  geom_point() +
stat_ellipse() +
  coord_cartesian(
    xlim = c(-5, 5),
    ylim = c(-5, 5)
  ) +
facet_grid(arm ~ visit)

# No real difference in resistome between arms or timepoints in this ax

# what about abx exposure - is probably not going to make *much* of a difference

# tomorrow - plot % of any beta lactamase gene
