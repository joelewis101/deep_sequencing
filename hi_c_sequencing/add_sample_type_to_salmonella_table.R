library(tidyverse)
library(blantyreESBL)

df <- read_tsv("Salmonella_plasmid_coverage_table.tab")


df |>
mutate(lab_id = gsub("^.*-","",Sample)) |>
mutate(lab_id = gsub("_.*$","", lab_id)) |>
left_join( select(btESBL_stoolESBL, lab_id, sample_type)) -> df

write_tsv(df, "Salmonella_plasmid_coverage_table_with_sample_type.tab")

