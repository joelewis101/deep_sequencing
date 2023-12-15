# plot shipped boxes to figure out NCs and give them unique identfifiers

library(tidyverse)
library(here)
library(janitor)
library(lubridate)


df_manifest <-
  read_csv(here("data_raw/sample_management/sample_manifests.csv")) %>%
  mutate(sample_id =
      if_else(grepl("NC", sample_id),
                 paste0(sample_id, "_", row,  column),
         sample_id)) %>%
  # correct clear manifest errors
  # if 1) sample id does not exist
  # 2) it can be clearly inferred from extraction list then
  # change
  mutate(sample_id =
         case_when(sample_id == "CAH12H" ~ "CAB12H",
                   TRUE ~ sample_id)
         )

df_extraction <-
  read_csv(here("data_raw/sample_management/dna_extraction_records.csv")) %>%
  # omit duplicated records
  filter(!extraction_id %in% 317:340) %>%
  mutate(extraction_date = dmy(extraction_date)) %>%
  # correct manually recognised errors
  # These are identfied by samples in the manifest that have
  # no corresponding sample in the extraction list then
  # manually checking them - if the patient id does not exist in
  # the sample database then change to the sample from the manifest
  # if that is easily identifiable
  mutate(sample_id =
         case_when(extraction_id == 830 ~ "CAI101",
                   extraction_id == 646 ~ "CAD114",
                   TRUE ~ sample_id)
         ) %>%
  # link NCs to manifest
  mutate(sample_id =
         case_when(extraction_id == 836 ~ "NC_A7",
                   extraction_id == 879 ~ "NC_E4",
                   extraction_id == 644 ~ "NC_I9",
                   extraction_id == 751 ~ "NC_H5",
                   extraction_id == 779 ~ "NC_D1",
                   extraction_id == 929 ~ "NC_B6",
                   extraction_id == 946 ~ "NC_B6",
                   extraction_id == 247 ~ "NC_C6",
                   TRUE ~ sample_id)
         )

# figure out NCs

df_manifest <-
  read_csv(here("data_raw/sample_management/sample_manifests.csv")) %>%
  mutate(sample_id =
      if_else(grepl("NC", sample_id),
                 paste0(sample_id, "_", row,  column),
         sample_id)) %>%
  # correct clear manifest errors
  # if 1) sample id does not exist
  # 2) it can be clearly inferred from extraction list then
  # change
  mutate(sample_id =
         case_when(sample_id == "CAH12H" ~ "CAB12H",
                   TRUE ~ sample_id)
         )

df_extraction <-
  read_csv(here("data_raw/sample_management/dna_extraction_records.csv")) %>%
  # omit duplicated records
  filter(!extraction_id %in% 317:340) %>%
  mutate(extraction_date = dmy(extraction_date)) %>%
  # correct manually recognised errors
  # These are identfied by samples in the manifest that have
  # no corresponding sample in the extraction list then
  # manually checking them - if the patient id does not exist in
  # the sample database then change to the sample from the manifest
  # if that is easily identifiable
  mutate(sample_id =
         case_when(extraction_id == 830 ~ "CAI101",
                   extraction_id == 646 ~ "CAD114",
                   TRUE ~ sample_id)
         ) %>%
  # link NCs to manifest
  mutate(sample_id =
         case_when(extraction_id == 836 ~ "NC_A7",
                   extraction_id == 879 ~ "NC_E4",
                   extraction_id == 644 ~ "NC_I9",
                   extraction_id == 751 ~ "NC_H5",
                   extraction_id == 779 ~ "NC_D1",
                   extraction_id == 929 ~ "NC_B6",
                   extraction_id == 946 ~ "NC_B6",
                   extraction_id == 247 ~ "NC_C6",
                   TRUE ~ sample_id)
         )

# figure out NCs

box_n <- 5
df_manifest %>%
  left_join(filter(df_extraction, !grepl("NC$", sample_id))) %>%
  group_by(sample_id) %>%
  mutate(labz =
         case_when(grepl("NC", sample_id) ~ sample_id,
                   TRUE ~ paste(unique(extraction_id), collapse = ","))
         ) %>%
  filter(box == box_n) %>%
  ggplot(aes(column,
             fct_rev(row),
             fill = as.character(extraction_date),
             shape = grepl("NC", sample_id))) +
  geom_tile() +
  geom_point() +
  scale_x_continuous(breaks = 1:9) +
  labs(title = paste0("Box ", box_n)) +
  geom_text(aes(label = labz), nudge_y = 0.2, size = 3)

library(blantyreESBL)
