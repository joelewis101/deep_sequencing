library(blantyreESBL)
library(tidyverse)
library(here)

df <- read_tsv(here("data_raw/ECO_all_metagenomic_read_stats_for_JL_11_02_25.tab"))

df <-
  df |>
  mutate(sample_id = str_extract(file, "(?<=^Sample_[0-9]{1,3}-).{6}"))
  
    semi_join(
    btESBL_stoolESBL,
    by = join_by(sample_id == lab_id)
  )


# days of hospitalisation
btESBL_exposures |>
  select(assess_type, hosp, pid) |>
  filter(hosp == 1) |>
  group_by(pid) |>
  filter(assess_type == max(assess_type)) |>
  transmute(
    pid = pid,
    days_of_hospitalisation = assess_type
  )


df_out <- 
  df |>
  select(file, sample_id) |>
      left_join( 
  btESBL_stoolESBL |>
      left_join(
        btESBL_exposures |>
          select(assess_type, hosp, pid) |>
          filter(hosp == 1) |>
          group_by(pid) |>
          filter(assess_type == max(assess_type)) |>
          transmute(
            pid = pid,
            days_of_hospitalisation = assess_type
          ),
          by = join_by(pid == pid)
          ) |>
        mutate(location_of_collection = 
          case_when(
            is.na(days_of_hospitalisation) ~ "community",
            !is.na(days_of_hospitalisation) & (t <= days_of_hospitalisation) ~ "hospital",
            !is.na(days_of_hospitalisation) & (t > days_of_hospitalisation) ~ "community")
            ) |>
      transmute(
        participant_id = pid,
        sample_id = lab_id,
        collection_date = data_date,
        sample_type = sample_type,
        location_of_collection = location_of_collection
      ) ,
    by = join_by(sample_id == sample_id)
  ) |>
  mutate(sample_type = if_else(is.na(sample_type), "unknown", sample_type))

write_csv(df_out, here("data_processed/metadata_for_ena_submission.csv"))
