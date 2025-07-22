# prepare spreasheet for upload to stinking CGR LIMS

library(tidyverse)
library(readxl)

dat <-   
  read_csv(
            "data_raw/sample_management/samples_for_deep_seq.csv"
  )

template <- read_xls("data_raw/sample_management/GeneSifter_InputGrid.xls")

n_plates <- ceiling(nrow(dat)/96)

dat %>% 
  transmute(
     `Sample number` = 1:nrow(dat),                                         
     `Tube label` = NA_character_,
     `Sample name` = seqid,
     `Plate name/number` = rep(1:n_plates, each = 96)[1:nrow(dat)],
     `Plate position (well A1, B1, C1 etc or A1, A2, A3 etc.)` = 
       rep(
         paste0(rep(LETTERS[1:8], 12), rep(1:12, each = 8)),
         n_plates)[1:nrow(dat)],
     `Organism` = NA_character_,
     `Quantification method` = "Nanodrop",
     `Concentration (ng/ul)` = if_else(
       !is.na(dna_quantity),
       as.character(dna_quantity),
       "-"),
     `Purification method` = "QIAGEN spin column",
     `What type of buffer is the sample suspended in?` = "QIAGEN AE",
     `NanoDrop 260/280 ratio` = if_else(
       !is.na(a260_280),
       as.character(a260_280),
       "-"),
     `NanoDrop 260/230 ratio` = "-",
     `Sample volume (ul)` = 10,
     `Possible contaminating organism` = NA_character_
  ) %>% 
  write_csv("data_raw/sample_management/2021-10-15_CGR-LIMS-submission.csv")
