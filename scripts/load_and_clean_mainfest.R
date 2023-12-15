# plot shpped boxes to figure out NCs and give them unique identfifiers

library(tidyverse)
library(here)
library(janitor)
library(lubridate)

cat("\nLoading data ...\n")

df_manifest <-
  read_csv(here("data_raw/sample_management/sample_manifests.csv")) %>%
  mutate(sample_id =
      if_else(grepl("NC", sample_id),
                 paste0(sample_id, "_", box, "_", row,  column),
         sample_id)) %>%
  # correct clear manifest errors
  # if 1) sample id does not exist
  # 2) it can be clearly inferred from extraction list then
  # change
  mutate(sample_id =
         case_when(sample_id == "CAH12H" ~ "CAB12H",
                   # based on extraction list
                   sample_id == "CAJ1-6" ~ "CAJ106",
                   sample_id == "CAH111" ~ "CAH11I",
                   TRUE ~ sample_id)
         )

df_extraction <-
  read_csv(here("data_raw/sample_management/dna_extraction_records.csv")) %>%
  select(!starts_with("X")) %>%
  # omit duplicated records
  filter(!extraction_id %in% c(317:340, 341:361, 269:292, 293:316)) %>%
  mutate(extraction_date = dmy(extraction_date)) %>%
  # correct manually recognised errors
  # These are identfied by samples in the manifest that have
  # no corresponding sample in the extraction list then
  # manually checking them - if the patient id does not exist in
  # the sample database then change to the sample from the manifest
  # if that is easily identifiable
  # I've been more lax about changing extraction ids -
  # the worst that can happen is the nandop / A280 ratio
  # will be wrong, right? right?
  mutate(sample_id =
         case_when(extraction_id == 830 ~ "CAI101",
                   extraction_id == 646 ~ "CAD114",
                   extraction_id == 613 ~ "CAC12D",
                   extraction_id == 626 ~ "CAD135",
                   extraction_id == 173 ~ "CAC13U",
                   extraction_id == 866 ~ "CAE10V",
                   extraction_id == 775 ~ "CAC11Y",
                   extraction_id == 738 ~ "CAL107",
                   extraction_id == 718 ~ "CAG11A",
                   extraction_id == 571 ~ "CAB15S",
                   extraction_id == 511 ~ "CAH11G",
                   TRUE ~ sample_id)
         ) %>%
  # link NCs to manifest
  mutate(sample_id =
         case_when(extraction_id == 836 ~ "NC_1_A7",
                   extraction_id == 879 ~ "NC_1_E4",
                   extraction_id == 644 ~ "NC_2_I9",
                   extraction_id == 751 ~ "NC_3_H5",
                   extraction_id == 779 ~ "NC_3_D1",
                   extraction_id == 929 ~ "NC_4_B6",
                   extraction_id == 946 ~ "NC_4_D5",
                   extraction_id == 247 ~ "NC_5_C6",
                   extraction_id == 105 ~ "NC_5_F6",
                   extraction_id == 81 ~ "NC_5_I3",
                   extraction_id == 635 ~ "NC_6_E3",
                   extraction_id == 611 ~ "NC_6_G9",
                   extraction_id == 587 ~ "NC_7_I6",
                   extraction_id == 151 ~ "NC_8_F3",
                   extraction_id == 268 ~ "NC_10_D4",
                   extraction_id == 515 ~ "NC_10_G2",
                   extraction_id == 491 ~ "NC_11_C5",
                   extraction_id == 64 ~ "NC_8_I2",
                   extraction_id == 22 ~ "NC_9_F1",
                   extraction_id == 223 ~ "NC_9_H3",
                   TRUE ~ sample_id))




df_out <-
  df_manifest %>%
  left_join(filter(df_extraction, !grepl("NC$", sample_id))) %>%
  # aim to manually match 1:1 maifest:extraction based on position on plates
  filter(
         !(box == 1 & row == "F" & column == 8 & extraction_id == 620),
         !(box == 1 & row == "H" & column == 3 &
           (extraction_id %in% c(217, 310))),
         !(box == 2 & row == "A" & column == 1 & extraction_id == 522),
         !(box == 2 & (row %in% c("B", "C", "D", "E")) &
            extraction_id > 709),
         !(box == 2 & row == "C" & column == 2 & extraction_id == 636),
         !(box == 2 & row == "H" & column == 1 &
           extraction_id %in% c(363, 678)),
         !(box == 2 & row == "H" & column == 8 &
           extraction_id %in% c(363, 671)),
         !(box == 2 & row == "I" & column == 1 & extraction_id %in% c(692)),
         !(box == 2 & row == "H" & column == 2 & extraction_id %in% c(362)),
         !(box == 3 & row == "C" & column == 6 & extraction_id %in% c(648)),
         !(!is.na(extraction_id) & box == 3 & (row %in% c("F", "G", "H")) &
            extraction_id < 710),
         !(box == 3 & row == "F" & column == 9 & extraction_id %in% c(838)),
         # these three are ambigous and nopt possible to resolve. I have picked
         # the ones with the least bad conc/A280 ratio
         !(box == 4 & row == "F" & column == 8 & extraction_id %in% c(881)),
         !(box == 4 & row == "F" & column == 9 & extraction_id %in% c(883)),
         !(box == 4 & row == "H" & column == 2 & extraction_id %in% c(897)),
         # back to business as usual with this lot(
         !(box == 5 & row == "F" & column == 1 & extraction_id %in% c(60)),
         !(box == 5 & row == "D" & column == 2 & extraction_id %in% c(59)),
         !(box == 6 & row == "A" & column == 2 & extraction_id %in% c(738)),
         !(box == 6 & row == "C" & column == 6 & extraction_id %in% c(866)),
         !(box == 6 & row == "D" & column == 3 & extraction_id %in% c(525)),
         !(box == 6 & row == "D" & column == 9 & extraction_id %in% c(634)),
         !(box == 6 & row == "E" & column == 2 & extraction_id %in% c(632)),
         !(box == 6 & row == "H" & column == 7 & extraction_id %in% c(723)),
         !(box == 6 & row == "I" & column == 1 & extraction_id %in% c(626)),
         # special case - this one is matched elsewhere more
         # reliably - restrict to one extraction id here then switch
         #that to NA below
         !(box == 7 & row == "C" & column == 7 & extraction_id %in% c(692)),
         # back to normal
         !(box == 7 & row == "E" & column == 5 & extraction_id %in% c(511)),
         !(box == 7 & row == "F" & extraction_id %in% c(565, 732)),
         !(box == 7 & row == "G" & column == 2 & extraction_id %in% c(557)),
         !(box == 8 & (row %in% c("B", "D", "E", "F")) &
         extraction_id < 100),
         !(box == 8 & (row %in% c("H", "I")) &
         !extraction_id %in% c(127, 47:64)),
         !(box == 10 & row == "E" & column == 3 & extraction_id %in% c(20)),
         !(box == 10 & (row %in% c("I")) &
         extraction_id < 460),
         !(box == 9 & row == "E" & column == 9 & extraction_id %in% c(157)),
         !(box == 9 & row == "I" & column == 5 & extraction_id %in% c(21)),
         ) %>%
  # remove the extraction metadata for those that are more reliably
  # matched elsewhere
  mutate(across(matches("extraction_id|dna|a260"), ~
         case_when(
         box == 7 & row == "C" & column == 7 ~ NA_real_,
         TRUE ~ .x)
         ),
         extraction_date =
           if_else(is.na(extraction_id), NA_Date_, extraction_date)
         ) %>%
  group_by(extraction_date) %>%
  mutate(negative_control =
         paste(sample_id[which(grepl("NC", sample_id))], collapse = ",")) %>%
  mutate(negative_control =
         if_else(is.na(extraction_date) |
                 negative_control == "", NA_character_, negative_control))

# get those matched outside box 9

 
df_out %>%
  filter(!is.na(extraction_id)) %>%
  group_by(extraction_id) %>%
  mutate(n = n()) %>%
  filter(box != 9) %>%
  pull(extraction_id) -> e_ids_matched_outside_box_9

# remove from box 9 those extraction_ids matched outside
# and manually remove those that dont fit in with the extraction

df_out <-
  df_out %>%
  filter(!(box == 9 & extraction_id %in%
           c(e_ids_matched_outside_box_9,
             404, 430,423)))

df_out %>%
  filter(!is.na(sample_id)) ->
    df_out

  cat("\nDone!\n")
  cat("linked manifest and extraction data in df_out\n")
  cat("REMEMBER: if there are two sample ids in the manifest\n")
  cat("and only one extraction id\n")
  cat("then they are both linked to the extraction id\n")
  cat("but who knows which one the extraction id\n")
  cat("*really* links to :-(\n")
