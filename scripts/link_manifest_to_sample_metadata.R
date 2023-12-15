library(blantyreESBL)
library(here)
source(here("scripts/load_and_clean_mainfest.R"))

df_out %>%
  filter(!is.na(sample_id)) %>%
  left_join(
            btESBL_stoolESBL,
            by = c("sample_id" = "lab_id")
            ) %>% 
  filter(!(is.na(pid) & !grepl("NC", sample_id))) ->
              df_final

# sort out unlinked lab ids

df_final %>% 
  filter(is.na(pid),
         !grepl("NC", sample_id)) %>% 
           View()


# Problems to solve:
#0) decide what are sensible cutoffs for QC
#1) where there are multiple extractions of the same sample, choose one
#2) decide what to do about samples without QC data
#3) decide on a strategy of which samples to select for sequencing

# Problem 0

# ratio > 1.7, dna > 5 ????
bind_rows(
  df_final %>%
    filter((a260_280 > 1.7 & dna_quantity >5) |
             is.na(a260_280))%>% 
    group_by(sample_id) %>% 
    slice(n = 1) %>% 
    group_by(pid, arm) %>% 
    tally() %>% 
    group_by(n, arm) %>% 
    tally() %>% 
    mutate(type = "inc. missing"),
df_final %>%
  filter((a260_280 > 1.7 & dna_quantity >5))%>% 
  group_by(sample_id) %>% 
  slice(n = 1) %>% 
  group_by(pid, arm) %>% 
  tally() %>% 
  group_by(n, arm) %>% 
  tally() %>% 
  mutate(type = "only passing QC") 
) %>% 
  ggplot(aes(n, nn, group= type, fill = type)) +
  geom_col(position = "dodge") + 
  facet_wrap(~ arm)

df_final %>%
  # filter((a260_280 > 1.7 & dna_quantity >5) |
  #          is.na(a260_280))%>% 
  mutate(qc = 
           case_when(a260_280 > 1.7 & dna_quantity >5 ~ "01good",
                     is.na(a260_280) ~ "03notdone",
                     TRUE ~ "02bad")
  ) %>% 
  group_by(sample_id) %>% 
  slice(n = 1) %>% 
  filter(!is.na(pid)) %>% 
  group_by(pid) %>% 
  mutate(n = n()) %>% 
  filter(n > 1) %>% 
  filter(qc != "02bad") %>% 
  ggplot(aes(visit, fct_reorder(pid,n), fill = qc)) +
  geom_tile() +
  facet_wrap(~ arm, scales = "free_y", ncol = 3)

# I think rectal swabs with no poo are part of the problem


df_final %>%
  filter(!is.na(sample_id)) %>%
  ggplot(aes(a260_280, dna_quantity, color = sample_type), alpha = 0.5) + 
  geom_point() +
  geom_vline(xintercept = c(1.7,2.5)) +
  geom_hline(yintercept = 5) +
  ylim(c(0,150)) +
  facet_wrap(~ visit)


df_final %>%
  # filter((a260_280 > 1.7 & dna_quantity >5) |
  #          is.na(a260_280))%>% 
  mutate(qc = 
           case_when(a260_280 > 1.7 & dna_quantity >5 ~ "01good",
                     is.na(a260_280) ~ "03notdone",
                     TRUE ~ "02bad")) %>% 
  filter(!is.na(pid)) %>% 
  ggplot(aes(visit, fill= qc)) +
  geom_bar() +
  facet_grid(arm ~ sample_type)

df_final %>%
  # filter((a260_280 > 1.7 & dna_quantity >5) |
  #          is.na(a260_280))%>% 
  mutate(qc = 
           case_when(a260_280 > 1.7 & dna_quantity >5 ~ "01good",
                     is.na(a260_280) ~ "03notdone",
                     TRUE ~ "02bad")) %>% 
  filter(!is.na(sample_type)) %>% 
  ungroup() %>% 
  count(sample_type, qc) %>% 
  group_by(sample_type) %>% 
  mutate(prop = n / sum(n)) %>% 
  filter(!grepl("notdone", qc)) %>% 
  pivot_wider(id_cols = sample_type,
              names_from = qc,
              values_from = n) %>%
  ungroup() %>% 
  select(-sample_type) %>% 
  fisher.test()

# supported by a cheeky fisher test.
# This at least suggests that the QC is worth paying attention to, I suppose

# SO - what samples to pick, aiming for ~ 450

df_final %>%
  filter((a260_280 > 1.7 & dna_quantity >5))%>% 
  group_by(sample_id) %>% 
  slice(n = 1) %>% 
  group_by(pid, arm) %>% 
  tally() %>% 
  group_by(n, arm) %>% 
  tally() %>% 
  pivot_wider(id_cols = n, names_from = arm, values_from = nn)
  

  
  df_final %>%
    filter((a260_280 > 1.7 & dna_quantity >5))%>% 
    group_by(sample_id) %>% 
    slice(n = 1) %>% 
    group_by(pid) %>% 
    mutate(n = n()) %>% 
    filter(n != 1) %>% 
    nrow()
  
  df_final %>%
    filter((a260_280 > 1.7 & dna_quantity >5) |
             is.na(a260_280)) %>% 
    group_by(sample_id) %>% 
    slice(n = 1) %>% 
    group_by(pid) %>% 
    mutate(n = n()) %>% 
    filter(n != 1) %>% 
    nrow()
  
  # how many visit 0,1,2 in arm 1 and 2?
  
  df_final %>%
      filter((a260_280 > 1.7 & dna_quantity >5) |
             is.na(a260_280)) %>% 
    group_by(sample_id) %>% 
    slice(n = 1) %>% 
    group_by(pid) %>% 
    mutate(n = n()) %>% 
    filter(n != 1) %>% 
    filter(arm %in% c(1,2), 
           visit %in% c(0,1,2),
           is.na(a260_280)) %>% 
    nrow()
  
  
  
  
# Problem 1

df_final %>% 
  group_by(sample_id) %>% 
  mutate(n = n()) %>% 
  filter(n > 1) %>% 
  arrange(pid, sample_id) %>% 
  View()

# making list for pia from here
# ultimately will need a unique seqid 
# where there are two samples linked
df_final %>%
  filter(
    !grepl("NC", sample_id),
    (a260_280 > 1.7 & dna_quantity > 5) |
           is.na(a260_280),
    status != "open") %>% 
  mutate(seqid = paste0(sample_id,"_", box, row, column)) %>% 
  # cut down to one sample
  group_by(sample_id) %>% 
  mutate(nanodrop_match_ambiguous = 
           length(unique(a260_280)) == 1 &
           length(unique(dna_quantity)) == 1) %>% 
  # filter(n > 1) %>% 
  arrange(sample_id, desc(dna_quantity)) %>% 
  # View()
  slice(n = 1) %>% 
  # now choose samples where there are more than one per participant
  group_by(pid) %>% 
  mutate(n = n()) %>% 
  filter(n != 1, 
         !all(visit >= 2)) %>% 
   filter(!(is.na(a260_280) & visit %in% c(3,4))) %>% 
  #          )
  #         ) %>%

   nrow()

  # View()
  # filter((a260_280 > 1.7 & dna_quantity >5) |
  #          is.na(a260_280))%>% 
  mutate(qc = 
           case_when(a260_280 > 1.7 & dna_quantity >5 ~ "01good",
                     is.na(a260_280) ~ "03notdone",
                     TRUE ~ "02bad")
  ) %>% 
  filter(!is.na(arm)) %>% 
  ggplot(aes(visit, fct_reorder(pid,n), fill = qc)) +
  geom_tile() +
  facet_wrap(~ arm, scales = "free_y", ncol = 3)
  
  
  group_by(pid, arm) %>% 
  tally() %>% 
  ggplot(aes(n)) +
  geom_bar() +
  facet_wrap(~ arm) +
  xlim(c(0,6))
nrow()


df_final %>%
  filter(
    grepl("NC", sample_id) |
      (a260_280 > 1.7 & dna_quantity >5) |
      is.na(a260_280)) %>% 
  mutate(seqid = paste0(sample_id,"_", box, row, column))%>% 
  # cut down to one sample
  group_by(sample_id) %>% 
  mutate(nanodrop_match_ambiguous = 
           length(unique(a260_280)) == 1 &
           length(unique(dna_quantity)) == 1) %>% 
  # filter(n > 1) %>% 
  arrange(pid, sample_id, desc(dna_quantity)) %>% 
  # View()
  slice(n = 1) %>% 
  # now choose samlpes where there are more than one per participant
  group_by(pid) %>% 
  mutate(n = n()) %>% 
  filter(n != 1) %>% # & any(visit == 0)) %>% 
  nrow()
  
  filter(!(is.na(a260_280) &
             (arm %in% c(3) |
                visit %in% c(3,4))
  )
  ) %>%
  nrow()

  
  # --- SAMPLES FOR PIA HERE
  
  df_final %>%
    filter(
      !grepl("NC", sample_id),
      (a260_280 > 1.7 & dna_quantity > 5) |
        is.na(a260_280),
      status != "open") %>% 
    mutate(seqid = paste0(sample_id,"_", box, row, column)) %>% 
    # cut down to one sample per ID
    group_by(sample_id) %>% 
    slice(n = 1) %>% 
    # now choose samples where there are more than one per participant
    # drop those where only samples are 28 days or later
    group_by(pid) %>% 
    mutate(n = n()) %>% 
    filter(n != 1, 
           !all(visit >= 2)) %>% 
    # add the un-qc'd samples - drop the 3/6 month ones
    filter(!(is.na(a260_280) & visit %in% c(3,4))) %>% 
    # now there are a couple that have on ly one sample - drop
    group_by(pid) %>%
    mutate(n = n()) %>%
    filter(n > 1) %>% 
    # add back in NCs
    bind_rows(
      df_final %>%
        filter(grepl("NC", sample_id)) %>% 
        mutate(seqid = sample_id)
    ) %>% 
    arrange(box,row,column) ->
    samples_for_deep_seq
  
  
  
  write_csv(samples_for_deep_seq,
    "data_raw/sample_management/samples_for_deep_seq.csv"
    )
  
  
  mutate(nanodrop_match_ambiguous = 
           length(unique(a260_280)) == 1 &
           length(unique(dna_quantity)) == 1) %>% 
    arrange(sample_id, desc(dna_quantity)) %>% 
    