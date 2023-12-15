source("load_and_clean_mainfest.R")

df_out %>% 
    group_by(box, row, column) %>%
    mutate(n = n()) %>%
    filter(n > 1) %>%
    View()

df_out %>%
  filter(box == 7 &
         row == "C" &
         column == 7)



  df_out %>%
    filter(extraction_id == 56
    )

  
  
  df_final %>%
    filter(!is.na(pid) | grepl("NC", sample_id),
           status == "ok") %>%
    nrow
  
  # how much qc
  
  df_final %>%
    filter(!is.na(pid) | grepl("NC", sample_id),
           status == "ok") %>%
    filter(!is.na(dna_quantity) & !is.na(a260_280)) %>%
    nrow
  
  
  df_final %>%
    filter(!is.na(pid) | grepl("NC", sample_id),
           status == "ok", 
           !is.na(a260_280),
           a260_280 > 1.7) %>%
    group_by(pid) %>%
    mutate(n = n()) %>%
    group_by(n) %>%
    tally()
  
  # how many participants/samples with > 1 sample?
  
  
  df_final %>%
    filter(!is.na(pid) | grepl("NC", sample_id),
           status == "ok") %>%
    group_by(sample_id) %>%
    slice(n=1) %>%
    group_by(pid) %>%
    mutate(n = n()) %>%
    filter(n > 1) %>%
    ungroup() %>%
    summarise(
      n_samples = n(),
      n_pid = length(unique(pid))
    )
  
  # what about passing qc defined by dna quantity and ratio
  
  df_final %>%
    filter(!is.na(pid) | grepl("NC", sample_id),
           status == "ok") %>%
    filter((a260_280 > 1.7) |# & dna_quantity >= 5) |
             #           is.na(a260_280) | is.na(dna_quantity) |
             grepl("NC", sample_id)) %>%
    group_by(sample_id) %>%
    arrange(sample_id, dna_quantity) %>%
    slice(n=1) %>%
    group_by(pid) %>%
    mutate(n = n()) %>%
    filter(n > 1) %>%
    ungroup() %>%
    summarise(
      n_samples = n(),
      n_pid = length(unique(pid))
    )
  
  # breakdown by arm
  
  df_final %>%
    filter(!is.na(pid) | grepl("NC", sample_id),
           status == "ok") %>%
    group_by(sample_id) %>%
    arrange(sample_id, dna_quantity) %>%
    slice(n = 1) %>%
    group_by(pid) %>%
    mutate(n = n()) %>%
    filter(n > 1) %>%
    group_by(pid, arm) %>%
    tally() %>%
    group_by(n,arm) %>%
    tally() %>%
    pivot_wider(names_from = "arm",
                values_from = "nn")
  
  # how many samples per participants
  
  
  df_final %>%
    filter(!is.na(pid) | grepl("NC", sample_id),
           status == "ok") %>%
    filter((a260_280) > 1.7 | #& dna_quantity >= 5) |
             is.na(a260_280) | is.na(dna_quantity) |
             grepl("NC", sample_id)) %>%
    group_by(sample_id) %>%
    arrange(sample_id, dna_quantity) %>%
    slice(n=1) %>%
    group_by(pid) %>%
    mutate(n = n()) %>%
    filter(n > 1) %>%
    ungroup() %>%
    group_by(pid, arm) %>%
    tally() %>%
    group_by(n,arm) %>%
    tally() %>%
    pivot_wider(names_from = "arm",
                values_from = "nn")
  
  
  df_final %>%
    filter(!is.na(sample_id)) %>%
    ggplot(aes(a260_280, dna_quantity)) + 
    geom_point() +
    geom_vline(xintercept = c(1.7,2.5)) +
    geom_hline(yintercept = 5) +
    ylim(c(0,150))
  
  
  df_final %>%
    filter(!is.na(pid) | grepl("NC", sample_id),
           status == "ok") %>%
    filter((a260_280 > 1.7 & dna_quantity >= 5) |
             is.na(a260_280) | is.na(dna_quantity) |
             grepl("NC", sample_id)) %>%
    group_by(sample_id) %>%
    arrange(sample_id, dna_quantity) %>%
    slice(n=1) %>%
    #   mutate(ratio = 
    #          case_when(grepl("NC", sample_id) ~ "NC",
    #                    is.na(a260_280) ~ "Not done",
    #                    a260_280 >= 1.7 & a260_280 <= 2.3 ~ "1.7-2.3",
    #                    TRUE ~ "bad"
    #                    )) %>%
    group_by(arm,pid) %>%
    #   mutate(all_good_or_not_done = all(ratio == "1.7-2.3" | ratio == "Not done"),
    #          n = n()) %>%
    tally() %>%
    group_by(n) %>%
    tally()
  
  # how many samples with different strategeies
  
  df_final %>%
    filter(!is.na(pid) | grepl("NC", sample_id),
           status == "ok") %>%
    #   filter((a260_280 > 1.7 & a260_280 < 2.3) | grepl("NC", sample_id)) %>%
    mutate(ratio = 
             case_when(grepl("NC", sample_id) ~ "NC",
                       is.na(a260_280) ~ "Not done",
                       a260_280 >= 1.7 & a260_280 <= 2.3 ~ "1.7-2.3",
                       TRUE ~ "bad"
             )) %>%
    group_by(arm,pid) %>%
    mutate(all_good_or_not_done = all(ratio == "1.7-2.3" | ratio == "Not done")) %>%
    group_by(arm,pid,all_good_or_not_done) %>%
    tally() %>%
    ggplot(aes(n, fill = all_good_or_not_done)) +
    geom_bar() +
    facet_wrap(~ arm, scales = "free_x")
  
  
  df_final %>%
    filter(!is.na(pid) | grepl("NC", sample_id),
           status == "ok") %>%
    filter((a260_280 > 1.7 & a260_280 < 2.3) | grepl("NC", sample_id)) %>%
    group_by(sample_id) %>%
    slice(n=1) %>%
    #   mutate(ratio = 
    #          case_when(grepl("NC", sample_id) ~ "NC",
    #                    is.na(a260_280) ~ "Not done",
    #                    a260_280 >= 1.7 & a260_280 <= 2.3 ~ "1.7-2.3",
    #                    TRUE ~ "bad"
    #                    )) %>%
    group_by(arm,pid) %>%
    #   mutate(all_good_or_not_done = all(ratio == "1.7-2.3" | ratio == "Not done"),
    #          n = n()) %>%
    tally() %>%
    group_by(n) %>%
    tally()
  
  df_final %>%
    filter(!is.na(pid) | grepl("NC", sample_id),
           status == "ok") %>%
    #   filter((a260_280 > 1.7 & a260_280 < 2.4) | grepl("NC", sample_id)) %>%
    mutate(ratio =
             case_when(grepl("NC", sample_id) ~ "NC",
                       is.na(a260_280) ~ "Not done",
                       a260_280 >= 1.7 & a260_280 <= 2.3 ~ "1.7-2.3",
                       TRUE ~ "bad"
             )) %>%
    ggplot(aes(pid, visit, color = ratio)) +
    geom_point() +
    coord_flip() +
    facet_wrap(~ arm, scale = "free")
  
  df_sliced %>%
    mutate(ratio =
             case_when(grepl("NC", sample_id) ~ "NC",
                       is.na(a260_280) ~ "Not done",
                       a260_280 >= 1.7 & a260_280 <= 2.3 ~ "1.7-2.3",
                       TRUE ~ "bad"
             )) %>%
    group_by(pid, ratio) %>%
    mutate(n = n()) %>%
    group_by(n, ratio) %>%
    tally()
  
  
  df_sliced %>%
    group_by(pid) %>%
    mutate(n = n()) %>%
    group_by(n) %>%
    tally()
  