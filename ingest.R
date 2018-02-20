library(tidyverse)



df_alka <- read_delim("data/sort_master_comb_annovar_annotate.txt", delim = " ") %>%
  select(PHENOTYPE, SNP, P, OR, NMISS, A1) %>%
  mutate(Analyst = "Alka") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "75% Caucasian independent samples, phenotypes from Sultan's File, ichip 1-5") %>%
  filter(P <= 0.05)




