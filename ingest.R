library(tidyverse)

## alka sub-clincial pheotypes
df_alka <- read_delim("data/alka_sort_master_comb_annovar_annotate.txt", delim = " ") %>%
  mutate(Population = if_else(CD_pheno == "CD_pheno", "Crohn's Disease",
                              if_else(CD_pheno == "UC_pheno", "Ulcerative Colitis", "wilson"))) %>%
  select(Population, PHENOTYPE, SNP, P, OR, NMISS, A1) %>%
  mutate(Analyst = "Alka") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "75% Caucasian independent samples, phenotypes from Sultan's File, ichip 1-5") %>%
  filter(P <= 0.05)

## dalin 
df_dalin_meta_ibd <- read_tsv("data/dalin_meta_res_with_columns_IBD.txt")
df_dalin_meta_cd <- read_tsv("data/dalin_meta_res_with_columns_CD.txt")
df <- read_tsv("data/dalin_meta_results_IBD.txt")



