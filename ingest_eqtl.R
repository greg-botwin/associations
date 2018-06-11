# ingest eqtl data
library(tidyverse)
library(readxl)

# SB139
list.files("data/eqtl/SB139")

sb139_trans <- read_xlsx("data/eqtl/SB139/trans_eqtl_sb139_eigen_full_annot.xlsx")
sb139_trans <- sb139_trans %>%
  mutate(Dataset = "SB 139") %>%
  mutate('Cis_Trans' = "Trans") %>%
  mutate(Population = "Crohn's Disease") %>%
  rename(eGENE = trans_eGENE, p = 'p-value') %>%
  select(SNP, eGENE, beta, p, Cis_Trans, Population, Dataset)

sb139_cis <- read_xlsx("data/eqtl/SB139/cis_eqtl_sb139_eigen_full_annot.xlsx")
sb139_cis <- sb139_cis %>%
  mutate(Dataset = "SB 139") %>%
  mutate('Cis_Trans' = "Cis") %>%
  mutate(Population = "Crohn's Disease") %>%
  rename(SNP = Illumina_ID, eGENE = cis_eGENE, p = 'p-value', priority = `Match_PRIORITY_cis-eqtl`) %>%
  select(SNP, eGENE, beta, p, priority, Cis_Trans, Population, Dataset)

eqtl_data <- bind_rows(sb139_cis, sb139_trans)

# read in alka's annotation files
ichip_v1_anno <- read_excel("data/annotate_ichip1_autosomes_noindels_nov2016_avsnp147.xlsx")
ichip_v1_anno[ichip_v1_anno == '.'] <- NA
ichip_v2_anno <- read_excel("data/ichipv2_annovar_annotation_hg19_avsnp147_basic_updated.xlsx")
ichip_v2_anno[ichip_v2_anno == '.'] <- NA

ichip_v1_anno <- ichip_v1_anno %>%
  select(Illumina_ichip_ID, Chr, Start, Ref, Alt, Func.knownGene, Gene.knownGene, avsnp147, Gene.refGene)

ichip_v2_anno <- ichip_v2_anno %>%
  select(Otherinfo_IlluminaName, Chr, Start, Ref, Alt, Func.knownGene, Gene.knownGene, avsnp147, Gene.refGene) %>%
  rename(Illumina_ichip_ID = Otherinfo_IlluminaName)

ichip_v1_v2_anno <- bind_rows(ichip_v1_anno, ichip_v2_anno)

ichip_v1_v2_anno <- ichip_v1_v2_anno %>%
  distinct(Illumina_ichip_ID, .keep_all = TRUE)

# read in snp list from talin without 0 
snp_list <- read_excel("data/ichip_snps.xlsx")

#find snps in all associations not already matched
zero_snps <- eqtl_data %>%
  filter(!SNP %in% ichip_v1_v2_anno$Illumina_ichip_ID) %>%
  filter(SNP %in% snp_list$`With zero instead of dash`) %>%
  select(SNP) %>%
  distinct() 

zero_snps <- left_join(zero_snps, snp_list, by = c("SNP" = 'With zero instead of dash'))

eqtl_data <- left_join(eqtl_data, zero_snps, by = "SNP")

eqtl_data <- eqtl_data %>%
  mutate(SNP = if_else(!is.na(Original), Original, SNP)) %>%
  select(-Original)

eqtl_data %>%
  filter(!SNP %in% ichip_v1_v2_anno$Illumina_ichip_ID) %>%
  filter(!SNP %in% ichip_v1_v2_anno$avsnp147) %>%
  View()
