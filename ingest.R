library(tidyverse)
library(readxl)

## alka sub-clincial pheotypes
df_alka1 <- read_delim("data/alka_sort_master_comb_annovar_annotate.txt", delim = " ") 
df_alka1 <- df_alka1 %>%
  mutate(Population = if_else(CD_pheno == "CD_pheno", "Crohn's Disease",
                              if_else(CD_pheno == "UC_pheno", "Ulcerative Colitis", "wilson"))) %>%
  select(Population, PHENOTYPE, SNP, P, OR, NMISS, A1) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Analyst = "Alka") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "75% Caucasian independent samples, phenotypes from Sultan's File, ichip 1-5") %>%
  filter(P <= 0.05)

df_alka2 <- read_excel("data/sort_merge_casecontrol__ichip1-5_trace_anotated_forGreg_041618.xlsx")
df_alka2 <- df_alka2 %>%
  mutate(P <= 0.05) %>%
  mutate(Population = if_else(Phenotype == "cd", "Crohn's Disease",
                              if_else(Phenotype == "uc", "Ulcerative Colitis",
                                      if_else(Phenotype == "ibd", "IBD", "wilson")))) %>%
  mutate(PHENOTYPE = paste(Population, "vs. non-IBD Controls", sep = " ")) %>%
  select(Population, PHENOTYPE, SNP, P, OR, NMISS, A1) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Analyst = "Alka") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "Case-control with ichip1-5 and DX 2016 phenotype file from Linda and using 75% Caucasian independent. Covariates including first 4 PCs")

df_alka3 <- read_excel("data/survival_CD_surgery_ichip1-5_allindep_cauc75CD_forGreg_041618.xlsx",
                       skip = 1)
df_alka3 <- df_alka3 %>%
  filter(`Pvalue(coxph)` <= 0.05) %>%
  rename(P = `Pvalue(coxph)`, OR_Z_B = HazardRatio, NMISS = N, A1 = MinorAllele) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(PHENOTYPE = "Time to Frist Surgery") %>%
  select(Population, PHENOTYPE, SNP, P, OR_Z_B, NMISS, A1) %>%
  mutate(Analyst = "Alka") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "Time to first CD surgery associations with ichip1-5 and Sultan’s phenotype file and using 75% Caucasian independent samples. For second surgery, phenotype file from DM. Covariates include first 4 PCs")

df_alka4 <- read_excel("data/survival_CD_surgery_ichip1-5_allindep_cauc75CD_forGreg_041618.xlsx",
                       skip = 1, sheet = 2)

df_alka4 <- df_alka4 %>%
  filter(`p-value(coxph)` <= 0.05) %>%
  rename(P = `p-value(coxph)`, OR_Z_B = HazardRatio, NMISS = N, A1 = MinorAllele) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(PHENOTYPE = "Time to Second Surgery") %>%
  select(Population, PHENOTYPE, SNP, P, OR_Z_B, NMISS, A1) %>%
  mutate(Analyst = "Alka") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "Time to Second CD surgery associations with ichip1-5 and Sultan’s phenotype file and using 75% Caucasian independent samples. For second surgery, phenotype file from DM. Covariates include first 4 PCs")

df_alka <- bind_rows(df_alka1, df_alka2, df_alka3, df_alka4)

## dalin 
df_dalin_meta_ibd <- read_tsv("data/dalin/meta/dalin_meta_results_IBD.txt")
df_dalin_meta_ibd <- df_dalin_meta_ibd %>%
  filter(P_meta_fixed_all <= 0.05) %>%
  mutate(NMISS = NMISS_cedars + NMISS_iibdgc + 59957) %>% # 59957 from de lange ng 3760 s2
  select(SNP, P_meta_fixed_all, beta_meta_fixed_all, NMISS, A1_cedars) %>%
  rename(P = P_meta_fixed_all, OR_Z_B = beta_meta_fixed_all, A1 = A1_cedars) %>%
  mutate(Population = "IBD") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "Meta of Cedars, IIBDGC, De Lange paper cohorts") %>%
  mutate(PHENOTYPE = "IBD vs. Ctrl") %>%
  distinct(SNP, .keep_all = TRUE) #metas have ~60 duplicate SNPs removing


df_dalin_meta_cd <- read_tsv("data/dalin/meta/dalin_meta_results_CD.txt")
df_dalin_meta_cd <- df_dalin_meta_cd %>%
  filter(P_meta_fixed_all <= 0.05) %>%
  mutate(NMISS = NMISS_cedars + NMISS_iibdgc + 40266) %>% # 40266 from de lange ng 3760 s2
  select(SNP, P_meta_fixed_all, beta_meta_fixed_all, NMISS, A1_cedars) %>%
  rename(P = P_meta_fixed_all, OR_Z_B = beta_meta_fixed_all, A1 = A1_cedars) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "Meta of Cedars, IIBDGC, De Lange paper cohorts") %>%
  mutate(PHENOTYPE = "CD vs. Ctrl") %>%
  distinct(SNP, .keep_all = TRUE) #metas have ~60 duplicate SNPs removing

df_dalin_meta_uc <- read_tsv("data/dalin/meta/dalin_meta_results_UC.txt")
df_dalin_meta_uc <- df_dalin_meta_uc %>%
  filter(P_meta_fixed_all <= 0.05) %>%
  mutate(NMISS = NMISS_cedars + NMISS_iibdgc + 45975) %>% # 45975 from de lange ng 3760 s2
  select(SNP, P_meta_fixed_all, beta_meta_fixed_all, NMISS, A1_cedars) %>%
  rename(P = P_meta_fixed_all, OR_Z_B = beta_meta_fixed_all, A1 = A1_cedars) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "Meta of Cedars, IIBDGC, De Lange paper cohorts") %>%
  mutate(PHENOTYPE = "UC vs. Ctrl") %>%
  distinct(SNP, .keep_all = TRUE) #metas have ~60 duplicate SNPs removing

df_anca_cd <- read_tsv("data/dalin/serology/anca_cd.txt")
df_anca_cd <- df_anca_cd %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "Anca CD vs. Anca Ctrl") 

df_asca_cd <- read_tsv("data/dalin/serology/asca_cd.txt")
df_asca_cd <- df_asca_cd %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "Asca CD vs. Asca Ctrl")

df_cbir_cd <- read_tsv("data/dalin/serology/cbir_cd.txt")
df_cbir_cd <- df_cbir_cd %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "Cbir CD vs. Cbir Ctrl")

df_i2_cd <- read_tsv("data/dalin/serology/i2_cd.txt")
df_i2_cd <- df_i2_cd %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "i2 CD vs. i2 Ctrl")

df_iga.asca_cd <- read_tsv("data/dalin/serology/iga.asca_cd.txt")
df_iga.asca_cd <- df_iga.asca_cd %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "IGA Asca CD vs. IGA Asca Ctrl")

df_igg.asca_cd <- read_tsv("data/dalin/serology/igg.asca_cd.txt")
df_igg.asca_cd <- df_igg.asca_cd %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "IGG Asca CD vs. IGA Asca Ctrl")

df_ompc_cd <- read_tsv("data/dalin/serology/ompc_cd.txt")
df_ompc_cd <- df_ompc_cd %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "OMPC CD vs. OMPC Ctrl")

df_anca_uc <- read_tsv("data/dalin/serology/anca_uc.txt")
df_anca_uc <- df_anca_uc %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "Anca UC vs. Anca Ctrl")

df_asca_uc <- read_tsv("data/dalin/serology/asca_uc.txt")
df_asca_uc <- df_asca_uc %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "Asca UC vs. Asca Ctrl")

df_cbir_uc <- read_tsv("data/dalin/serology/cbir_uc.txt")
df_cbir_uc <- df_cbir_uc %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "Cbir UC vs. Cbir Ctrl")

df_i2_uc <- read_tsv("data/dalin/serology/i2_uc.txt")
df_i2_uc <- df_i2_uc %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "i2 UC vs. i2 Ctrl")

df_iga.asca_uc <- read_tsv("data/dalin/serology/iga.asca_uc.txt")
df_iga.asca_uc <- df_iga.asca_uc %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "IGA Asca UC vs. IGA Asca Ctrl")

df_igg.asca_uc <- read_tsv("data/dalin/serology/igg.asca_uc.txt")
df_igg.asca_uc <- df_igg.asca_uc %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "IGG Asca UC vs. IGA Asca Ctrl")

df_ompc_uc <- read_tsv("data/dalin/serology/ompc_uc.txt")
df_ompc_uc <- df_ompc_uc %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = n.case + n.ctrl) %>%
  select(-n.case, -n.ctrl) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Dalin") %>%
  mutate(Year = "2017") %>%
  mutate(Notes = "iChip") %>%
  mutate(PHENOTYPE = "OMPC UC vs. OMPC Ctrl")

df_dalin <- bind_rows(df_dalin_meta_cd, df_dalin_meta_ibd, df_dalin_meta_uc, 
                      df_i2_cd, df_i2_uc, df_igg.asca_uc, df_iga.asca_uc,
                      df_igg.asca_cd, df_ompc_uc, df_anca_uc, df_ompc_cd,
                      df_asca_cd, df_cbir_cd, df_cbir_uc, df_anca_cd,
                      df_asca_uc)

## talin
df_talin_1 <- read_table2("data/talin/aTNF/IBDichip12345TOP_TNF_CD_PNRvSNRDR_2PC.assoc.logistic.txt") #error is okay adding additional column 
df_talin_1 <- df_talin_1 %>%
  filter(P <= 0.05) %>%
  filter(TEST == "ADD") %>%
  select(SNP, P, OR, NMISS, A1) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Cedars ichip1-5 CD anti-TNF all races but mostly Caucasian, adjusted for 2PC, smoking, disease location, combination therapy") %>%
  mutate(PHENOTYPE = "anti-TNF primary non-response vs primary response")

df_talin_2 <- read_table2("data/talin/aTNF/IBDichip12345TOP_TNF_CD_coxph_2PCclin_autoR2.txt")
df_talin_2 <- df_talin_2 %>%
  rename(SNP = snp, P = pvalue_snp, OR_Z_B = zscore_snp, NMISS = number_subjects, A1  = allele) %>%
  filter(P <= 0.05) %>%
  select(SNP, P, OR_Z_B, NMISS, A1) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Cedars ichip1-5 CD anti-TNF all races but mostly Caucasian, adjusted for 2PC, anca level, family history IBD") %>%
  mutate(PHENOTYPE = "anti-TNF time to loss of response") 

df_talin_3 <- read_table2("data/talin/aTNF/IBDichip12345TOP_TNF_UC_2pc.PNR2_vs_SNRDR1.assoc.logistic.txt")
df_talin_3 <- df_talin_3 %>%
  filter(P <= 0.05) %>%
  filter(TEST == "ADD") %>%
  select(SNP, P, OR, NMISS, A1) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Ulcerative Colitits") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Cedars ichip1-5 UC anti-TNF all races but mostly Caucasian, adjusted for 2PC") %>%
  mutate(PHENOTYPE = "anti-TNF primary non-response vs primary response")

df_talin_4 <- read_table2("data/talin/aTNF/IBDichip12345TOP_TNF_UC_coxph_2PCclin_autoR2.txt")
df_talin_4 <- df_talin_4 %>%
  select(snp, pvalue_snp, zscore_snp, number_subjects, allele) %>%
  rename(SNP = snp, P = pvalue_snp, OR_Z_B = zscore_snp, NMISS = number_subjects, A1  = allele) %>%
  filter(P <= 0.05) %>%
  select(SNP, P, OR_Z_B, NMISS, A1) %>%
  mutate(Population = "Ulcerative Colitits") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Cedars ichip1-5 UC anti-TNF all races but mostly Caucasian, adjusted for 2PC, ompc positivity, anca positivity") %>%
  mutate(PHENOTYPE = "anti-TNF time to loss of response") 

df_talin_5 <- read_table2("data/talin/aTNF/IBDichip12345TOP_TNF_IBD_PNRvPR_2PCAgeDx.assoc.logistic.txt")
df_talin_5 <- df_talin_5 %>%
  filter(P <= 0.05) %>%
  filter(TEST == "ADD") %>%
  select(SNP, P, OR, NMISS, A1) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "IBD") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Cedars ichip1-5 all IBD anti-TNF all races but mostly Caucasian, adjusted for 2PC, age dx") %>%
  mutate(PHENOTYPE = "anti-TNF primary non-response vs primary response")

df_talin_6 <- read_table2("data/talin/aTNF/IBDichip12345TOP_TNF_IBD_coxph_2PCFamHx.autoR.txt")
df_talin_6 <- df_talin_6 %>%
  select(snp, pvalue_snp, zscore_snp, number_subjects, allele) %>%
  rename(SNP = snp, P = pvalue_snp, OR_Z_B = zscore_snp, NMISS = number_subjects, A1  = allele) %>%
  filter(P <= 0.05) %>%
  select(SNP, P, OR_Z_B, NMISS, A1) %>%
  mutate(Population = "IBD") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Cedars ichip1-5 all IBD anti-TNF all races but mostly Caucasian, adjusted for 2PC, fam dx") %>%
  mutate(PHENOTYPE = "anti-TNF time to loss of response") 

df_talin_7 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_linearD0.assoc.linear.txt")
df_talin_7 <- df_talin_7 %>%
  filter(P <= 0.05) %>%
  filter(TEST == "ADD") %>%
  select(SNP, P, BETA, NMISS, A1) %>%
  rename(OR_Z_B = BETA) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "Cedars ichip1-6 + Stappenbeck samples, Caucasian only, adjusted for 2PC; with corresponding permuted results for selected top snps [permuted due to highly skewed phenotype; need EMP1 pvalues from permuted") %>%
  mutate(PHENOTYPE = "Paneth-D0 phenotype")

df1 <- read_table2("data/talin/Paneth/IBDichip123456aThadTIOP_panethCauc_D0_1M.assoc.linear.mperm.txt")
df10 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D0_10M.assoc.linear.mperm.txt")
df100 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D0_100M.assoc.linear.mperm.txt")
permuts <- bind_rows(df1, df10, df100)

df_talin_7 <- left_join(df_talin_7, permuts, by = c("SNP" = "SNP"))
df_talin_7 <- df_talin_7 %>%
  mutate(orig_P = ifelse(!is.na(EMP1), P, NA)) %>%
  mutate(P = if_else(!is.na(EMP1), EMP1, P)) %>%
  select(-CHR, -EMP1, -EMP2, -X5) 

df_talin_8 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_linear.D1.assoc.linear.txt")
df_talin_8 <- df_talin_8 %>%
  filter(P <= 0.05) %>%
  filter(TEST == "ADD") %>%
  select(SNP, P, BETA, NMISS, A1) %>%
  rename(OR_Z_B = BETA) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "Cedars ichip1-6 + Stappenbeck samples, Caucasian only, adjusted for 2PC; with corresponding permuted results for selected top snps [permuted due to highly skewed phenotype; need EMP1 pvalues from permuted") %>%
  mutate(PHENOTYPE = "Paneth-D1 phenotype")

df1 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D1_10M.assoc.linear.mperm.txt")
df100 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D1_100M.assoc.linear.mperm.txt")
permuts <- bind_rows(df1, df100)

df_talin_8 <- left_join(df_talin_8, permuts, by = c("SNP" = "SNP"))
df_talin_8 <- df_talin_8 %>%
  mutate(orig_P = ifelse(!is.na(EMP1), P, NA)) %>%
  mutate(P = if_else(!is.na(EMP1), EMP1, P)) %>%
  select(-CHR, -EMP1, -EMP2, -X5) 

df_talin_9 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_linear.D2.assoc.linear.txt")
df_talin_9 <- df_talin_9 %>%
  filter(P <= 0.05) %>%
  filter(TEST == "ADD") %>%
  select(SNP, P, BETA, NMISS, A1) %>%
  rename(OR_Z_B = BETA) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "Cedars ichip1-6 + Stappenbeck samples, Caucasian only, adjusted for 2PC; with corresponding permuted results for selected top snps [permuted due to highly skewed phenotype; need EMP1 pvalues from permuted") %>%
  mutate(PHENOTYPE = "Paneth-D2 phenotype")

df1 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D2_1M.assoc.linear.mperm.txt")
df10 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D2_10M.assoc.linear.mperm.txt")
df100 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D2_100M.assoc.linear.mperm.txt")
permuts <- bind_rows(df1, df10, df100)

df_talin_9 <- left_join(df_talin_9, permuts, by = c("SNP" = "SNP"))
df_talin_9 <- df_talin_9 %>%
  mutate(orig_P = ifelse(!is.na(EMP1), P, NA)) %>%
  mutate(P = if_else(!is.na(EMP1), EMP1, P)) %>%
  select(-CHR, -EMP1, -EMP2, -X5) 

df_talin_10 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_linear.D3.assoc.linear.txt")
df_talin_10 <- df_talin_10 %>%
  filter(P <= 0.05) %>%
  filter(TEST == "ADD") %>%
  select(SNP, P, BETA, NMISS, A1) %>%
  rename(OR_Z_B = BETA) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "Cedars ichip1-6 + Stappenbeck samples, Caucasian only, adjusted for 2PC; with corresponding permuted results for selected top snps [permuted due to highly skewed phenotype; need EMP1 pvalues from permuted") %>%
  mutate(PHENOTYPE = "Paneth-D3 phenotype")

df10 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D3_10M.assoc.linear.mperm.txt")
df100 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D3_100M.assoc.linear.mperm.txt")
permuts <- bind_rows(df10, df100)

df_talin_10 <- left_join(df_talin_10, permuts, by = c("SNP" = "SNP"))
df_talin_10 <- df_talin_10 %>%
  mutate(orig_P = ifelse(!is.na(EMP1), P, NA)) %>%
  mutate(P = if_else(!is.na(EMP1), EMP1, P)) %>%
  select(-CHR, -EMP1, -EMP2, -X5) 

df_talin_11 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_linear.D4.assoc.linear.txt")
df_talin_11 <- df_talin_11 %>%
  filter(P <= 0.05) %>%
  filter(TEST == "ADD") %>%
  select(SNP, P, BETA, NMISS, A1) %>%
  rename(OR_Z_B = BETA) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "Cedars ichip1-6 + Stappenbeck samples, Caucasian only, adjusted for 2PC; with corresponding permuted results for selected top snps [permuted due to highly skewed phenotype; need EMP1 pvalues from permuted") %>%
  mutate(PHENOTYPE = "Paneth-D4 phenotype")

df10 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D4_10M.assoc.linear.mperm.txt")
df100 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D4_100M.assoc.linear.mperm.txt")
permuts <- bind_rows(df10, df100)

df_talin_11 <- left_join(df_talin_11, permuts, by = c("SNP" = "SNP"))
df_talin_11 <- df_talin_11 %>%
  mutate(orig_P = ifelse(!is.na(EMP1), P, NA)) %>%
  mutate(P = if_else(!is.na(EMP1), EMP1, P)) %>%
  select(-CHR, -EMP1, -EMP2, -X5) 

df_talin_12 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_linear.D5.assoc.linear.txt")
df_talin_12 <- df_talin_12 %>%
  filter(P <= 0.05) %>%
  filter(TEST == "ADD") %>%
  select(SNP, P, BETA, NMISS, A1) %>%
  rename(OR_Z_B = BETA) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "Cedars ichip1-6 + Stappenbeck samples, Caucasian only, adjusted for 2PC; with corresponding permuted results for selected top snps [permuted due to highly skewed phenotype; need EMP1 pvalues from permuted") %>%
  mutate(PHENOTYPE = "Paneth-D5 phenotype")

df10 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D5_10M.assoc.linear.mperm.txt")
df50 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D5_50Ma.assoc.linear.mperm.txt")
df100 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D5_100Mb.assoc.linear.mperm.txt")
permuts <- bind_rows(df10, df50, df100)

df_talin_12 <- left_join(df_talin_12, permuts, by = c("SNP" = "SNP"))
df_talin_12 <- df_talin_12 %>%
  mutate(orig_P = ifelse(!is.na(EMP1), P, NA)) %>%
  mutate(P = if_else(!is.na(EMP1), EMP1, P)) %>%
  select(-CHR, -EMP1, -EMP2, -X5) 

df_talin_13 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_linear.D1234.assoc.linear.txt")
df_talin_13 <- df_talin_13 %>%
  filter(P <= 0.05) %>%
  filter(TEST == "ADD") %>%
  select(SNP, P, BETA, NMISS, A1) %>%
  rename(OR_Z_B = BETA) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "Cedars ichip1-6 + Stappenbeck samples, Caucasian only, adjusted for 2PC; with corresponding permuted results for selected top snps [permuted due to highly skewed phenotype; need EMP1 pvalues from permuted") %>%
  mutate(PHENOTYPE = "Paneth-D1234 phenotype")

df1 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D1234_1M.assoc.linear.mperm.txt")
df10 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D1234_10M.assoc.linear.mperm.txt")
df100 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_D1234_100M.assoc.linear.mperm.txt")
permuts <- bind_rows(df1, df10, df100)

df_talin_13 <- left_join(df_talin_13, permuts, by = c("SNP" = "SNP"))
df_talin_13 <- df_talin_13 %>%
  mutate(orig_P = ifelse(!is.na(EMP1), P, NA)) %>%
  mutate(P = if_else(!is.na(EMP1), EMP1, P)) %>%
  select(-CHR, -EMP1, -EMP2, -X5) 

df_talin_14 <- read_table2("data/talin/Paneth/IBDichip123456aThadTOP_panethCauc_highlow.assoc.logistic.txt")
df_talin_14 <- df_talin_14 %>%
  filter(P <= 0.05) %>%
  filter(TEST == "ADD") %>%
  select(SNP, P, OR, NMISS, A1) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2016 - 2017") %>%
  mutate(Notes = "Cedars ichip1-6 + Stappenbeck samples, Caucasian only, adjusted for 2PC") %>%
  mutate(PHENOTYPE = "Paneth high (>20%) vs low (<20%) abnormal paneth phenotype")

df_talin_15 <- read_table2("data/talin/MRUC/meta060415_mrucVnon_logistic1.txt")
df_talin_15 <- df_talin_15 %>%
  rename(SNP = MarkerName, P = `P-value`, OR_Z_B = Effect, A1 = Allele1) %>%
  select(SNP, P, OR_Z_B, A1) %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = 10060 + 962) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Meta-analysis of Cedars ichip1-4 and international IIBDGC ichip, Caucasian only") %>%
  mutate(PHENOTYPE = "MRUC vs non-MRUC")

df_talin_16 <- read_table2("data/talin/MRUC/meta063015_snp2hla_mrucVnon1_logistic.txt")
df_talin_16 <- df_talin_16 %>%
  rename(SNP = MarkerName, P = `P-value`, OR_Z_B = Effect, A1 = Allele1) %>%
  select(SNP, P, OR_Z_B, A1) %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = 10060 + 962) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Meta-analysis of Cedars iChip1-4 and international IIBDGC ichip using SNP2HLA (MHC/HLA) imputed, Caucasian only") %>%
  mutate(PHENOTYPE = "MRUC vs non-MRUC")

df_talin_17 <- read_table2("data/talin/MRUC/meta061715_mrucVnon_coxphALLyrs1.txt")
df_talin_17 <- df_talin_17 %>%
  rename(SNP = MarkerName, P = `P-value`, OR_Z_B = Effect, A1 = Allele1) %>%
  select(SNP, P, OR_Z_B, A1) %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = 947 + 5379) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Meta-analysis of Cedars iChip1-4 and international IIBDGC ichip, Caucasian only") %>%
  mutate(PHENOTYPE = "MRUC time to colectomy")
  
df_talin_18 <- read_table2("data/talin/MRUC/meta062915_mruc_snp2hla_coxphALLyrs1.txt")
df_talin_18 <- df_talin_18 %>%
  rename(SNP = MarkerName, P = `P-value`, OR_Z_B = Effect, A1 = Allele1) %>%
  select(SNP, P, OR_Z_B, A1) %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = 947 + 5379) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Meta-analysis of Cedars iChip1-4 and international IIBDGC ichip using SNP2HLA (MHC/HLA) imputed, Caucasian only") %>%
  mutate(PHENOTYPE = "MRUC time to colectomy")

df_talin_19 <- read_table2("data/talin/MRUC/meta060415_mrucVnon_coxph5yrs_b1.txt")
df_talin_19 <- df_talin_19 %>%
  rename(SNP = MarkerName, P = `P-value`, OR_Z_B = Effect, A1 = Allele1) %>%
  select(SNP, P, OR_Z_B, A1) %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = 596 + 3335) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Meta-analysis of Cedars iChip1-4 and international IIBDGC ichip, Caucasian only") %>%
  mutate(PHENOTYPE = "MRUC Time to Colectomy <60 months vs >60 months")

df_talin_20 <- read_table2("data/talin/MRUC/meta062515_mrucVnon_snp2hla_coxph5yrs1.txt")
df_talin_20 <- df_talin_20 %>%
  rename(SNP = MarkerName, P = `P-value`, OR_Z_B = Effect, A1 = Allele1) %>%
  select(SNP, P, OR_Z_B, A1) %>%
  filter(P <= 0.05) %>%
  mutate(Population = "Ulcerative Colitis") %>%
  mutate(NMISS = 596 + 3335) %>%
  mutate(Analyst = "Talin") %>%
  mutate(Year = "2015 - 2016") %>%
  mutate(Notes = "Meta-analysis of Cedars iChip1-4 and international IIBDGC ichip using SNP2HLA (MHC/HLA) imputed, Caucasian only") %>%
  mutate(PHENOTYPE = "MRUC Time to Colectomy <60 months vs >60 months")

df_talin <- bind_rows(df_talin_1, df_talin_2, df_talin_3, df_talin_4,
                      df_talin_5, df_talin_6, df_talin_7, df_talin_8,
                      df_talin_9, df_talin_10, df_talin_11, df_talin_12,
                      df_talin_13, df_talin_14, df_talin_15, df_talin_16,
                      df_talin_17, df_talin_18, df_talin_19, df_talin_20)

# ingest shishir
df_shishir_1 <- read_table2("data/shishir/B1vCtrl_L1_MAF.Status.assoc.logistic")
df_shishir_1 <- df_shishir_1 %>%
  filter(P <= 0.05) %>%
  select(SNP, A1, NMISS, OR, L95, U95, P) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Shishir") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "iChip 1-7, Caucasian") %>%
  mutate(PHENOTYPE = "L1 B1 vs non-IBD ctrl")

df_shishir_2 <- read_table2("data/shishir/B1vCtrl_L2_MAF.Status.assoc.logistic")
df_shishir_2 <- df_shishir_2 %>%
  filter(P <= 0.05) %>%
  select(SNP, A1, NMISS, OR, L95, U95, P) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Shishir") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "iChip 1-7, Caucasian") %>%
  mutate(PHENOTYPE = "L2 B1 vs non-IBD ctrl")

df_shishir_3 <- read_table2("data/shishir/B1vCtrl_L3_MAF.Status.assoc.logistic")
df_shishir_3 <- df_shishir_3 %>%
  filter(P <= 0.05) %>%
  select(SNP, A1, NMISS, OR, L95, U95, P) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Shishir") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "iChip 1-7, Caucasian") %>%
  mutate(PHENOTYPE = "L3 B1 vs non-IBD ctrl")

df_shishir_4 <- read_table2("data/shishir/B2a+B2bvB1_L1_MAF.Status.assoc.logistic")
df_shishir_4 <- df_shishir_4 %>%
  filter(P <= 0.05) %>%
  select(SNP, A1, NMISS, OR, L95, U95, P) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Shishir") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "iChip 1-7, Caucasian") %>%
  mutate(PHENOTYPE = "L1 B2a+B2b vs non-IBD ctrl")

df_shishir_5 <- read_table2("data/shishir/B2a+B2bvB1_L2_MAF.Status.assoc.logistic")
df_shishir_5 <- df_shishir_5 %>%
  filter(P <= 0.05) %>%
  select(SNP, A1, NMISS, OR, L95, U95, P) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Shishir") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "iChip 1-7, Caucasian") %>%
  mutate(PHENOTYPE = "L2 B2a+B2b vs non-IBD ctrl")

df_shishir_6 <- read_table2("data/shishir/B2a+B2bvB1_L3_MAF.Status.assoc.logistic")
df_shishir_6 <- df_shishir_6 %>%
  filter(P <= 0.05) %>%
  select(SNP, A1, NMISS, OR, L95, U95, P) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Shishir") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "iChip 1-7, Caucasian") %>%
  mutate(PHENOTYPE = "L3 B2a+B2b vs non-IBD ctrl")

df_shishir_7 <- read_table2("data/shishir/B2a+B2bvCtrl_L1_MAF.Status.assoc.logistic")
df_shishir_7 <- df_shishir_7 %>%
  filter(P <= 0.05) %>%
  select(SNP, A1, NMISS, OR, L95, U95, P) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Shishir") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "iChip 1-7, Caucasian") %>%
  mutate(PHENOTYPE = "L1 B2a+B2b vs B1")

df_shishir_8 <- read_table2("data/shishir/B2a+B2bvCtrl_L2_MAF.Status.assoc.logistic")
df_shishir_8 <- df_shishir_8 %>%
  filter(P <= 0.05) %>%
  select(SNP, A1, NMISS, OR, L95, U95, P) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Shishir") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "iChip 1-7, Caucasian") %>%
  mutate(PHENOTYPE = "L2 B2a+B2b vs B1")

df_shishir_9 <- read_table2("data/shishir/B2a+B2bvCtrl_L3_MAF.Status.assoc.logistic")
df_shishir_9 <- df_shishir_9 %>%
  filter(P <= 0.05) %>%
  select(SNP, A1, NMISS, OR, L95, U95, P) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Shishir") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "iChip 1-7, Caucasian") %>%
  mutate(PHENOTYPE = "L3 B2a+B2b vs B1")
  
df_shishir <- bind_rows(df_shishir_1, df_shishir_2, df_shishir_3, df_shishir_4,
                        df_shishir_5, df_shishir_6, df_shishir_7, df_shishir_8,
                        df_shishir_9)

# marcy associations
df_marcy_1 <- read_xlsx("data/marcy/International_pCD+ver.pCD-.xlsx")
df_marcy_1 <- df_marcy_1 %>%
  filter(P <= 0.05) %>%
  mutate(NMISS = 2974 + 7764) %>%
  select(SNP, A1, NMISS, OR, P) %>%
  rename(OR_Z_B = OR) %>%
  mutate(Population = "Crohn's Disease") %>%
  mutate(Analyst = "Marcy") %>%
  mutate(Year = "2018") %>%
  mutate(Notes = "International Cohort") %>%
  mutate(PHENOTYPE = "Perianal CD+ vs.  Perianal CD- ")


# df all
df_all <- bind_rows(df_alka, df_dalin, df_talin, df_shishir, df_marcy_1)

# ensure A1 always capital
df_all <- df_all %>%
  mutate(A1 = toupper(A1))


# 1354 SNPs without an A1 from Dalin's IBD meta 
# 7117 SNPS without NMISS from Talin's MRUCs and Dalins metas 

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
zero_snps <- df_all %>%
  filter(!SNP %in% ichip_v1_v2_anno$Illumina_ichip_ID) %>%
  filter(SNP %in% snp_list$`With zero instead of dash`) %>%
  select(SNP) %>%
  distinct() 

zero_snps <- left_join(zero_snps, snp_list, by = c("SNP" = 'With zero instead of dash'))

df_all <- left_join(df_all, zero_snps, by = "SNP")

df_all <- df_all %>%
  mutate(SNP = if_else(!is.na(Original), Original, SNP)) %>%
  select(-Original)

df_all %>%
  filter(!SNP %in% ichip_v1_v2_anno$Illumina_ichip_ID) %>%
  filter(!SNP %in% ichip_v1_v2_anno$avsnp147) %>%
  group_by(Analyst) %>%
  summarise(n = n())

# the remaining SNPS primarily come from Talin's HLA analaysis

df_annotated1 <- df_all %>%
  filter(SNP %in% ichip_v1_v2_anno$Illumina_ichip_ID)
df_annotated1 <- left_join(df_annotated1, ichip_v1_v2_anno, by = c("SNP" = "Illumina_ichip_ID"))

df_annotated2 <- df_all %>%
  filter(!SNP %in% ichip_v1_v2_anno$Illumina_ichip_ID)
df_annotated2 <- left_join(df_annotated2, ichip_v1_v2_anno, by = c("SNP" = "avsnp147")) 
df_annotated2 <- df_annotated2 %>%
  select(-Illumina_ichip_ID) %>%
  filter(!is.na(Gene.knownGene))

df_annotated3 <- df_all %>%
  filter(!SNP %in% ichip_v1_v2_anno$Illumina_ichip_ID & !SNP %in% ichip_v1_v2_anno$avsnp147) %>%
  mutate(Gene.knownGene = ifelse(Analyst == "Talin" & str_detect(Notes, "MHC/HLA"), "MHC/HLA", NA)) %>%
  mutate(Gene.refGene = ifelse(Analyst == "Talin" & str_detect(Notes, "MHC/HLA"), "MHC/HLA", NA)) %>%
  mutate(Func.knownGene = ifelse(Analyst == "Talin" & str_detect(Notes, "MHC/HLA"), "MHC/HLA Perm", NA)) %>%
  mutate(Chr = ifelse(Analyst == "Talin" & str_detect(Notes, "MHC/HLA"), 6, NA))

df_all_annotated <- bind_rows(df_annotated1, df_annotated2, df_annotated3)


## add sb 139 dataset as an additional annotation
sb_139_cis <- read_excel("data/cis_eqtl_sb139_eigen_full_annot.xlsx")
sb_139_cis <- sb_139_cis %>%
  group_by(Illumina_ID) %>% 
  summarise(cis_eGENE = paste(cis_eGENE, collapse = ", "),
            eqtl_beta = paste(beta, collapse = ", "),
            eqtl_p = paste(`p-value`, collapse = ", "), 
            priority = paste(`Match_PRIORITY_cis-eqtl`, collapse = ", "))


zero_snps_eqtl <- sb_139_cis %>%
  filter(!Illumina_ID %in% ichip_v1_v2_anno$Illumina_ichip_ID) %>%
  filter(Illumina_ID %in% snp_list$`With zero instead of dash`) %>%
  select(Illumina_ID) %>%
  distinct() 

zero_snps_eqtl <- left_join(zero_snps_eqtl, snp_list, by = c("Illumina_ID" = 'With zero instead of dash'))

sb_139_cis  <- left_join(sb_139_cis , zero_snps_eqtl, by = "Illumina_ID")

sb_139_cis <- sb_139_cis %>%
  mutate(Illumina_ID = if_else(!is.na(Original), Original, Illumina_ID)) %>%
  select(-Original)

df_all_annotated <- left_join(df_all_annotated, sb_139_cis, by = c("SNP" = "Illumina_ID"))

# the filter criteria needs to be compelte for all markers otherwise even searching the exact SNP 
# will not return the desired result. e.g. if rsxyz does not have a Func location it will not be visible
# replace unknown Func with "UNK"
# replace unknown Chr with 999
# replace unknown start with 0

df_all_annotated <- df_all_annotated %>%
  mutate(Func.knownGene = if_else(is.na(Func.knownGene), "UNK", Func.knownGene)) %>%
  mutate(Chr = if_else(is.na(Chr), 999, Chr)) %>%
  mutate(Start = if_else(is.na(Start), 0, Start))

  
df_all_annotated %>%
  filter(P != 0) %>% # 40  makres from Dalin have P = 0 need to fix removing for now
  write_csv("df_all_annotated.csv")

sb_139_cis %>%
  filter(!Illumina_ID %in% df_all_annotated$SNP) %>%
  View()

