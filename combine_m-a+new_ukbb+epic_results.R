setwd("~/Dropbox/PRS_clean_run/PRS_meta-anlysis/2022march/")
`%nin%` = Negate(`%in%`)
# first, combine the files in the same way as for meta-analysis (the chunk of code below is from meta-analysis script)

df_epic = fread("~/Dropbox/PRS_clean_run/PRS_meta-anlysis/2022march/epic_prs_for_m-a.csv")
df_ukbb_cancer_prs = fread("~/Dropbox/PRS_clean_run/4_results/2022march_reanalysis/cancer_prs_ukbb_2022-03-23.csv")
df_ukbb_t2d_groups_prs = read.csv("~/Dropbox/PRS_clean_run/4_results/2022march_reanalysis/t2d_and_groups_prs_ukbb_2022-03-23.csv")

all(colnames(df_ukbb_cancer_prs) == colnames(df_ukbb_t2d_groups_prs))
df_ukbb = rbind(df_ukbb_cancer_prs, df_ukbb_t2d_groups_prs)

# which columns from epic are not in UKBB files
colnames(df_epic)[colnames(df_epic) %nin% colnames(df_ukbb)]
df_ukbb$Cohort = "UKBB"

colnames(df_ukbb)[colnames(df_ukbb) %nin% colnames(df_epic)]
df_ukbb[, colnames(df_ukbb)[colnames(df_ukbb) %nin% colnames(df_epic)]] = list(NULL)

df = rbind(df_ukbb, df_epic)

df %>% select(where(is.character)) %>% apply(., 2, function(x){levels(as.factor(x))})
df$Weighting[df$Weighting == "Betas"] = "Weighted"

df = df %>% mutate(across(is.character, factor))

levels(df$Outcome)[levels(df$Outcome)=="bc2"] <- "Breast Cancer"
levels(df$Outcome)[levels(df$Outcome)=="crc"] <- "Colorectal Cancer"
levels(df$Outcome)[levels(df$Outcome)=="panc"] <- "Pancreatic Cancer"
levels(df$Outcome)[levels(df$Outcome)=="prc"] <- "Prostate Cancer"
levels(df$Outcome)[levels(df$Outcome)=="T2D.all"] <- "T2D"

levels(df$SNPs)[levels(df$SNPs)=="brc"] <- "Breast cancer"
levels(df$SNPs)[levels(df$SNPs)=="crc"] <- "Colorectal cancer"
levels(df$SNPs)[levels(df$SNPs)=="panc"] = "Pancreatic cancer"
levels(df$SNPs)[levels(df$SNPs)=="prc"] = "Prostate cancer"
levels(df$SNPs)[levels(df$SNPs)=="group1_464"] = "1a. Higher Adiposity and decreased Sex Hormones (N = 464)​"
levels(df$SNPs)[levels(df$SNPs)=="group5_36"] = "1b. Higher Adiposity, Insulin Resistance, Lower Age at Menarche and decreased Sex Hormones (N = 36)"
levels(df$SNPs)[levels(df$SNPs)=="group4_78"] = "2. Reduced b-cell function (N=78)"
levels(df$SNPs)[levels(df$SNPs)=="group2_409"] = "3. Decreased lipids and Blood Pressure (N=409)"
levels(df$SNPs)[levels(df$SNPs)=="group3_96"] = "4. Metabolic Syndrome and Insulin Resistance (N=96)"
levels(df$SNPs)[levels(df$SNPs)=="t2d"] = "T2D"

df %>% select(where(is.factor)) %>% apply(., 2, function(x){levels(as.factor(x))})

df = df %>% filter(Adjustment == "6 PCs", Weighting == "Weighted") %>% select(-N_cases, -N_controls, -Weighting, -Adjustment)
colnames(df)

df_epic_harmonised = df %>% 
  filter(Cohort == "EPIC") %>% 
  select(-Cohort) %>%
  rename(epic_OR = OR,
         epic_lower_CI = lower_CI,
         epic_upper_CI = upper_CI,
         epic_P_value = P_value)

df_ukbb_harmonised = df %>% 
  filter(Cohort == "UKBB") %>% 
  select(-Cohort) %>%
  rename(ukbb_OR = OR,
         ukbb_lower_CI = lower_CI,
         ukbb_upper_CI = upper_CI,
         ukbb_P_value = P_value)

ma = fread("PRS_meta_22-03-31.txt")
colnames(ma)
ma = ma %>% rename(SNPs = V1,
                   Outcome = PRS_Outcome)

# make "wide"
ma_epic_ukbb = full_join(ma, df_epic_harmonised) %>% full_join(., df_ukbb_harmonised)
write.table(ma_epic_ukbb, "all_prs_m-a+new_ukbb+epic_22-04-03.csv ", row.names = F, quote = T, sep = ",")
  