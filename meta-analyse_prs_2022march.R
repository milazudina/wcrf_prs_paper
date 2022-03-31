# meta analysis
library(meta)
library(readxl)
library(tidyr)
library(dplyr)
library(data.table)

`%nin%` = Negate(`%in%`)
setwd("~/Dropbox/PRS_clean_run/PRS_meta-anlysis/2022march/")

# first, combine the files
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

df2 <- df %>% 
  filter(Adjustment == "6 PCs") %>% 
  select(SNPs, Cohort, Weighting, Outcome, OR, lower_CI, upper_CI, P_value) %>% 
  mutate(Beta = log(OR), 
         SE = (upper_CI - lower_CI)/3.92,
         CI = upper_CI - lower_CI)

wei <- df2 %>% filter(Weighting == "Weighted") 
unwei <- df2 %>% filter(Weighting == "Unweighted")

wei$studlab <- paste(wei$SNPs, wei$Outcome, sep = " | ")
unwei$studlab <- paste(unwei$SNPs, unwei$Outcome, sep = " | ")

wei_list <- split(as.data.frame(wei), f = wei$studlab)
wei_list <- wei_list[1:(length(wei_list)-1)]

meta <- data.frame(PRS_Outcome = c(rep(NA, length(wei_list))),
                   fixed_OR = c(rep(NA, length(wei_list))), 
                   fixed_lower_CI = c(rep(NA, length(wei_list))), 
                   fixed_upper_CI = c(rep(NA, length(wei_list))), 
                   fixed_P_value = c(rep(NA, length(wei_list))),
                   fixed_UKBB_weight = c(rep(NA, length(wei_list))),
                   fixed_EPIC_weight = c(rep(NA, length(wei_list))),
                   Q = c(rep(NA, length(wei_list))),
                   Q_pval = c(rep(NA, length(wei_list))),
                   random_OR = c(rep(NA, length(wei_list))), 
                   random_lower_CI = c(rep(NA, length(wei_list))), 
                   random_upper_CI = c(rep(NA, length(wei_list))), 
                   random_P_value = c(rep(NA, length(wei_list))))

for (i in 1:length(wei_list)) {
  
  m.gen <- metagen(TE = Beta,
                   seTE = SE,
                   studlab = studlab,
                   data = wei_list[[i]],
                   sm = "SMD",
                   fixed = TRUE,
                   random = TRUE,
                   method.tau = "REML",
                   hakn = TRUE,
                   title = "EPIC and UKBB PRS Meta-analysis")
  
  meta$PRS_Outcome[i] <- m.gen$studlab[1]
  meta$fixed_OR[i] <- exp(m.gen$TE.fixed)
  meta$fixed_lower_CI[i] <- exp(m.gen$lower.fixed)
  meta$fixed_upper_CI[i] <- exp(m.gen$upper.fixed)
  meta$fixed_P_value[i] <- m.gen$pval.fixed
  meta$fixed_UKBB_weight[i] <- m.gen$w.fixed[1]
  meta$fixed_EPIC_weight[i] <- m.gen$w.fixed[2]
  meta$Q[i] <- m.gen$Q
  meta$Q_pval[i] <- m.gen$pval.Q
  meta$random_OR[i] <- exp(m.gen$TE.random)
  meta$random_lower_CI[i] <- exp(m.gen$lower.random)
  meta$random_upper_CI[i] <- exp(m.gen$upper.random)
  meta$random_P_value[i] <- m.gen$pval.random
  
}

wei_meta <- meta

write.table(wei_meta, "PRS_meta_22-03-31_test.txt", row.names = F, sep = "|", quote = F)


unwei_list <- split(as.data.frame(unwei), f = unwei$studlab)
unwei_list <- unwei_list[1:(length(unwei_list)-1)]

meta <- data.frame(PRS_Outcome = c(rep(NA, length(unwei_list))),
                   fixed_OR = c(rep(NA, length(unwei_list))), 
                   fixed_lower_CI = c(rep(NA, length(unwei_list))), 
                   fixed_upper_CI = c(rep(NA, length(unwei_list))), 
                   fixed_P_value = c(rep(NA, length(unwei_list))),
                   fixed_UKBB_weight = c(rep(NA, length(unwei_list))),
                   fixed_EPIC_weight = c(rep(NA, length(unwei_list))),
                   Q = c(rep(NA, length(unwei_list))),
                   Q_pval = c(rep(NA, length(unwei_list))),
                   random_OR = c(rep(NA, length(unwei_list))), 
                   random_lower_CI = c(rep(NA, length(unwei_list))), 
                   random_upper_CI = c(rep(NA, length(unwei_list))), 
                   random_P_value = c(rep(NA, length(unwei_list))))

for (i in 1:length(unwei_list)) {
  
  m.gen <- metagen(TE = Beta,
                   seTE = SE,
                   studlab = studlab,
                   data = unwei_list[[i]],
                   sm = "SMD",
                   fixed = TRUE,
                   random = TRUE,
                   method.tau = "REML",
                   hakn = TRUE,
                   title = "EPIC and UKBB PRS Meta-analysis")
  
  meta$PRS_Outcome[i] <- m.gen$studlab[1]
  meta$fixed_OR[i] <- exp(m.gen$TE.fixed)
  meta$fixed_lower_CI[i] <- exp(m.gen$lower.fixed)
  meta$fixed_upper_CI[i] <- exp(m.gen$upper.fixed)
  meta$fixed_P_value[i] <- m.gen$pval.fixed
  meta$fixed_UKBB_weight[i] <- m.gen$w.fixed[1]
  meta$fixed_EPIC_weight[i] <- m.gen$w.fixed[2]
  meta$Q[i] <- m.gen$Q
  meta$Q_pval[i] <- m.gen$pval.Q
  meta$random_OR[i] <- exp(m.gen$TE.random)
  meta$random_lower_CI[i] <- exp(m.gen$lower.random)
  meta$random_upper_CI[i] <- exp(m.gen$upper.random)
  meta$random_P_value[i] <- m.gen$pval.random
  
}

unwei_meta <- meta

write.table(unwei_meta, "unweighted_PRS_meta_22-03-31.txt", row.names = F, sep = "|", quote = F)

