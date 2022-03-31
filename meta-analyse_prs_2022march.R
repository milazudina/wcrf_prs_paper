# meta analysis
library(meta)
library(readxl)
library(tidyr)
library(dplyr)
library(data.table)

`%nin%` = Negate(`%in%`)

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

# this is for t2d & groups
# here we need to exclude all with the Weighting == UKBB and then rename UKBB without ambiguous as UKBB
df2 <- df %>% 
  filter(Adjustment == "6 PCs", Cohort != "UKBB") %>% 
  select(SNPs, Cohort, Weighting, Outcome, OR, lower_CI, upper_CI, P_value) %>% 
  mutate(Beta = log(OR), 
         SE = (upper_CI - lower_CI)/3.92,
         CI = upper_CI - lower_CI) %>%
  mutate(Cohort=recode_factor(Cohort, `UKBB without ambiguous SNPs` = "UKBB"))

df3 <- read_excel("~/Dropbox/WCRF_Manuscript1/PRS_meta-anlysis/ukbb_and_epic_prs+new_21-10-06.xlsx", sheet = 2)

# this is for cancer prss
df4 <- df3 %>% 
  filter(Adjustment == "6 PCs") %>%
  select(SNPs, Cohort, Weighting, Outcome, OR, lower_CI, upper_CI, P_value) %>% 
  mutate(Beta = log(OR), 
         SE = (upper_CI - lower_CI)/3.92,
         CI = upper_CI - lower_CI)

df5 <- rbind(df2, df4)

df6 <- df5 %>% filter(!grepl("M-A", Cohort)) %>% 
  mutate(Cohort = factor(Cohort, levels = c("UKBB", "EPIC"))) %>%
  mutate(Outcome = factor(Outcome, levels = c("T2D", "Pancreatic Cancer", "Colorectal Cancer", "Prostate Cancer", "Breast Cancer")))
 
wei <- df6 %>% filter(Weighting == "Betas") 
unwei <- df6 %>% filter(Weighting == "Unweighted")

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
                   comb.fixed = FALSE,
                   comb.random = TRUE,
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

write.table(wei_meta, "~/Dropbox/PRS_meta_21-10-07.txt", row.names = F, sep = "|", quote = F)


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
                   comb.fixed = FALSE,
                   comb.random = TRUE,
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