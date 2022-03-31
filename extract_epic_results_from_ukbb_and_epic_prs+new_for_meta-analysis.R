library(readxl)
library(dplyr)
`%nin%` = Negate(`%in%`)
setwd("~/Dropbox/PRS_clean_run/4_results/2022march_reanalysis/")

df1 = read_excel("../ukbb_and_epic_prs+new_21-10-06.xlsx", sheet = 1)
df2 = read_excel("../ukbb_and_epic_prs+new_21-10-06.xlsx", sheet = 2)
df2$N_cases = NA
df2$N_controls = NA

str(df1)
df1 %>% select(where(is.character)) %>% apply(., 2, function(x){levels(as.factor(x))})

epic_t2d_and_group_prs = df1 %>% 
  mutate(across(is.character, factor)) %>% 
  filter(Cohort == "EPIC") %>% 
  select(-c(colnames(df1)[colnames(df1) %nin% colnames(df2)], "Aligned to"))


df2 %>% select(where(is.character)) %>% apply(., 2, function(x){levels(as.factor(x))})

epic_cancer = df2 %>% 
  mutate(across(is.character, factor)) %>% 
  filter(Cohort == "EPIC") %>%
  select(-c(colnames(df2)[colnames(df2) %nin% colnames(df1)], "Aligned to"))

epic_prss = rbind(epic_t2d_and_group_prs, epic_cancer)

write.table(epic_prss, "../../PRS_meta-anlysis/2022march/epic_prs_for_m-a.csv", row.names = F, quote = T, sep = ",")
