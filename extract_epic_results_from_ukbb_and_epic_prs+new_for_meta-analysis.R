library(readxl)
library(dplyr)
`%nin%` = Negate(`%in%`)
setwd("~/Dropbox/PRS_clean_run/4_results/2022march_reanalysis/")

df1 = read_excel("../ukbb_and_epic_prs+new_21-10-06.xlsx", sheet = 1)
df2 = read_excel("../ukbb_and_epic_prs+new_21-10-06.xlsx", sheet = 2)

str(df1)
df1 %>% select(where(is.character)) %>% apply(., 2, function(x){levels(as.factor(x))})

epic_t2d_and_group_prs = df1 %>% 
  mutate(across(is.character, factor)) %>% 
  filter(Cohort == "EPIC") %>% 
  select(-c(colnames(epic_t2d_and_group_prs)[colnames(epic_t2d_and_group_prs) %nin% colnames(epic_cancer)], "Aligned to"))


df2 %>% select(where(is.character)) %>% apply(., 2, function(x){levels(as.factor(x))})

epic_cancer = df2 %>% 
  mutate(across(is.character, factor)) %>% 
  filter(Cohort == "EPIC") %>%
  select(-c(colnames(epic_cancer)[colnames(epic_cancer) %nin% colnames(epic_t2d_and_group_prs)], "Aligned to"))

epic_prss = rbind(epic_t2d_and_group_prs, epic_cancer)

write.table(epic_prss, "../../PRS_meta-anlysis/2022march/epic_prs_for_m-a.csv", row.names = F, quote = T, sep = ",")
