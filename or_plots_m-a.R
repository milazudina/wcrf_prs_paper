library(dplyr)
library(tidyr)
library(ggrepel)
library(readxl)
library(forcats)

setwd("~/Dropbox/PRS_clean_run/PRS_meta-anlysis/2022march/")
ukbb_colour <- "#DC602E"
epic_colour <- "#05A8AA"
meta_colour <- "#7779a5"

#### T2D PRS ####

p = scores_for_plots %>%
  filter(PRS == "T2D", Outcome != "T2D") %>%
  select(Study, PRS, Outcome, OR, lower_CI, upper_CI, P_value, weight) %>%
  mutate(across(c(Study, Outcome, PRS), factor)) %>%
  mutate(Study = recode_factor(Study, fixed = "Meta-analysis", epic = "EPIC", ukbb = "UKBB")) %>% 
  mutate(Outcome = fct_relevel(Outcome, c("Pancreatic Cancer", "Colorectal Cancer", "Prostate Cancer", "Breast Cancer"))) %>%
  ggplot(., aes(x = Outcome, y = OR, ymin = lower_CI, ymax = upper_CI, color = Study, size = weight*100))+
  geom_point(position = position_dodge(0.7),
             shape = 16,
             cex = 1.5)+
  geom_point(position = position_dodge(0.7),
             shape = 15,
             alpha = 1)+
  geom_errorbar(lwd = 1.25,
                width = 0.3,
                position = position_dodge(0.7))+
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey")+
  geom_text(aes(label = ifelse(P_value <= 5e-2, formatC(P_value, format = "e", digits = 1), ""), group = Study), 
            hjust=-0.3, vjust=0.3, colour = 'black', size = 3, position = position_dodge(width = 0.7))+
  #geom_text_repel(aes(label = ifelse(P_value <= 5e-2, formatC(P_value, format = "e", digits = 2), ""), group = Study),
  #                size = 3, box.padding = 0.5, min.segment.length = 0.8, segment.alpha = 0, colour = "black", position = position_dodge(width = 0.7)) +
  scale_color_manual(values = c(ukbb_colour, epic_colour, meta_colour))+
  scale_size_continuous(name = "Weight, %",
                        breaks = c(10,50,90),
                        limits = c(0, 100),
                        labels = c(10,50,90),
                        range = c(0,6))+
  guides(size = "none")+
  theme_minimal()+
  theme(axis.title.x = element_text(vjust=-0.5, size = 12.5),
        axis.text.x.bottom = element_text(size = 12.5),
        axis.title.y = element_text(size = 12.5),
        axis.text.y.left = element_text(size = 12.5))

ggsave(p, file = "t2d_prs_cancer_outcome_meta_21-10-08_fin.png", height = 6, width = 9)

#### Cancer PRS ####

p = scores_for_plots %>%
  filter(!grepl("T2D|N", PRS)) %>%
  select(Study, PRS, Outcome, OR, lower_CI, upper_CI, P_value, weight) %>%
  mutate(across(c(Study, Outcome, PRS), factor)) %>%
  mutate(Study = recode_factor(Study, fixed = "Meta-analysis", epic = "EPIC", ukbb = "UKBB")) %>% 
  mutate(PRS = fct_relevel(PRS, c("Pancreatic cancer", "Colorectal cancer", "Prostate cancer", "Breast cancer"))) %>%
  ggplot(., aes(x = PRS, y = OR, ymin = lower_CI, ymax = upper_CI, color = Study, size = weight*100))+
  geom_point(position = position_dodge(0.7),
             shape = 16,
             cex = 1.5)+
  geom_point(position = position_dodge(0.7),
             shape = 15,
             alpha = 1)+
  geom_errorbar(lwd = 1.25,
                width = 0.3,
                position = position_dodge(0.7))+
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey")+
  geom_text(aes(label = ifelse(P_value <= 5e-2, formatC(P_value, format = "e", digits = 1), ""), group = Study), 
            hjust=-0.3, vjust=0.3, colour = 'black', size = 3, position = position_dodge(width = 0.7))+
  #geom_text_repel(aes(label = ifelse(P_value <= 5e-2, formatC(P_value, format = "e", digits = 2), ""), group = Study),
  #                size = 3, box.padding = 0.5, min.segment.length = 0.8, segment.alpha = 0, colour = "black", position = position_dodge(width = 0.7)) +
  scale_color_manual(values = c(ukbb_colour, epic_colour, meta_colour))+
  scale_size_continuous(name = "Weight, %",
                        breaks = c(10,50,90),
                        limits = c(0, 100),
                        labels = c(10,50,90),
                        range = c(0,6))+
  guides(size = "none")+
  theme_minimal()+
  theme(axis.title.x = element_text(vjust=-0.5, size = 12.5),
        axis.text.x.bottom = element_text(size = 12.5),
        axis.title.y = element_text(size = 12.5),
        axis.text.y.left = element_text(size = 12.5))

ggsave(p, file = "cancer_prs_t2d_outcome_meta_21-10-08_fin.png", height = 6, width = 9)

#### Plots with all groups together - META_ANALYSIS ####

plot_groups_betas = scores_for_plots %>%
  filter(!grepl("cancer|T2D", PRS), Study == "fixed") %>%
  select(Study, PRS, Outcome, OR, lower_CI, upper_CI, P_value, weight) %>%
  mutate(across(c(Study, Outcome, PRS), factor)) %>%
  mutate(Study = recode_factor(Study, fixed = "Meta-analysis"),
         PRS = recode_factor(PRS, "1b. Higher Adiposity, Insulin Resistance, Lower Age at Menarche and decreased Sex Hormones (N = 36)" = "1b. Higher Adiposity, Insulin Resistance,\nLower Age at Menarche and decreased Sex Hormones (N = 36)",
                             "3. Decreased lipids and Blood Pressure (N=409)" = "3. Higher Blood Pressure and lower TG, TC, LDL-C (N=409)")) %>% 
  mutate(Outcome = fct_relevel(Outcome, c("T2D", "Pancreatic Cancer", "Colorectal Cancer", "Prostate Cancer", "Breast Cancer")),
         PRS = fct_relevel(PRS, sort)) %>%
  ggplot(., aes(x = Outcome, y = OR, ymin = lower_CI, ymax = upper_CI, color = PRS))+
  geom_point(position = position_dodge(0.7),
             cex = 2)+
  geom_errorbar(lwd = 1.25,
                width = 0.3,
                position = position_dodge(0.7))+
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey")+
  #geom_text(aes(label = ifelse(P_value <= 5e-2, formatC(P_value, format = "e", digits = 1), ""), group = PRS), 
  #          hjust=-0.3, vjust=0.3, colour = 'black', size = 3, position = position_dodge(width = 0.7))+
  geom_text_repel(aes(label = ifelse(P_value <= 5e-2, formatC(P_value, format = "e", digits = 2), ""), group = PRS),
                  size = 3, box.padding = 0.5, min.segment.length = 0.8, segment.alpha = 0, colour = "black", position = position_dodge(width = 0.7)) +
  scale_colour_manual(name = "Group", values = c("#BF1A2F", "#454E9E", "#F00699", "#018E42", "#F7D002"), guide = guide_legend())+
  theme_minimal()+
  ggtitle("Fixed Effects Meta-analysis")+
  theme(axis.title.x = element_text(vjust=-0.5, size = 9),
        axis.text.x.bottom = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.y.left = element_text(size = 9))

ggsave(plot_groups_betas, file = "all_groups_weighted_meta_21-10-11.png", height = 6, width = 14)

#### Plots with all groups together - EPIC ####

scores_for_plots_groups_betas <- scores_for_plots %>%
  filter(PRS != "T2D", Weighting == "Betas", Cohort == "EPIC", Outcome != "T2D") %>%
  select(Cohort, PRS, Outcome, OR, lower_CI, upper_CI, P_value) %>%
  mutate(Outcome = fct_relevel(Outcome, c("T2D", "Pancreatic Cancer", "Colorectal Cancer", "Prostate Cancer", "Breast Cancer")))

#plot_groups_betas <- 
ggplot(scores_for_plots_groups_betas, aes(x = Outcome, y = OR, ymin = lower_CI, ymax = upper_CI, color = PRS))+
  geom_point(position = position_dodge(0.7),
             cex = 2)+
  geom_errorbar(lwd = 1.25,
                width = 0.3,
                position = position_dodge(0.7))+
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey")+
  geom_text(aes(label = ifelse(P_value <= 5e-2, formatC(P_value, format = "e", digits = 1), ""), group = PRS), 
            hjust=-0.3, vjust=0.3, colour = 'black', size = 3, position = position_dodge(width = 0.7))+
  #geom_text_repel(aes(label = ifelse(P_value <= 5e-2, formatC(P_value, format = "e", digits = 2), ""), group = Weighting),
  #                size = 3, box.padding = 0.5, min.segment.length = 0.8, segment.alpha = 0, colour = "black", position = position_dodge(width = 0.7)) +
  scale_colour_manual(name = "Group", values = c("#BF1A2F", "#454E9E", "#F00699", "#018E42", "#F7D002"), guide = guide_legend())+
  theme_minimal()+
  ggtitle("EPIC")+
  theme(axis.title.x = element_text(vjust=-0.5, size = 12.5),
        axis.text.x.bottom = element_text(size = 12.5),
        axis.title.y = element_text(size = 12.5),
        axis.text.y.left = element_text(size = 12.5))

ggsave(plot_groups_betas, file = "all_groups_weighted_meta_21-09-09.png", height = 6, width = 12)
