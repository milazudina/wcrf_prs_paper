setwd("~/Dropbox/PRS_clean_run/4_results/2022march_reanalysis/")
library(readxl)
library(data.table)
library(dplyr)
source("~/Dropbox/PRS_clean_run/scripts/association_testing_functions.R")

snpset_array <- c("group1_464", "group5_36", "group4_78", "group2_409", "group3_96", "t2d")
assoc_test_output_list <- list()

for (snpset in snpset_array){
  
  print(snpset)
  
  #### Pancreatic Cancer ####
  
  pheno_path <- "~/Dropbox/UKBB_files/Jareds_UKBB_cancer_pheno_files/panc_hosp_icd_ukbb.txt" # phenotype file
  # Weighted by Betas
  scores_path <- paste("~/Dropbox/PRS_clean_run/3_scores/ukbb_", snpset, "_weighted_score.profile", sep = "") # scores file
  outcome <- "panc"
  
  df <- merge_scores_and_pheno(pheno_path, scores_path)
  
  cat("Outcome:", outcome, "\n")
  cat("N cases:", sum(df[, outcome] == 1), "\n")
  cat("N controls:", sum(df[, outcome] == 0), "\n")
  
  assoc_test_output <- data.frame(SNPs = c(rep(snpset, 6)),
                                  N_cases = c(rep(sum(df[, outcome] == 1), 6)),
                                  N_controls = c(rep(sum(df[, outcome] == 0), 6)),
                                  Weighting = c("Weighted", "Weighted","Unweighted", "Unweighted", "Weighted", "Unweighted"),
                                  Outcome = c(rep(outcome, 6)),
                                  Adjustment = c("Unadjusted", "6 PCs", "Unadjusted", "6 PCs", "6 PCs + Age", "6 PCs + Age"),
                                  OR = c(rep(NA, 6)), 
                                  lower_CI = c(rep(NA, 6)), 
                                  upper_CI = c(rep(NA, 6)),
                                  Beta = c(rep(NA,6)),
                                  StdErr = c(rep(NA,6)),
                                  P_value = c(rep(NA, 6)))
  
  m1 <- glm(panc ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m1, assoc_test_output, 1)
  print(assoc_test_output[1, ])
  
  m2 <- glm(panc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m2, assoc_test_output, 2)
  print(assoc_test_output[2, ])
  
  ## Weighted with age
  m5 <- glm(panc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6+AgeBaseline, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m5, assoc_test_output, 5)
  print(assoc_test_output[5, ])
  
  ## Unweighted
  scores_path <- paste("~/Dropbox/PRS_clean_run/3_scores/ukbb_", snpset, "_unweighted_score.profile", sep = "") # scores file # scores file
  df <- merge_scores_and_pheno(pheno_path, scores_path)
  
  m3 <- glm(panc ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m3, assoc_test_output, 3)
  print(assoc_test_output[3, ])
  
  m4 <- glm(panc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m4, assoc_test_output, 4)
  print(assoc_test_output[4, ])
  
  ## Unweighted with age
  m6 <- glm(panc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6+AgeBaseline, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m6, assoc_test_output, 6)
  print(assoc_test_output[6, ])
  
  assoc_test_output$OR = exp(assoc_test_output$Beta)
  assoc_test_output$lower_CI = exp(assoc_test_output$Beta - 1.96*assoc_test_output$StdErr)
  assoc_test_output$upper_CI = exp(assoc_test_output$Beta + 1.96*assoc_test_output$StdErr)
  
  assoc_test_output_panc <- assoc_test_output
  
  #### Colorectal Cancer ####
  
  pheno_path <- "~/Dropbox/UKBB_files/Jareds_UKBB_cancer_pheno_files/crc_hosp_icd_ukbb.txt" # phenotype file
  # Weighted by Betas
  scores_path <- paste("~/Dropbox/PRS_clean_run/3_scores/ukbb_", snpset, "_weighted_score.profile", sep = "") # scores file
  outcome <- "crc"
  
  df <- merge_scores_and_pheno(pheno_path, scores_path)

  cat("Outcome:", outcome, "\n")
  cat("N cases:", sum(df[, outcome] == 1), "\n")
  cat("N controls:", sum(df[, outcome] == 0), "\n")
  
  assoc_test_output <- data.frame(SNPs = c(rep(snpset, 6)),
                                  N_cases = c(rep(sum(df[, outcome] == 1), 6)),
                                  N_controls = c(rep(sum(df[, outcome] == 0), 6)),
                                  Weighting = c("Weighted", "Weighted","Unweighted", "Unweighted", "Weighted", "Unweighted"),
                                  Outcome = c(rep(outcome, 6)),
                                  Adjustment = c("Unadjusted", "6 PCs", "Unadjusted", "6 PCs", "6 PCs + Age", "6 PCs + Age"),
                                  OR = c(rep(NA, 6)), 
                                  lower_CI = c(rep(NA, 6)), 
                                  upper_CI = c(rep(NA, 6)),
                                  Beta = c(rep(NA,6)),
                                  StdErr = c(rep(NA,6)),
                                  P_value = c(rep(NA, 6)))
  
  m1 <- glm(crc ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m1, assoc_test_output, 1)
  print(assoc_test_output[1, ])
  
  m2 <- glm(crc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m2, assoc_test_output, 2)
  print(assoc_test_output[2, ])
  
  ## Weighted with age
  m5 <- glm(crc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6+AgeBaseline, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m5, assoc_test_output, 5)
  print(assoc_test_output[5, ])
  
  ## Unweighted
  scores_path <- paste("~/Dropbox/PRS_clean_run/3_scores/ukbb_", snpset, "_unweighted_score.profile", sep = "") # scores file # scores file
  df <- merge_scores_and_pheno(pheno_path, scores_path)
  
  m3 <- glm(crc ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m3, assoc_test_output, 3)
  print(assoc_test_output[3, ])
  
  m4 <- glm(crc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m4, assoc_test_output, 4)
  print(assoc_test_output[4, ])
  
  ## Unweighted with age
  m6 <- glm(crc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6+AgeBaseline, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m6, assoc_test_output, 6)
  print(assoc_test_output[6, ])
  
  assoc_test_output$OR = exp(assoc_test_output$Beta)
  assoc_test_output$lower_CI = exp(assoc_test_output$Beta - 1.96*assoc_test_output$StdErr)
  assoc_test_output$upper_CI = exp(assoc_test_output$Beta + 1.96*assoc_test_output$StdErr)
  
  assoc_test_output_crc <- assoc_test_output
  
  #### Prostate Cancer ####
  pheno_path <- "~/Dropbox/UKBB_files/Jareds_UKBB_cancer_pheno_files/prc_hosp_icd_ukbb.txt" # phenotype file
  # Weighted by Betas
  scores_path <- paste("~/Dropbox/PRS_clean_run/3_scores/ukbb_", snpset, "_weighted_score.profile", sep = "")
  outcome <- "prc"
  
  df <- merge_scores_and_pheno(pheno_path, scores_path)
  #df <- remove_t2d_cases_from_controls(df, t2d_pheno_path, outcome)
  cat("Outcome:", outcome, "\n")
  cat("N cases:", sum(df[, outcome] == 1), "\n")
  cat("N controls:", sum(df[, outcome] == 0), "\n")
  
  assoc_test_output <- data.frame(SNPs = c(rep(snpset, 6)),
                                  N_cases = c(rep(sum(df[, outcome] == 1), 6)),
                                  N_controls = c(rep(sum(df[, outcome] == 0), 6)),
                                  Weighting = c("Weighted", "Weighted","Unweighted", "Unweighted", "Weighted", "Unweighted"),
                                  Outcome = c(rep(outcome, 6)),
                                  Adjustment = c("Unadjusted", "6 PCs", "Unadjusted", "6 PCs", "6 PCs + Age", "6 PCs + Age"),
                                  OR = c(rep(NA, 6)), 
                                  lower_CI = c(rep(NA, 6)), 
                                  upper_CI = c(rep(NA, 6)),
                                  Beta = c(rep(NA,6)),
                                  StdErr = c(rep(NA,6)),
                                  P_value = c(rep(NA, 6)))
  
  m1 <- glm(prc ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m1, assoc_test_output, 1)
  print(assoc_test_output[1, ])
  
  m2 <- glm(prc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m2, assoc_test_output, 2)
  print(assoc_test_output[2, ])
  
  ## Weighted with age
  m5 <- glm(prc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6+AgeBaseline, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m5, assoc_test_output, 5)
  print(assoc_test_output[5, ])
  
  ## Unweighted
  scores_path <- paste("~/Dropbox/PRS_clean_run/3_scores/ukbb_", snpset, "_unweighted_score.profile", sep = "") # scores file # scores file
  df <- merge_scores_and_pheno(pheno_path, scores_path)
  
  m3 <- glm(prc ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m3, assoc_test_output, 3)
  print(assoc_test_output[3, ])
  
  m4 <- glm(prc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m4, assoc_test_output, 4)
  print(assoc_test_output[4, ])
  
  ## Unweighted with age
  m6 <- glm(prc ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6+AgeBaseline, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m6, assoc_test_output, 6)
  print(assoc_test_output[6, ])
  
  assoc_test_output$OR = exp(assoc_test_output$Beta)
  assoc_test_output$lower_CI = exp(assoc_test_output$Beta - 1.96*assoc_test_output$StdErr)
  assoc_test_output$upper_CI = exp(assoc_test_output$Beta + 1.96*assoc_test_output$StdErr)
  
  assoc_test_output_prc <- assoc_test_output
  
  #### Breast Cancer ####
  
  pheno_path <- "~/Dropbox/UKBB_files/Jareds_UKBB_cancer_pheno_files/brc_hosp_icd_ukbb.txt" # phenotype file
  # Weighted by Betas
  scores_path <- paste("~/Dropbox/PRS_clean_run/3_scores/ukbb_", snpset, "_weighted_score.profile", sep = "") # scores file
  outcome <- "bc2"
  
  df <- merge_scores_and_pheno(pheno_path, scores_path)
  #df <- remove_t2d_cases_from_controls(df, t2d_pheno_path, outcome)
  cat("Outcome:", outcome, "\n")
  cat("N cases:", sum(df[, outcome] == 1), "\n")
  cat("N controls:", sum(df[, outcome] == 0), "\n")
  
  assoc_test_output <- data.frame(SNPs = c(rep(snpset, 6)),
                                  N_cases = c(rep(sum(df[, outcome] == 1), 6)),
                                  N_controls = c(rep(sum(df[, outcome] == 0), 6)),
                                  Weighting = c("Weighted", "Weighted","Unweighted", "Unweighted", "Weighted", "Unweighted"),
                                  Outcome = c(rep(outcome, 6)),
                                  Adjustment = c("Unadjusted", "6 PCs", "Unadjusted", "6 PCs", "6 PCs + Age", "6 PCs + Age"),
                                  OR = c(rep(NA, 6)), 
                                  lower_CI = c(rep(NA, 6)), 
                                  upper_CI = c(rep(NA, 6)),
                                  Beta = c(rep(NA,6)),
                                  StdErr = c(rep(NA,6)),
                                  P_value = c(rep(NA, 6)))
  
  m1 <- glm(bc2 ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m1, assoc_test_output, 1)
  print(assoc_test_output[1, ])
  
  m2 <- glm(bc2 ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m2, assoc_test_output, 2)
  print(assoc_test_output[2, ])
  
  ## Weighted with age
  m5 <- glm(bc2 ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6+AgeBaseline, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m5, assoc_test_output, 5)
  print(assoc_test_output[5, ])
  
  ## Unweighted
  scores_path <- paste("~/Dropbox/PRS_clean_run/3_scores/ukbb_", snpset, "_unweighted_score.profile", sep = "") # scores file # scores file
  df <- merge_scores_and_pheno(pheno_path, scores_path)
  
  m3 <- glm(bc2 ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m3, assoc_test_output, 3)
  print(assoc_test_output[3, ])
  
  m4 <- glm(bc2 ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m4, assoc_test_output, 4)
  print(assoc_test_output[4, ])
  
  ## Unweighted with age
  m6 <- glm(bc2 ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6+AgeBaseline, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m6, assoc_test_output, 6)
  print(assoc_test_output[6, ])
  
  assoc_test_output_brc <- assoc_test_output
  
  #### T2D ####

  pheno_path <- "~/Dropbox/UKBB_files/T2D_pheno/t2d_pheno_with_pcs_and_ethnicity_2022-03-23.txt" # phenotype file
  # Weighted by Betas
  scores_path <- paste("~/Dropbox/PRS_clean_run/3_scores/ukbb_", snpset, "_weighted_score.profile", sep = "") # scores file
  outcome <- "T2D.all"
  
  df <- merge_scores_and_pheno(pheno_path, scores_path)
  cat("Outcome:", outcome, "\n")
  cat("N cases:", sum(df[, outcome] == 1), "\n")
  cat("N controls:", sum(df[, outcome] == 0), "\n")

  assoc_test_output <- data.frame(SNPs = c(rep(snpset, 6)),
                                  N_cases = c(rep(sum(df[, outcome] == 1), 6)),
                                  N_controls = c(rep(sum(df[, outcome] == 0), 6)),
                                  Weighting = c("Weighted", "Weighted","Unweighted", "Unweighted", "Weighted", "Unweighted"),
                                  Outcome = c(rep(outcome, 6)),
                                  Adjustment = c("Unadjusted", "6 PCs", "Unadjusted", "6 PCs", "6 PCs + Age", "6 PCs + Age"),
                                  OR = c(rep(NA, 6)), 
                                  lower_CI = c(rep(NA, 6)), 
                                  upper_CI = c(rep(NA, 6)),
                                  Beta = c(rep(NA,6)),
                                  StdErr = c(rep(NA,6)),
                                  P_value = c(rep(NA, 6)))

  m1 <- glm(T2D.all ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m1, assoc_test_output, 1)
  print(assoc_test_output[1, ])

  m2 <- glm(T2D.all ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m2, assoc_test_output, 2)
  print(assoc_test_output[2, ])
  
  m5 <- glm(T2D.all ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6+AgeBaseline, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m5, assoc_test_output, 5)
  print(assoc_test_output[5, ])

  ## Unweighted
  scores_path <- paste("~/Dropbox/PRS_clean_run/3_scores/ukbb_", snpset, "_unweighted_score.profile", sep = "") # scores file # scores file
  df <- merge_scores_and_pheno(pheno_path, scores_path)

  m3 <- glm(T2D.all ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m3, assoc_test_output, 3)
  print(assoc_test_output[3, ])

  m4 <- glm(T2D.all ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m4, assoc_test_output, 4)
  print(assoc_test_output[4, ])
  
  ## Unweighted with age
  m6 <- glm(T2D.all ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6+AgeBaseline, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m6, assoc_test_output, 6)
  print(assoc_test_output[6, ])
  
  assoc_test_output$OR = exp(assoc_test_output$Beta)
  assoc_test_output$lower_CI = exp(assoc_test_output$Beta - 1.96*assoc_test_output$StdErr)
  assoc_test_output$upper_CI = exp(assoc_test_output$Beta + 1.96*assoc_test_output$StdErr)

  assoc_test_output_t2d <- assoc_test_output
  
  #### Concatenate them ####
  assoc_test_output <- rbind(assoc_test_output_panc, assoc_test_output_crc, assoc_test_output_prc, assoc_test_output_brc, assoc_test_output_t2d)

  assoc_test_output_list[[snpset]] <- assoc_test_output
  View(assoc_test_output_list[[snpset]])
  
}

prs_assoc_table <- rbindlist(assoc_test_output_list)
write.table(prs_assoc_table, "t2d_and_groups_prs_ukbb_2022-03-23.csv", sep = ",", row.names = F)
