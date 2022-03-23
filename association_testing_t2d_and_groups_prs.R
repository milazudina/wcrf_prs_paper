setwd("~/Dropbox/PRS_clean_run/4_results/2022march_reanalysis/")
library(readxl)
library(data.table)
library(dplyr)

############## START OF FUNCTIONS #######################

# removing T2D cases from controls
summarise_model <- function(m, outfile = assoc_test_output, row_number) {
  
  # Get Model Summary and OR
  model_summary <- summary(m)
  #model_summary
  #exp(coef(model_summary))
  
  # Get Confidence Interval
  #conf_ints <- exp(confint(m))
  #conf_ints
  
  # Get P-values
  #p_vals <- coef(model_summary)[,4] # when p values are below <2e-16, you need to get the precise p-value using this
  #p_vals
  
  #assoc_test_output$OR[row_number] <- exp(coef(model_summary))[2,1]
  #assoc_test_output$lower_CI[row_number] <- conf_ints[2, 1]
  #assoc_test_output$upper_CI[row_number] <- conf_ints[2, 2]
  assoc_test_output$Beta[row_number] <- coef(model_summary)[2,1]
  assoc_test_output$StdErr[row_number] <- coef(model_summary)[2,2]
  assoc_test_output$P_value[row_number] <- coef(model_summary)[2,4]
  
  return(assoc_test_output)
  
}

merge_scores_and_pheno <- function(pheno_path, scores_path, pc_ethnicity_path = "~/Dropbox/UKBB_files/T2D_pheno/t2d_pheno_with_pcs_and_ethnicity_2022-03-23.txt"){
  
  pc_ethnicity = fread(pc_ethnicity_path) %>% select(IID, GeneticEthnicGrouping, AgeBaseline)
  pheno <- fread(pheno_path)
  scores <- read.table(scores_path, header=TRUE, sep="",stringsAsFactors = F)
  cat("N snps in the score:", scores$CNT[1], "\n")
  scores$SCORESUM_norm <- scale(scores$SCORESUM) # normalise the scoresum
  df <- merge(scores, pheno, by.x = "IID", by.y = "IID")
  
  if (pheno_path != pc_ethnicity_path){
    
    df = merge(df, pc_ethnicity, by.x = "IID", by.y = "IID")
    df <- df[df$IID > 0 & !is.na(df$GeneticEthnicGrouping), ]
    
  }
  
  return(df)
}

# don't remove T2D cases from controls in UKBB!
remove_t2d_cases_from_controls <- function(df, t2d_pheno_path, outcome) {
  t2d_pheno <- fread(t2d_pheno_path)
  dm <- merge(df, t2d_pheno[, c(1,15)], by = "IID", all.x = F, all.y = F)
  colnames(dm)[colnames(dm) == outcome] <- "temp"
  cat("Outcome:", outcome, "\n")
  cat("N cases:", sum(dm$temp == 1), "\n")
  cat("N controls:", sum(dm$temp == 0), "\n")
  dm2 <- dm[!(dm$temp == 0 & dm$T2D.all == 1), ]
  cat("N controls after the exclusion of T2D cases from controls:", sum(dm2$temp == 0), "\n")
  colnames(dm2)[colnames(dm2) == "temp"] <- outcome
  return(dm2)
}

############## END OF FUNCTIONS #######################

snpset_array <- c("group1_464", "group5_36", "group4_78", "group2_409", "group3_96", "t2d")
assoc_test_output_list <- list()

for (snpset in snpset_array){
  
  print(snpset)
  
  #### Pancreatic Cancer ####
  
  pheno_path <- "/Users/zudina_work/Dropbox/UKBB_files/Jareds_UKBB_cancer_pheno_files/panc_hosp_icd_ukbb.txt" # phenotype file
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
write.table(prs_assoc_table, "t2d_and_groups_prs_ukbb_2022-03-23.csv", sep = ",")
