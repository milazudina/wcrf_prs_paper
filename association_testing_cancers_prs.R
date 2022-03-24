setwd("~/Dropbox/PRS_clean_run/4_results/2022march_reanalysis")
library(readxl)
library(data.table)
source("~/Dropbox/PRS_clean_run/scripts/association_testing_functions.R")

pheno_path <- "/Users/zudina_work/Dropbox/prs_ukbb_replication/phenotype_files/t2d_ukbb_with_PCs.txt"
snpset_array <- c("brc", "crc", "panc", "prc")
assoc_test_output_list <- list()

for (snpset in snpset_array){
  
  print(snpset)
  
  scores_path <- paste("step6_out_24_08_21/ukbb_", snpset, "_snps_", snpset, "_betas_score.profile", sep = "") # scores file
  df <- merge_scores_and_pheno(pheno_path, scores_path)

  assoc_test_output <- data.frame(OR = c(rep(NA, 6)), lower_CI = c(rep(NA, 6)), upper_CI = c(rep(NA, 6)), P_value = c(rep(NA, 6)),
                                  row.names = c("Betas unadjusted", 
                                                "Betas 6PCs", 
                                                "Z-scores unadjusted",
                                                "Z-scores 6PCs",
                                                "Unweighted unadjusted", 
                                                "Unweighted 6PCs"))


  m1 <- glm(T2D.all ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m1, assoc_test_output, 1)
  print(assoc_test_output[1, ])

  m2 <- glm(T2D.all ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m2, assoc_test_output, 2)
  print(assoc_test_output[2, ])

  ## Weighted by Z-scores
  scores_path <- paste("step6_out_24_08_21/ukbb_", snpset, "_snps_", snpset, "_zscores_score.profile", sep = "") # scores file # scores file
  df <- merge_scores_and_pheno(pheno_path, scores_path)

  m3 <- glm(T2D.all ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m3, assoc_test_output, 3)
  print(assoc_test_output[3, ])

  m4 <- glm(T2D.all ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m4, assoc_test_output, 4)
  print(assoc_test_output[4, ])

  ## Unweighted
  scores_path <- paste("step6_out_24_08_21/ukbb_", snpset, "_snps_unweighted_score.profile", sep = "") # scores file # scores file
  df <- merge_scores_and_pheno(pheno_path, scores_path)

  m5 <- glm(T2D.all ~ SCORESUM_norm, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m5, assoc_test_output, 5)
  print(assoc_test_output[5, ])

  m6 <- glm(T2D.all ~ SCORESUM_norm+PC1+PC2+PC3+PC4+PC5+PC6, family="binomial"(link="logit"), data=df)
  assoc_test_output <- summarise_model(m6, assoc_test_output, 6)
  print(assoc_test_output[6, ])

  assoc_test_output_list[[snpset]] <- assoc_test_output
  View(assoc_test_output_list[[snpset]])
    
}
