setwd("~/Dropbox/prs_24_08_21/")
library(readxl)
library(data.table)

# removing T2D cases from controls

summarise_model <- function(m, outfile = assoc_test_output, row_number) {
  
  # Get Model Summary and OR
  model_summary <- summary(m)
  model_summary
  exp(coef(model_summary))
  
  # Get Confidence Interval
  conf_ints <- exp(confint(m))
  conf_ints
  
  # Get P-values
  p_vals <- coef(model_summary)[,4] # when p values are below <2e-16, you need to get the precise p-value using this
  p_vals
  
  assoc_test_output$OR[row_number] <- exp(coef(model_summary))[2,1]
  assoc_test_output$lower_CI[row_number] <- conf_ints[2, 1]
  assoc_test_output$upper_CI[row_number] <- conf_ints[2, 2]
  assoc_test_output$P_value[row_number] <- coef(model_summary)[2,4]
  
  return(assoc_test_output)
  
}
merge_scores_and_pheno <- function(pheno_path, scores_path){
  
  pheno <- fread(pheno_path)
  scores <- read.table(scores_path, header=TRUE, sep="",stringsAsFactors = F)
  cat("N snps in the score:", scores$CNT[1], "\n")
  scores$SCORESUM_norm <- scale(scores$SCORESUM) # normalise the scoresum
  df <- merge(scores, pheno, by.x = "IID", by.y = "IID")
  df <- df[df$IID > 0, ]
  
  return(df)
}

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
