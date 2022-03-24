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
  
  # the commented lines - this is a slow way to calculate OR & confints - instead, I do it for all of them together
  
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
    
  }
  
  df = df[df$IID > 0 & !is.na(df$GeneticEthnicGrouping), ]
  
  return(df)
}

# don't remove T2D cases from controls in UKBB! this was just for the sensitivity analysis
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
