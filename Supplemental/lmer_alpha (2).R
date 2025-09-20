library(lme4)
library(lmerTest)
library(tidyverse)
library(MuMIn)

data <- read.csv('log_ratio_foods_analysis.csv', 
                 row.names=1)
cols_to_exclude <- c('thdmi_cohort', 'covid_level_of_wellbeing_coded', 
                     'bmi_cat_coded', 'antibiotic_history_coded', 
                     'faecali_clr', 
                     'prevotella_clr', 'Amount_Energy_in_kcal', 
                     'sex', 'Amount_Energy_in_kj')
variables_to_test <- colnames(data)[!(colnames(data) %in% cols_to_exclude)]
#variables_to_test <- c('Processed_Calories_Nova_ultra_processed_foods_per1000kcal')
log_ratios <- c('faecali_clr', 
                'prevotella_clr')


out_df <- data.frame()
for(x in 1:length(variables_to_test)){
  # remove values that are NA 
  no_nans = data %>% drop_na(variables_to_test[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(log_ratios)){
    form = paste(c(log_ratios[y], '~', variables_to_test[x], '+', 
                   'covid_level_of_wellbeing_coded', '+', 'bmi_cat_coded', '+', 
                   'antibiotic_history_coded', '+', '(1 | thdmi_cohort)'), 
                 collapse = " ")
    ex <- lmer(formula(form), data=data)
    out <- summary(ex)
    r2 <- r.squaredGLMM(ex)
      
    # Define the coefficient name to check, assuming we're interested in the effect of the tested variable
    coefficient_name <- variables_to_test[x]
    
    # Check if the coefficient exists and handle negative values of the estimate
    if (coefficient_name %in% rownames(out$coefficients) && !is.na(out$coefficients[coefficient_name, "Estimate"])) {
      if (out$coefficients[coefficient_name, "Estimate"] < 0) {
        marginal_r2 <- -1 * r2[1]
        conditional_r2 <- -1 * r2[2]
      } else {
        marginal_r2 <- r2[1]
        conditional_r2 <- r2[2]
      }
      p_value <- out$coefficients[coefficient_name, "Pr(>|t|)"]
    } else {
      marginal_r2 <- NA
      conditional_r2 <- NA
      p_value <- NA
      warning(paste("Coefficient", coefficient_name, "is missing or NA in model for variable", log_ratios[y]))
    }

    out_df <- rbind(out_df, c(variables_to_test[x], log_ratios[y], p_value, marginal_r2, conditional_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'marginal_r2', 'conditional_r2')
write.csv(out_df, file='log_ratio_analysis.csv')

