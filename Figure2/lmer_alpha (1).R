library(lme4)
library(lmerTest)
library(tidyverse)
library(MuMIn)

data <- read.csv('log_ratio_foods_analysis.csv', 
                 row.names=1)
cols_to_exclude <- c('thdmi_cohort', 'covid_level_of_wellbeing_coded', 
                     'bmi_cat_coded', 'antibiotic_history_coded', 
                     'log_Faecalibacterium_to_Bacteroides', 
                     'log_Prevotella_to_Bacteroides', 'Amount_Energy_in_kcal', 
                     'sex', 'Amount_Energy_in_kj')
variables_to_test <- colnames(data)[!(colnames(data) %in% cols_to_exclude)]
#variables_to_test <- c('Processed_Calories_Nova_ultra_processed_foods_per1000kcal')
log_ratios <- c('log_Faecalibacterium_to_Bacteroides', 
                'log_Prevotella_to_Bacteroides')

prev_df <- read.csv('prevotella_cols.csv', 
                    row.names=2)

# prev data
prev_lmer <- read.csv('prevotella_lmer.csv', 
                      row.names=1)
prev_high <- read.csv('prevotella_high.csv', 
                      row.names=1)
prev_low <- read.csv('prevotella_low.csv', 
                     row.names=1)
prev_mx <- read.csv('prevotella_mexico.csv', 
                    row.names=1)
prev_uk <- read.csv('prevotella_uk.csv', 
                    row.names=1)
prev_us <- read.csv('prevotella_us.csv', 
                    row.names=1)

faec_lmer <- read.csv('faecalibacterium_lmer.csv', 
                      row.names=1)
faec_cols <- colnames(faec_lmer)[49:67]
faec_sum <- rowSums(faec_lmer[ ,faec_cols])

faec_mx <- faec_lmer[faec_lmer['thdmi_cohort'] == 'Mexico', ]
faec_uk <- faec_lmer[faec_lmer['thdmi_cohort'] == 'UK', ]
faec_us <- faec_lmer[faec_lmer['thdmi_cohort'] == 'US', ]

prevotella_strains <- rownames(prev_df)
variables_prev <- colnames(data_prev)[!(colnames(data_prev) %in% prevotella_strains)]
variables_prev <- variables_prev[!(variables_prev) %in% cols_to_exclude]

variables_to_test <- c('DP_Core', 'DP_All')
out_df <- data.frame()
for(x in 1:length(variables_to_test)){
  # remove values that are NA 
  no_nans = data %>% drop_na(variables_to_test[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(log_ratios)){
    form = paste(c(log_ratios[y], '~', variables_to_test[x], '+', 
                   'covid_level_of_wellbeing_coded', '+', 
                   'Amount_Energy_in_kcal', '+', 'bmi_cat_coded', '+', 
                   'antibiotic_history_coded', '+', '(1 | thdmi_cohort)'), 
                 collapse = " ")
    ex <- lmer(formula(form), data=data)
    out <- summary(ex)
    r2 <- r.squaredGLMM(ex)
    if (out$coefficients[23] < 0) {
      marginal_r2 <- -1* r2[1]
      conditional_r2 <- -1* r2[2]  
    } else {
      marginal_r2 <- r2[1]
      conditional_r2 <- r2[2]
    }
    p_value <- out$coefficients[30]
    out_df <- rbind(out_df, c(variables_to_test[x], log_ratios[y], p_value, marginal_r2, conditional_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'marginal_r2', 'conditional_r2')
write.csv(out_df, file='log_ratio_analysis_calories_noenergy.csv')

# all prevotella
out_df <- data.frame()
for(x in 1:length(variables_prev)){
  # remove values that are NA 
  no_nans = prev_lmer %>% drop_na(variables_prev[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(prevotella_strains)){
    form = paste(c(prevotella_strains[y], '~', variables_prev[x], '+', 
                   'sex', '+', 'covid_level_of_wellbeing_coded',
                   '+', 'bmi_cat_coded', '+', 'antibiotic_history_coded', 
                   '+', '(1 | thdmi_cohort)'), 
                 collapse = " ")
    ex <- lmer(formula(form), data=no_nans)
    out <- summary(ex)
    r2 <- r.squaredGLMM(ex)
    if (out$coefficients[20] < 0) {
      marginal_r2 <- -1* r2[1]
      conditional_r2 <- -1* r2[2] 
    } else {
      marginal_r2 <- r2[1]
      conditional_r2 <- r2[2]
    }
    p_value_var <- out$coefficients[26]
    out_df <- rbind(out_df, c(variables_prev[x], prevotella_strains[y], 
                              p_value_var, marginal_r2, conditional_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'marginal_r2', 'conditional_r2')
write.csv(out_df, file='prevotella_full_log_ratio_analysis.csv')

# all faecalibacterium 
out_df <- data.frame()
for(x in 1:length(variables_prev)){
  # remove values that are NA 
  no_nans = faec_lmer %>% drop_na(variables_prev[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(faec_cols)){
    form = paste(c(faec_cols[y], '~', variables_prev[x], '+', 
                   'covid_level_of_wellbeing_coded', '+', 'sex',
                   '+', 'bmi_cat_coded', '+', 'antibiotic_history_coded', 
                   '+', '(1 | thdmi_cohort)'), 
                 collapse = " ")
    ex <- lmer(formula(form), data=no_nans)
    out <- summary(ex)
    r2 <- r.squaredGLMM(ex)
    if (out$coefficients[20] < 0) {
      marginal_r2 <- -1* r2[1]
      conditional_r2 <- -1* r2[2] 
    } else {
      marginal_r2 <- r2[1]
      conditional_r2 <- r2[2]
    }
    p_value_var <- out$coefficients[26]
    out_df <- rbind(out_df, c(variables_prev[x], faec_cols[y], 
                              p_value_var, marginal_r2, conditional_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'marginal_r2', 'conditional_r2')
write.csv(out_df, file='faecalibacterium_full_log_ratio_analysis.csv')



# prev_high 
out_df <- data.frame()
for(x in 1:length(variables_prev)){
  # remove values that are NA 
  no_nans = prev_high %>% drop_na(variables_prev[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(prevotella_strains)){
    form = paste(c(prevotella_strains[y], '~', variables_prev[x], '+', 'sex', 
                   '+','covid_level_of_wellbeing_coded',
                   '+', 'bmi_cat_coded', '+', 'antibiotic_history_coded', 
                   '+', '(1 | thdmi_cohort)'), 
                 collapse = " ")
    ex <- lmer(formula(form), data=no_nans)
    out <- summary(ex)
    r2 <- r.squaredGLMM(ex)
    if (out$coefficients[20] < 0) {
      marginal_r2 <- -1* r2[1]
      conditional_r2 <- -1* r2[2] 
    } else {
      marginal_r2 <- r2[1]
      conditional_r2 <- r2[2]
    }
    p_value_var <- out$coefficients[26]
    out_df <- rbind(out_df, c(variables_prev[x], prevotella_strains[y], 
                              p_value_var, marginal_r2, conditional_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'marginal_r2', 'conditional_r2')
write.csv(out_df, file='prevotella_high_log_ratio_analysis.csv')

# prev_low 
out_df <- data.frame()
for(x in 1:length(variables_prev)){
  # remove values that are NA 
  no_nans = prev_low %>% drop_na(variables_prev[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(prevotella_strains)){
    form = paste(c(prevotella_strains[y], '~', variables_prev[x], '+', 'sex',
                   '+', 'covid_level_of_wellbeing_coded',
                   '+', 'bmi_cat_coded', '+', 'antibiotic_history_coded', 
                   '+', '(1 | thdmi_cohort)'), 
                 collapse = " ")
    ex <- lmer(formula(form), data=no_nans)
    out <- summary(ex)
    r2 <- r.squaredGLMM(ex)
    if (out$coefficients[20] < 0) {
      marginal_r2 <- -1* r2[1]
      conditional_r2 <- -1* r2[2] 
    } else {
      marginal_r2 <- r2[1]
      conditional_r2 <- r2[2]
    }
    p_value_var <- out$coefficients[26]
    out_df <- rbind(out_df, c(variables_prev[x], prevotella_strains[y], 
                              p_value_var, marginal_r2, conditional_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'multiple_r2', 'adjusted_r2')
write.csv(out_df, file='prevotella_low_log_ratio_analysis.csv')

#variables_prev <- c(variables_prev, 'log_Prevotella_to_Bacteroides')

# prev_mx 
out_df <- data.frame()
for(x in 1:2){
  # remove values that are NA 
  no_nans = prev_mx %>% drop_na(variables_prev[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(prevotella_strains)){
    form = paste(c(prevotella_strains[y], '~', variables_prev[x], '+', 
                   'covid_level_of_wellbeing_coded', '+', 'sex',
                   '+', 'bmi_cat_coded', '+', 'antibiotic_history_coded'), 
                 collapse = " ")
    ex <- lm(formula(form), data=no_nans)
    out <- summary(ex)
    if (out$coefficients[14] < 0) {
      multiple_r2 <- -1* out$r.squared[1]
      adjusted_r2 <- -1* out$adj.r.squared[1]
    } else {
      multiple_r2 <- out$r.squared[1]
      adjusted_r2 <- out$adj.r.squared[1]
    }
    p_value_var <- out$coefficients[20]
    out_df <- rbind(out_df, c(variables_prev[x], prevotella_strains[y], 
                              p_value_var, multiple_r2, adjusted_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'multiple_r2', 'adjusted_r2')
write.csv(out_df, file='prevotella_mexico_log_ratio_analysis.csv')

# faecalibacterium mx
out_df <- data.frame()
for(x in 1:length(variables_prev)){
  # remove values that are NA 
  no_nans = faec_mx %>% drop_na(variables_prev[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(faec_cols)){
    form = paste(c(faec_cols[y], '~', variables_prev[x], '+', 'sex', '+',
                   'covid_level_of_wellbeing_coded',
                   '+', 'bmi_cat_coded', '+', 'antibiotic_history_coded'), 
                 collapse = " ")
    ex <- lm(formula(form), data=no_nans)
    out <- summary(ex)
    if (out$coefficients[14] < 0) {
      multiple_r2 <- -1* out$r.squared[1]
      adjusted_r2 <- -1* out$adj.r.squared[1]
    } else {
      multiple_r2 <- out$r.squared[1]
      adjusted_r2 <- out$adj.r.squared[1]
    }
    p_value_var <- out$coefficients[20]
    out_df <- rbind(out_df, c(variables_prev[x], faec_cols[y], 
                              p_value_var, multiple_r2, adjusted_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'multiple_r2', 'adjusted_r2')
write.csv(out_df, file='faecalibacterium_mexico_log_ratio_analysis.csv')

# prev_uk 
out_df <- data.frame()
for(x in 1:length(variables_prev)){
  # remove values that are NA 
  no_nans = prev_uk %>% drop_na(variables_prev[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(prevotella_strains)){
    form = paste(c(prevotella_strains[y], '~', variables_prev[x], '+', 
                   'covid_level_of_wellbeing_coded', '+', 'sex', 
                   '+', 'bmi_cat_coded', '+', 'antibiotic_history_coded'), 
                 collapse = " ")
    ex <- lm(formula(form), data=no_nans)
    out <- summary(ex)
    if (out$coefficients[14] < 0) {
      multiple_r2 <- -1* out$r.squared[1]
      adjusted_r2 <- -1* out$adj.r.squared[1]
    } else {
      multiple_r2 <- out$r.squared[1]
      adjusted_r2 <- out$adj.r.squared[1]
    }
    p_value_var <- out$coefficients[20]
    out_df <- rbind(out_df, c(variables_prev[x], prevotella_strains[y], 
                              p_value_var, multiple_r2, adjusted_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'multiple_r2', 'adjusted_r2')
write.csv(out_df, file='prevotella_uk_log_ratio_analysis.csv')

# faecalibacterium uk
out_df <- data.frame()
for(x in 1:length(variables_prev)){
  # remove values that are NA 
  no_nans = faec_uk %>% drop_na(variables_prev[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(faec_cols)){
    form = paste(c(faec_cols[y], '~', variables_prev[x], '+', 
                   'covid_level_of_wellbeing_coded', '+', 'sex', 
                   '+', 'bmi_cat_coded', '+', 'antibiotic_history_coded'), 
                 collapse = " ")
    ex <- lm(formula(form), data=no_nans)
    out <- summary(ex)
    if (out$coefficients[14] < 0) {
      multiple_r2 <- -1* out$r.squared[1]
      adjusted_r2 <- -1* out$adj.r.squared[1]
    } else {
      multiple_r2 <- out$r.squared[1]
      adjusted_r2 <- out$adj.r.squared[1]
    }
    p_value_var <- out$coefficients[20]
    out_df <- rbind(out_df, c(variables_prev[x], faec_cols[y], 
                              p_value_var, multiple_r2, adjusted_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'multiple_r2', 'adjusted_r2')
write.csv(out_df, file='faecalibacterium_uk_log_ratio_analysis.csv')

# prev_us
out_df <- data.frame()
for(x in 1:length(variables_prev)){
  # remove values that are NA 
  no_nans = prev_us %>% drop_na(variables_prev[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(prevotella_strains)){
    form = paste(c(prevotella_strains[y], '~', variables_prev[x], '+', 
                   'covid_level_of_wellbeing_coded', '+', 'sex', 
                   '+', 'bmi_cat_coded', '+', 'antibiotic_history_coded'), 
                 collapse = " ")
    ex <- lm(formula(form), data=no_nans)
    out <- summary(ex)
    if (out$coefficients[14] < 0) {
      multiple_r2 <- -1* out$r.squared[1]
      adjusted_r2 <- -1* out$adj.r.squared[1]
    } else {
      multiple_r2 <- out$r.squared[1]
      adjusted_r2 <- out$adj.r.squared[1]
    }
    p_value_var <- out$coefficients[20]
    out_df <- rbind(out_df, c(variables_prev[x], prevotella_strains[y], 
                              p_value_var, multiple_r2, adjusted_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'multiple_r2', 'adjusted_r2')
write.csv(out_df, file='prevotella_us_log_ratio_analysis.csv')

# faecalibacterium us
out_df <- data.frame()
for(x in 1:length(variables_prev)){
  # remove values that are NA 
  no_nans = faec_us %>% drop_na(variables_prev[x])
  
  # calculate shannon & faith's PD 
  for(y in 1:length(faec_cols)){
    form = paste(c(faec_cols[y], '~', variables_prev[x], '+', 'sex', '+', 
                   'covid_level_of_wellbeing_coded',
                   '+', 'bmi_cat_coded', '+', 'antibiotic_history_coded'), 
                 collapse = " ")
    ex <- lm(formula(form), data=no_nans)
    out <- summary(ex)
    if (out$coefficients[14] < 0) {
      multiple_r2 <- -1* out$r.squared[1]
      adjusted_r2 <- -1* out$adj.r.squared[1]
    } else {
      multiple_r2 <- out$r.squared[1]
      adjusted_r2 <- out$adj.r.squared[1]
    }
    p_value_var <- out$coefficients[20]
    out_df <- rbind(out_df, c(variables_prev[x], faec_cols[y], 
                              p_value_var, multiple_r2, adjusted_r2))
  }
  
}
colnames(out_df) <- c('variable', 'log_ratio', 'p_value', 'multiple_r2', 'adjusted_r2')
write.csv(out_df, file='faecalibacterium_us_log_ratio_analysis.csv')

