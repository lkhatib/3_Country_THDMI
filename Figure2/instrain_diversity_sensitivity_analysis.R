---
title: "instrain_diversity"
output: html_document
date: "2023-02-23"
---

```{r}
require(ggplot2)
require(dplyr)
require(cluster)
require(tidyr)
require(ComplexHeatmap)
require(RColorBrewer)
require(circlize)
require(grid)
require(glmnet)
require(data.table)
```

## Load metadata

```{r}

md_fp = '../consolidated_metadata_with_kcal_normalized_subset.tsv'
md_df = read.csv(md_fp, sep='\t', quote="")
colnames(md_df)[1] = 'Sample'

```

Target variables

```{r}
target_variables_fp = '../selected_variables.txt'
variables = read.table(target_variables_fp, sep='\t')$V1

variables %in% colnames(md_df)
```

```{r}
variable_categories <- read.table('../VariableCategory.txt', sep='\t', header = TRUE)

variable_subcat <- read.table('../VariableSubCategory.tsv', sep='\t', header=TRUE)
variable_subcat[variable_subcat$panel_name == '' & !is.na(variable_subcat$panel_name), 'panel_name'] = NA
variable_subcat <- variable_subcat %>% drop_na(panel_name)


variable_cats <- merge(x = variable_categories, 
      y = variable_subcat,
      by.x = 'Variable',
      by.y = 'variable_name',
      all.x = TRUE)

variable_cats[is.na(variable_cats$panel_name) & !is.na(variable_cats$Category) & variable_cats$Category != 'Nutrients', 'panel_name'] = 
variable_cats[is.na(variable_cats$panel_name) & !is.na(variable_cats$Category) & variable_cats$Category != 'Nutrients', 'Category']

variable_cats[variable_cats$panel_name == 'Food_Groups' & !is.na(variable_cats$panel_name), 'panel_name'] = 'food groups'
variable_cats[variable_cats$panel_name == 'Overall_Patterns' & !is.na(variable_cats$panel_name), 'panel_name'] = 'overall diet'
variable_cats[variable_cats$panel_name == 'Food_Groups' & !is.na(variable_cats$panel_name), 'panel_name'] = 'food groups'
variable_cats[variable_cats$panel_name =='Fermented' & !is.na(variable_cats$panel_name), 'panel_name'] = 'fermented'

final_variable_cats <- variable_cats[variable_cats$panel_name != 'remove' & !is.na(variable_cats$panel_name), c('Variable', 'panel_name')]

```

## Load instrain output file

```{r}

#instrain_df_in <- read.csv('../instrain_genome_df.csv')
instrain_df_in <- read.csv('../instrain_genome_df.high_med.csv')

instrain_df <- subset(instrain_df_in, sample %in% md_df$Sample)

instrain_df$Genus_orig <- instrain_df$Genus
instrain_df$Genus <- gsub("-", "_", instrain_df$Genus_orig)
```

## Strain diversity vs various variables

Lasso regression to do variable selection

```{r}

# get variables

predictors <- as.matrix(md_df[, variables])

# graft nucleotide diversity per sGB onto metadata table

nucl_values <- instrain_df[,c('genome',
               'nucl_diversity_rarefied',
               'sample')] %>% pivot_wider(names_from = genome,
                                          values_from = nucl_diversity_rarefied)

nucl_data <- merge(x = md_df[, c(c('Sample'), variables)],
      y = nucl_values,
      by.x = 'Sample',
      by.y = 'sample')

```

```{r}

exclude <- c('DP_All',
                             'DP_Core',
                             'diet_type_coded',
                             'specialized_diet_exclude_dairy',
                             'specialized_diet_exclude_refined_sugars')
sub_variables <- variables[-which(variables %in% exclude )]
```

## Genus

```{r}

taxon = 'Genus'

run_lasso <- function(instrain_df, md_df, variables, sub_variables, force_variables, taxon) {
  
  
nucl_values_genus <- instrain_df[,c('genome',
               taxon,
               'nucl_diversity_rarefied',
               'sample')] %>% pivot_wider(names_from = {{taxon}},
                                          values_from = nucl_diversity_rarefied)


nucl_data_genus <- merge(x = md_df[, c(c('Sample', 'thdmi_cohort'), variables)],
      y = nucl_values_genus,
      by.x = 'Sample',
      by.y = 'sample')

# create empty dataframe for pvals

pvals_taxon = data.frame(matrix(nrow=length(unique(instrain_df$Genus)), ncol=(length(sub_variables) + 1)))
colnames(pvals_taxon) = c('thdmi_cohort', sub_variables)
pvals_taxon$thdmi_cohortUK = NA
pvals_taxon$thdmi_cohortUS = NA
taxon_list = unique(instrain_df[,taxon])
rownames(pvals_taxon) = taxon_list
tvals_taxon <- pvals_taxon

for(target in taxon_list) {
  

  subset = drop_na(nucl_data_genus[,c(c(target, 'thdmi_cohort'), sub_variables)])
  # Skip if too few rows (e.g. < 10), which glmnet might fail on
  if (nrow(subset) < 13) {
    message(paste("Skipping", target, "- too few samples"))
    next
  }

  #use 5-fold cross validation to pick lambda
  tryCatch(expr = {
        print(paste("Running LASSO for:", target))  # Debugging print
        print(dim(subset))  # Print dimensions of subset to verify if it's empty
        cv_lasso_fit <- cv.glmnet(x = as.matrix(subset[, c(sub_variables)]),
                                  y = subset[,target],
                                  alpha = 1,
                                  nfolds = 5)
        lasso_coef = predict(cv_lasso_fit, type = "coefficients", s = cv_lasso_fit$lambda.1se)
        if(length(lasso_coef[lasso_coef != 0]) > 1) {
          print(target)
          plot(cv_lasso_fit)
          mod1.1se <- which(cv_lasso_fit$glmnet.fit$beta[, which(cv_lasso_fit$lambda == cv_lasso_fit$lambda.1se)] != 0)

          pred <- paste(c(names(mod1.1se), force_variables), collapse = ' + ')
         
        }
        else {
          if( length(force_variables) > 0 ){
            pred <- paste(force_variables, collapse = ' + ')
          }
          else {
            pred <- "1"
          }
        }
        

        form <- as.formula(paste(target, pred, sep = " ~ "))
        cmod1.1se <- lm(form, subset)
        modsum <- summary(cmod1.1se)
        for(x in names(modsum$coefficients[,4])) {
          print(x)
          p = modsum$coefficients[x,4]
          t = modsum$coefficients[x,3]
          if (x == "(Intercept)") { next }
          #if (startsWith(x, 'thdmi_cohort')) {
          #  x = 'thdmi_cohort'
          #}
          tvals_taxon[target, x] = t
          if (p < 0.05) {
            pvals_taxon[target, x] = p
          }
          
         }
        },
        
        error = function(e){
            message('Caught an error!')
            print(e)
            
        },
        warning = function(w){
            #message('Caught an warning!')
            #print(w)
        },
        finally = {
            #message('All done, quitting.')
        })
  
  bestlam = cv_lasso_fit$lambda.min # Select lambda that minimizes training MSE
  
  lasso_coef = predict(cv_lasso_fit, type = "coefficients", s = cv_lasso_fit$lambda.1se) 
  

}

return(list(tvals_taxon, pvals_taxon))
}

```

```{r}

library(ggplot2)
library(dplyr)
library(tidyr)

# Function to process LASSO results and count significant p-values
count_significant <- function(pvals_df, threshold = 0.05) {
  pvals_df <- pvals_df[, -1, drop = FALSE]  # Drop cohort column
  pvals_df[is.na(pvals_df)] <- 1  # Replace NA with 1 (non-significant)
  
  # Count the number of significant p-values per genus
  significant_counts <- rowSums(pvals_df < threshold)
  
  # Convert to data frame
  result_df <- data.frame(Genus = rownames(pvals_df), Significant_Count = significant_counts)
  result_df <- result_df[result_df$Significant_Count != 0, ]
  return(result_df)
}

# Function to create bar plots
# Function to create bar plots with cohort-specific colors
plot_significant_counts <- function(results_df, cohort_name, cohort_colors) {
  
  # Select color based on cohort_name, with a default fallback
  fill_color <- cohort_colors[[cohort_name]]
  if (is.null(fill_color)) fill_color <- "gray50"  # default if name not found

  ggplot(results_df, aes(x = reorder(Genus, -Significant_Count), y = Significant_Count)) +
    geom_bar(stat = "identity", fill = fill_color) +
    theme_classic() +
    labs(
      title = paste("Significant Dietary Associations in", cohort_name),
      x = "Genus",
      y = "Number of Significant Correlations"
    ) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 11, angle = 90, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold")
    )
}
```

```{r}
lasso_correlations <- function(md_dict, force_variables, cohort_colors) {
  # Create an empty list to store results
 lasso_cohort_results <- list()
  
  # Run LASSO for each country and store the result
for (name in names(md_dict)) {
   print(paste("Running LASSO for cohort:", name))  # Debugging message
    lasso_cohort_results[[name]] <- run_lasso(instrain_df, md_dict[[name]], variables, sub_variables, force_variables, 'Genus')
 }
  
  
  # Create and store plots for each cohort
  significant_results <- list()
  plots <- list()
  
  for (name in names(lasso_cohort_results)) {
    print(paste("Processing results for:", name))
    
      # Extract p-values data frame
   pvals_df <- lasso_cohort_results[[name]][[2]]  # Second element is p-values
    
    # Count significant correlations
   significant_results[[name]] <- count_significant(pvals_df)
    
    # Generate bar plot
    plots[[name]] <- plot_significant_counts(significant_results[[name]], name, cohort_colors)
    
    # Print plot
    print(plots[[name]])
  }
  
  output_dir <- file.path(getwd(), "lasso_plots")  # Ensure full path
  
  # Create directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  
  for (name in names(plots)) {
    file_path <- file.path(output_dir, paste0("lasso_results_", name, ".png"))
    ggsave(file_path, plot = plots[[name]], width = 10, height = 6, dpi = 300)
    print(paste("Saved plot for", name, "to", file_path))
  }
  }
```

```{r}
# Create a dictionary of metadata subsets
md_dict <- list(
  "Full" = md_df,
#  "US" = md_df[md_df$thdmi_cohort == "US", ], 
  "UK" = md_df[md_df$thdmi_cohort == "UK", ],
  "Mexico" = md_df[md_df$thdmi_cohort == "Mexico", ]
)

  # Define a named vector of cohort colors
  cohort_colors <- c(
    "US" = "steelblue",
    "UK" = "orange",
    "Mexico" = "forestgreen",
    "Full" = "purple"
  )

lasso_correlations(md_dict, c('antibiotic_history_coded', 'covid_level_of_wellbeing_coded', 'bmi_cat_coded'))

```


```{r}

# Create a dictionary of metadata subsets
md_dict <- list(
#  "3.5" = md_df[md_df$antibiotic_history_coded == 3.5, ]
#  "15" = md_df[md_df$antibiotic_history_coded == 15, ],
#   "90" = md_df[md_df$antibiotic_history_coded == 90, ],
#  "180" = md_df[md_df$antibiotic_history_coded == 180, ],
  "365" = md_df[md_df$antibiotic_history_coded == 365, ]
)

  # Define a named vector of cohort colors
  cohort_colors <- c(
    "3.5" = "#a50f15",
    "15" = "#de2d26",
    "90" = "#fb6a4a",
    "180" = "#fcae91",
    "365" = "#fee5d9"
  )
lasso_correlations(md_dict, c('thdmi_cohort', 'covid_level_of_wellbeing_coded', 'bmi_cat_coded'), cohort_colors)
```


```{r}
# Create a dictionary of metadata subsets
md_dict <- list(
  "1" = md_df[md_df$covid_level_of_wellbeing_coded == 1, ],
  "2" = md_df[md_df$covid_level_of_wellbeing_coded == 2, ],
  "3" = md_df[md_df$covid_level_of_wellbeing_coded == 3, ],
  "4" = md_df[md_df$covid_level_of_wellbeing_coded == 4, ],
  "5" = md_df[md_df$covid_level_of_wellbeing_coded == 5, ]
)

# Define a named vector of cohort colors
cohort_colors <- c(
  "1" = "#2171b5",
  "2" = "#6baed6",
  "3" = "#9ecae1",
  "4" = "#c6dbef",
  "5" = "#deebf7"
)

lasso_correlations(md_dict, c('thdmi_cohort', 'antibiotic_history_coded', 'bmi_cat_coded'), cohort_colors)
```


```{r}

# Create a dictionary of metadata subsets
md_dict <- list(
#  "Underweight" = md_df[md_df$bmi_cat == "Underweight", ], 
  "Normal" = md_df[md_df$bmi_cat == "Normal", ]
#  "Overweight" = md_df[md_df$bmi_cat == "Overweight", ],
#  "Obese" = md_df[md_df$bmi_cat == "Obese", ]
)

# Define a named vector of cohort colors
cohort_colors <- c(
  "Underweight" = "#f2e5ff",
  "Normal" = "#d9b3ff",
  "Overweight" = "#bf80ff",
  "Obese" = "#9933ff"
)

lasso_correlations(md_dict, c('thdmi_cohort', 'antibiotic_history_coded', 'covid_level_of_wellbeing_coded'), cohort_colors)
```