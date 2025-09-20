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

instrain_df <- subset(instrain_df, sample %in% md_df$Sample)

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

```{r}
target = 'sGB_00223'

subset = drop_na(nucl_data[,c(c(target), sub_variables)])



#use 5-fold cross validation to pick lambda
cv_lasso_fit <- cv.glmnet(x = as.matrix(subset[, c(sub_variables)]),
                          y = subset[,target],
                          alpha = 1,
                          nfolds = 5)
bestlam = cv_lasso_fit$lambda.min # Select lamda that minimizes training MSE
plot(cv_lasso_fit)
```

```{r}

lasso_coef = predict(cv_lasso_fit, type = "coefficients", s = cv_lasso_fit$lambda.min) # Display coefficients using lambda chosen by CV
lasso_coef[lasso_coef != 0]

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


  #use 5-fold cross validation to pick lambda
  tryCatch(expr = {
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
          if( length(force_variables) > 1 ){
            pred <- paste(force_variables, '+')
          }
          else {
            pred <- force_variables
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


```{r}
nucl_data_genus <- merge(x = md_df[, c(c('Sample', 'thdmi_cohort'), variables)],
      y = nucl_values_genus,
      by.x = 'Sample',
      by.y = 'sample')

```


## clustered heatmap

```{r}

# An R function to save pheatmap figure into pdf
# This was copied from Stackflow: https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}
```

```{r}

vals_main <- run_lasso(instrain_df, md_df, variables, sub_variables, c('thdmi_cohort'), 'Genus')
```

```{r}
vals_force_vars <-  run_lasso(instrain_df,
                              md_df,
                              variables,
                              sub_variables,
                              c('thdmi_cohort',
                                'antibiotic_history_coded',
                                'covid_level_of_wellbeing_coded'),
                              'Genus')
```

```{r}

us_cohort <- md_df[md_df$thdmi_cohort == 'US',]
instrain_df_us <- instrain_df[instrain_df$sample %in% us_cohort$Sample,]

vals_us <- run_lasso(instrain_df_us,
                              us_cohort,
                              variables,
                              sub_variables,
                              c('thdmi_cohort',
                                'antibiotic_history_coded',
                                'covid_level_of_wellbeing_coded',
                                'DP_All'),
                              'Genus')
```

## Main heatmap analysis

Only forcing THDMI_cohort

```{r}

run_plot <- function(vals, plotname) {



# T values

# clip and scale values
tvals_scaled <- subset(vals[[1]], select = -c(thdmi_cohort, thdmi_cohortUS, thdmi_cohortUK))
tvals_scaled[is.na(tvals_scaled)] <- 0

#vals_scaled[vals_scaled > 10] <- 10

tvals_scaled <- t(tvals_scaled[rowSums(tvals_scaled) != 0,
                            colSums(tvals_scaled) != 0])


# P values

# clip and scale values
pvals_scaled <- -log10(subset(vals[[2]], select = -c(thdmi_cohort, thdmi_cohortUS, thdmi_cohortUK)))
pvals_scaled[is.na(pvals_scaled)] <- 0

pvals_scaled[pvals_scaled > 10] <- 10

# select only p values for selected t values
pvals_scaled <- t(pvals_scaled[rownames(pvals_scaled) %in% colnames(tvals_scaled),
                               colnames(pvals_scaled) %in% rownames(tvals_scaled)])


genus_phylum <- distinct(instrain_df[,c('Genus','Phylum')])
genus_phylum <- genus_phylum[genus_phylum$Genus %in% colnames(pvals_scaled),]
rownames(genus_phylum) <- genus_phylum$Genus

phyla <- genus_phylum[colnames(pvals_scaled),'Phylum']
names(phyla) <- colnames(pvals_scaled)
  
panel_names <-  unique(final_variable_cats$panel_name)
panel_names <- c("food groups",
                 "fermented",
                 "sweeteners",
                 "overall diet",
                 "covariates",
                 "major macronutrients",
                 "macronutrients",
                 "major micronutrients",
                 "micronutrients",
                 "phytochemicals")
panel_dfs_p <- vector("list", length(panel_names))
panel_dfs_t <- vector("list", length(panel_names))
panel_vars <- vector("list", length(panel_names))
panel_hts_p <- vector("list", length(panel_names))
panel_hts_t <- vector("list", length(panel_names))
panel_ht_str_p <- character(length = length(panel_names))
panel_ht_str_t <- character(length = length(panel_names))
names(panel_dfs_p) <- panel_names
names(panel_dfs_t) <- panel_names
names(panel_vars) <- panel_names
names(panel_hts_p) <- panel_names
names(panel_hts_t) <- panel_names

virids_col = colorRamp2(seq(0, 10, length.out = 100), viridis::viridis(n = 100))

redblu_col = colorRamp2(seq(-10, 10, length.out = 11), brewer.pal(n = 11, name = "RdBu"))

phyla_col = brewer.pal(n = length(unique(phyla)),
                       name = "Set3")
names(phyla_col) <- sort(unique(phyla))

i = 1
colclust = agnes(t(tvals_scaled))

for(panel in panel_names) {
  this_panel_vars = final_variable_cats[final_variable_cats$panel_name == panel, 'Variable']
  panel_vars[[panel]] = final_variable_cats[final_variable_cats$panel_name == panel, 'Variable']
  
  panel_dfs_p[[panel]] = pvals_scaled[this_panel_vars[this_panel_vars %in% rownames(pvals_scaled)],]
  panel_dfs_t[[panel]] = tvals_scaled[this_panel_vars[this_panel_vars %in% rownames(tvals_scaled)],]
  
  rowclust = diana(tvals_scaled[rownames(tvals_scaled) %in% final_variable_cats[final_variable_cats$panel_name == panel, 'Variable'],])
  
  panel_hts_p[[panel]] <- Heatmap(as.matrix(panel_dfs_p[[panel]]),
                                col = virids_col,
                                name = panel,
                                cluster_columns = colclust,
                                cluster_rows = rowclust,
                                row_title = '',
                                row_title_rot = 0,
                                column_title = "P values (-log10)",
                                column_title_gp = gpar(fontsize = 24),
                                show_row_dend = FALSE,
                                show_row_names = TRUE,
                                row_names_gp = gpar(fontsize = 7),
                                show_heatmap_legend = FALSE)
  panel_hts_t[[panel]] <- Heatmap(as.matrix(panel_dfs_t[[panel]]),
                                col = redblu_col,
                                name = panel,
                                cluster_columns = colclust,
                                cluster_rows = rowclust,
                                row_title = panel,
                                row_title_rot = 0,
                                column_title = "T values",
                                column_title_gp = gpar(fontsize = 24),
                                show_row_dend = TRUE,
                                show_row_names = FALSE,
                                row_names_gp = gpar(fontsize = 7),
                                show_heatmap_legend = FALSE)
  
  panel_ht_str_p[i] = paste("panel_hts_p[['", panel, "']]", sep="")
  panel_ht_str_t[i] = paste("panel_hts_t[['", panel, "']]", sep="")
  i = i + 1
}

ht_phyla <- Heatmap(t(phyla),
                    col=phyla_col,
                    height=unit(1, "cm"),
                    show_heatmap_legend = FALSE)

ht_str_p <- paste(panel_ht_str_p, collapse= " %v% ")
ht_p <- eval(parse(text=paste(ht_str_p, '%v% ht_phyla',  sep='')))
ht_str_t <- paste(panel_ht_str_t, collapse= " %v% ")
ht_t <- eval(parse(text=paste(ht_str_t, '%v% ht_phyla',  sep='')))

# legends
lgd_phy = Legend(labels = names(phyla_col), 
                legend_gp = gpar(fill = phyla_col),
                title = "Phyla")
lgd_p = Legend(col_fun = virids_col,
               title = "P values (-log10)", 
               at = c(0, 2, 4, 6, 8, 10),
               direction = "horizontal")
lgd_t = Legend(col_fun = redblu_col, 
               title = "T values",
               at = c(-10, 0, 10),
               direction = "horizontal")

pd = packLegend(lgd_t, lgd_p, lgd_phy, direction = "horizontal")


pdf(paste("./plots/nucdiv_heatmap_", plotname ,"_vars_p.pdf", sep=''), width=8, height=24)
draw(ht_p,
     heatmap_legend_side = "bottom")
dev.off()
pdf(paste("./plots/nucdiv_heatmap_", plotname ,"_vars_t.pdf", sep=''), width=8, height=24)
draw(ht_t,
     heatmap_legend_side = "bottom")
dev.off()

#combined plot

pdf(paste("./plots/nucdiv_heatmap_", plotname ,"_vars_combined.pdf", sep=''), width=16, height=24)

gb_p <- grid.grabExpr(draw(ht_p,
                      heatmap_legend_side = "bottom"))
gb_t <- grid.grabExpr(draw(ht_t, 
                      heatmap_legend_side = "bottom"))

gb_l <- grid.grabExpr(draw(pd))

grid.newpage()
pushViewport(viewport(x = 0.5, y = 1,
            width=unit(8, "inches"), 
            height = unit(22, "inches"),
            just = c("left", "top")))
grid.draw(gb_p)
popViewport()

pushViewport(viewport(x = 0.0, y = 1,
            width=unit(8, "inches"), 
            height = unit(22, "inches"),
            just = c("left", "top")))
grid.draw(gb_t)
popViewport()

pushViewport(viewport(x = 0.5,
                      y = 0,
                      width = unit(12, "inches"),
                      height = unit(2, "inches"),
                      just = c("center","bottom")))
grid.draw(gb_l)
popViewport()

dev.off()


return(list(panel_dfs_p, panel_dfs_t))
}
```

```{r}
run_plot_binary_p <- function(vals, plotname, threshold=0.001) {

  # T values
  tvals_scaled <- subset(vals[[1]], select = -c(thdmi_cohort, thdmi_cohortUS, thdmi_cohortUK))
  tvals_scaled[is.na(tvals_scaled)] <- 0
  tvals_scaled <- t(tvals_scaled[rowSums(tvals_scaled) != 0, colSums(tvals_scaled) != 0])

  # P values
  pvals_scaled <- -log10(subset(vals[[2]], select = -c(thdmi_cohort, thdmi_cohortUS, thdmi_cohortUK)))
  pvals_scaled[is.na(pvals_scaled)] <- 0
  pvals_scaled[pvals_scaled > 10] <- 10

  # Convert P values to binary (0 or 1) based on a threshold
  pvals_scaled <- ifelse(pvals_scaled > -log10(threshold), 1, 0)

  # select only p values for selected t values
  pvals_scaled <- t(pvals_scaled[rownames(pvals_scaled) %in% colnames(tvals_scaled),
                                 colnames(pvals_scaled) %in% rownames(tvals_scaled)])

  genus_phylum <- distinct(instrain_df[,c('Genus','Phylum')])
  genus_phylum <- genus_phylum[genus_phylum$Genus %in% colnames(pvals_scaled),]
  rownames(genus_phylum) <- genus_phylum$Genus

  phyla <- genus_phylum[colnames(pvals_scaled),'Phylum']
  names(phyla) <- colnames(pvals_scaled)

  panel_names <- c("food groups",
                   "fermented",
                   "sweeteners",
                   "overall diet",
                   "covariates",
                   "major macronutrients",
                   "macronutrients",
                   "major micronutrients",
                   "micronutrients",
                   "phytochemicals")

  panel_dfs_p <- vector("list", length(panel_names))
  panel_dfs_t <- vector("list", length(panel_names))
  panel_vars <- vector("list", length(panel_names))
  panel_hts_p <- vector("list", length(panel_names))
  panel_hts_t <- vector("list", length(panel_names))
  panel_ht_str_p <- character(length = length(panel_names))
  panel_ht_str_t <- character(length = length(panel_names))
  names(panel_dfs_p) <- panel_names
  names(panel_dfs_t) <- panel_names
  names(panel_vars) <- panel_names
  names(panel_hts_p) <- panel_names
  names(panel_hts_t) <- panel_names

  # Adjust color scale to be binary for P values
  binary_col = colorRamp2(c(0, 1), viridis::viridis(n = 2))

  redblu_col = colorRamp2(seq(-10, 10, length.out = 11), brewer.pal(n = 11, name = "RdBu"))

  phyla_col = brewer.pal(n = length(unique(phyla)), name = "Set3")
  names(phyla_col) <- sort(unique(phyla))

  i = 1
  colclust = agnes(t(tvals_scaled))

  for(panel in panel_names) {
    this_panel_vars = final_variable_cats[final_variable_cats$panel_name == panel, 'Variable']
    panel_vars[[panel]] = final_variable_cats[final_variable_cats$panel_name == panel, 'Variable']

    panel_dfs_p[[panel]] = pvals_scaled[this_panel_vars[this_panel_vars %in% rownames(pvals_scaled)],]
    panel_dfs_t[[panel]] = tvals_scaled[this_panel_vars[this_panel_vars %in% rownames(tvals_scaled)],]

    rowclust = diana(tvals_scaled[rownames(tvals_scaled) %in% final_variable_cats[final_variable_cats$panel_name == panel, 'Variable'],])

    # Set fixed width for both heatmaps
    fixed_width <- unit(13, "cm")

    panel_hts_p[[panel]] <- Heatmap(as.matrix(panel_dfs_p[[panel]]),
                                  col = binary_col,  # Use binary color scale
                                  name = panel,
                                  cluster_columns = colclust,
                                  cluster_rows = rowclust,
                                  row_title = '',
                                  row_title_rot = 0,
                                  column_title = "P values",
                                  column_title_gp = gpar(fontsize = 24),
                                  show_row_dend = FALSE,
                                  show_row_names = FALSE,  # Hide row names (covariate labels)
                                  row_names_gp = gpar(fontsize = 7),
                                  width = fixed_width,  # Set fixed width for consistent layout
                                  show_heatmap_legend = FALSE)  # Hide legend

    panel_hts_t[[panel]] <- Heatmap(as.matrix(panel_dfs_t[[panel]]),
                                  col = redblu_col,
                                  name = panel,
                                  cluster_columns = colclust,
                                  cluster_rows = rowclust,
                                  row_title = panel,
                                  row_title_rot = 0,
                                  column_title = "t statistic",
                                  column_title_gp = gpar(fontsize = 24),
                                  show_row_dend = TRUE,
                                  show_row_names = FALSE,  # Hide covariate labels
                                  row_names_gp = gpar(fontsize = 7),
                                  width = fixed_width,  # Set fixed width for consistent layout
                                  show_heatmap_legend = FALSE)  # Keep the T values legend

    panel_ht_str_p[i] = paste("panel_hts_p[['", panel, "']]", sep="")
    panel_ht_str_t[i] = paste("panel_hts_t[['", panel, "']]", sep="")
    i = i + 1
  }

  ht_phyla <- Heatmap(t(phyla),
                      col=phyla_col,
                      height=unit(1, "cm"),
                      show_heatmap_legend = FALSE)

  ht_str_p <- paste(panel_ht_str_p, collapse= " %v% ")
  ht_p <- eval(parse(text=paste(ht_str_p, '%v% ht_phyla',  sep='')))
  ht_str_t <- paste(panel_ht_str_t, collapse= " %v% ")
  ht_t <- eval(parse(text=paste(ht_str_t, '%v% ht_phyla',  sep='')))

  # Reduce padding between heatmaps
  ht_layout <- HeatmapAnnotation(padding = unit(c(0, 0, 0, 0), "mm"))  # Adjust padding (left, top, right, bottom)

  # legends
  lgd_phy = Legend(labels = names(phyla_col), 
                   legend_gp = gpar(fill = phyla_col),
                   title = "Phyla")
  lgd_t = Legend(col_fun = redblu_col, 
                 title = "T values",
                 at = c(-10, 0, 10),
                 direction = "horizontal")
  
  pd = packLegend(lgd_t, lgd_phy, direction = "horizontal")
  
  
  # Combine remaining legends (T values and phyla)
  lgd_phy = Legend(labels = names(phyla_col), legend_gp = gpar(fill = phyla_col), title = "Phyla")
  lgd_t = Legend(col_fun = redblu_col, title = "T values", at = c(-10, 0, 10), direction = "horizontal")
  pd = packLegend(lgd_t, lgd_phy, direction = "horizontal")

  # combined plot
  pdf(paste("./plots/nucdiv_heatmap_", plotname ,"_vars_combined.binary_p.pdf", sep=''), width=16, height=24)

  gb_p <- grid.grabExpr(draw(ht_p,
                             heatmap_legend_side = "bottom"))
  gb_t <- grid.grabExpr(draw(ht_t, 
                             heatmap_legend_side = "bottom"))
  
  gb_l <- grid.grabExpr(draw(pd))
  
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 1,
                        width=unit(8, "inches"), 
                        height = unit(22, "inches"),
                        just = c("left", "top")))
  grid.draw(gb_p)
  popViewport()
  
  pushViewport(viewport(x = 0.0, y = 1,
                        width=unit(8, "inches"), 
                        height = unit(22, "inches"),
                        just = c("left", "top")))
  grid.draw(gb_t)
  popViewport()
  
  pushViewport(viewport(x = 0.5,
                        y = 0,
                        width = unit(12, "inches"),
                        height = unit(2, "inches"),
                        just = c("center","bottom")))
  grid.draw(gb_l)
  popViewport()
  
  dev.off()

  return(list(panel_dfs_p, panel_dfs_t))
}

```

```{r}

vals_main_output <- run_plot(vals_main, 'cohort')
vals_main_output_2 <- run_plot_binary_p(vals_main, 'cohort', threshold=0.01)

run_plot(vals_force_vars, 'cohort_covid_abx')

run_plot(vals_us, 'us_cohort')
```

