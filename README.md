# 3_Country_THDMI  
**Code for manuscript:** *“A three-country analysis of the gut microbiome indicates taxon associations with diet vary by taxon resolution and population”*  

## 📋 Overview  
This repository contains the code and supporting files used in the above-named manuscript (Khatib et al.; https://doi.org/10.1128/msystems.00544-25). The analysis investigates how gut microbial taxa relate to dietary variables across three countries, and examines how associations may vary by taxon resolution and population.

## 🧬 Aims  
- To explore taxon-diet associations in a multi-country cohort.  
- To assess how resolution (e.g., genus vs. strains) affects associations.  
- To highlight population (country)-specific differences in microbiome-diet relationships.  
- To provide reproducible code for the analyses reported in the manuscript.

## 📊 Data Availability
The data used in this study will be made available through:
- European Bioinformatics Institute (EBI) under accession number PRJEB11419
- Qiita under Study ID 10317

## 📂 Repository Structure  
```text

3_Country_THDMI/
├── Data_Cleaning/
│ ├── 0_calculating_coverages_micov-THDMI.sh ← Shell script to compute microbial coverages for THDMI samples
│ └── 1_sample_feature_filtering.ipynb ← Notebook for sample- and feature-level preprocessing and filtering
│
├── Figure1/
│ ├── adonis/ ← Alpha/beta diversity and PERMANOVA (adonis) pipelines for Figure 1
│ │ ├── 1_run_alpha_beta.sh ← Runs alpha and beta diversity calculations
│ │ ├── 2_run_adonis.sh ← Pipeline for running PERMANOVA analyses
│ │ ├── adonis.R ← R implementation of PERMANOVA and model setup
│ │ ├── adonis.py ← Python implementation of adonis-like testing
│ │ ├── adonis_and_variables.ipynb ← Notebook exploring adonis models and metadata variables
│ │ ├── adonis_run.py ← Main Python runner for PERMANOVA workflows
│ │ └── alpha_beta_calc.py ← Script to compute alpha/beta diversity metrics
│ ├── figure1c_radar_plots.ipynb ← Generates radar plots for Figure 1C
│ ├── figure1d_RF_Classifiers.ipynb ← Random forest classification models for Figure 1D
│ └── figure1d_rpca.sh ← Shell script to run RPCA for Figure 1D
│
├── Figure2/
│ ├── figure2a_instrain_diversity.R ← InStrain-based microdiversity calculations for Figure 2A
│ ├── figure2b_instrain_diversity_sensitivity_analysis.R ← Sensitivity analyses for microdiversity (Figure 2B)
│ ├── figure2c_1_log_ratio_analysis.ipynb ← Log-ratio analysis of taxa for Figure 2C
│ ├── figure2c_2_MAG_prevotella.ipynb ← Prevotella-focused MAG analysis for Figure 2C
│ ├── figure2c_3_lmer_alpha.R ← Linear mixed-effects models for MAGs for Figure 2C
│ └── main_diet_variables.txt ← Primary diet variable list used across analyses
│
├── Supplemental/
│ ├── figureS2a_b_log_ratio_analysis_clr.ipynb ← Supplemental CLR-based log-ratio analyses
│ ├── figureS2c_d_MAG_per_country.ipynb ← Supplemental MAG analyses stratified by country
│ └── lmer_alpha.R ← Supplemental linear mixed-effects model code
└── STORMS_Excel_1.03.xlsx ← STORMS metadata template used for data organization
 
