#!/bin/bash

qiime gemelli rpca \
        --i-table /home/lakhatib/3country/final_scripts/data_cleaning/filtered_ft.qzaa \
        --o-biplot rpca_ordination.qza \
        --o-distance-matrix rpca_distance_matrix.qza
        
        
qiime emperor biplot \
    --i-biplot rpca_ordination.qza \
    --m-sample-metadata-file /home/lakhatib/3country/final_scripts/data_cleaning/subsetted_md.tsv \
    --o-visualization biplot.qzv \
    
    
qiime diversity beta-group-significance \
    --i-distance-matrix rpca_distance_matrix.qza \
    --m-metadata-file /home/lakhatib/3country/final_scripts/data_cleaning/subsetted_md.tsv \
    --m-metadata-column thdmi_cohort \
    --p-method permanova \
    --o-visualization cohort_significance.qzv