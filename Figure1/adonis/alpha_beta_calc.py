import qiime2 as q2
import numpy as np
import pandas as pd
import biom
from qiime2.plugins import diversity, gemelli, feature_table

ft = q2.Artifact.load('/home/lakhatib/3country/final_scripts/data/filtered_ft.qza')

tree = q2.Artifact.import_data('Phylogeny[Rooted]', '/projects/wol/qiyun/wol2/phylogeny/tree.nwk') 

# alpha diversity (shannon, faith's PD)
alpha_results = pd.read_csv('/home/lakhatib/3country/final_scripts/data/md_for_adonis.tsv', sep='\t')

ft = feature_table.methods.rarefy(table=ft, sampling_depth=1000000).rarefied_table

# alpha diversity (shannon) 
shannon = diversity.pipelines.alpha(table=ft, metric='shannon').alpha_diversity
alpha_results['shannon'] = shannon.view(pd.Series)

# beta diversity (rpca) 
rpca = gemelli.methods.rpca(table=ft).distance_matrix
rpca.save('results/rpca_dm.qza') 

# calculate faith's pd too
faith = diversity.pipelines.alpha_phylogenetic(table=ft, phylogeny=tree, metric='faith_pd').alpha_diversity
alpha_results['faith_pd'] = faith.view(pd.Series)

# beta diversity (uwuf, wuf, phylo-rpca)
uwuf = diversity.pipelines.beta_phylogenetic(table=ft, phylogeny=tree, metric='unweighted_unifrac').distance_matrix
uwuf.save('results/uwuf_dm.qza')
wuf = diversity.pipelines.beta_phylogenetic(table=ft, phylogeny=tree, metric='weighted_unifrac').distance_matrix
wuf.save('results/wuf_dm.qza')
prpca = gemelli.methods.phylogenetic_rpca_without_taxonomy(table=ft, phylogeny=tree).distance_matrix
prpca.save('results/prpca_dm.qza')

alpha_results.to_csv('results/alpha_results_3country.tsv', sep='\t') 