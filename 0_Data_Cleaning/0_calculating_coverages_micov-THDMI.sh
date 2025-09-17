#!/bin/bash
#SBATCH --chdir=/projects/thdmi/pangenome_filtered/coverages
#SBATCH --output=/projects/thdmi/pangenome_filtered/coverages/slurm_out/calculating_coverages.out
#SBATCH --time 2:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user="lakhatib@ucsd.edu"
#SBATCH --mem 120G
#SBATCH -N 1
#SBATCH -c 8

#source activate micov

micov qiita-coverage \
    --lengths /projects/wol/qiyun/wol2/genomes/length.map \
    --output /projects/thdsmi/pangenome_filtered/coverages/coverages \
    --qiita-coverages /qmounts/qiita_data/BIOM/198496/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/198709/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/198544/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/198506/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/198737/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/197725/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/197270/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/197129/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/197167/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/198637/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/197644/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/197653/coverages.tgz \
    --qiita-coverages /qmounts/qiita_data/BIOM/198415/coverages.tgz \
    --samples-to-keep /projects/thdmi/pangenome_filtered/thdmi_metadata.tsv








