#!/bin/bash

#SBATCH --job-name=ldsc_annot
#SBATCH --output=ldsc_annot_%J.out
#SBATCH --error=ldsc_annot_%J.err
#SBATCH --partition=broadwl
#SBATCH --mem=10G

## Compute annotation-specific LD scores
## https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial

prefix_annot=$1
echo "prefix_annot: ${prefix_annot}"

dir_LDSC=/project2/xinhe/kevinluo/ldsc/

mkdir -p ${dir_LDSC}/annot/ldscores/${prefix_annot}

conda activate ldsc

module load R/3.5.1

cd ${dir_LDSC}/annot/ldscores/

for chrom in {1..22}
do
  echo ${chrom}

  ## Step 1: Creating an annot file
  echo "Make ldsc-friendly annotation files for ${prefix_annot}.bed"
  Rscript ~/projects/analysis_pipelines/code/make_ldsc_binary_annot.R \
  ${dir_LDSC}/annot/annot_bed/${prefix_annot}.bed \
  ${dir_LDSC}/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
  ${dir_LDSC}/annot/ldscores/${prefix_annot}/${prefix_annot}.${chrom}.annot.gz "full-annot"

  ## Step 2: Computing LD scores with an annot file
  echo "Computing LD scores with the annot file ${prefix_annot}.${chrom}.annot.gz"
  python $HOME/softwares/ldsc/ldsc.py \
  --l2 \
  --bfile ${dir_LDSC}/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
  --print-snps ${dir_LDSC}/LDSCORE/listHM3.txt \
  --ld-wind-cm 1 \
  --annot ${dir_LDSC}/annot/ldscores/${prefix_annot}/${prefix_annot}.${chrom}.annot.gz \
  --out ${dir_LDSC}/annot/ldscores/${prefix_annot}/${prefix_annot}.${chrom}

done
