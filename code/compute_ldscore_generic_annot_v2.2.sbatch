#!/bin/bash

#SBATCH --job-name=ldsc_generic_annot
#SBATCH --output=ldsc_generic_annot_%J.out
#SBATCH --error=ldsc_generic_annot_%J.err
#SBATCH --partition=broadwl
#SBATCH --mem=20G

source activate ldsc

module load R/4.2.0

prefix_annot=$1

dir_LDSC=/project2/xinhe/kevinluo/ldsc
dir_ldscore=${dir_LDSC}/annot/ldscores/

mkdir -p ${dir_ldscore}/${prefix_annot}
cd ${dir_ldscore}/${prefix_annot}

echo "annot: ${dir_ldscore}/${prefix_annot}"

for chrom in {1..22}
do
  echo ${chrom}

  Rscript ~/projects/analysis_pipelines/code/extract_baselineLDv2.2_generic_annot.R ${prefix_annot} ${dir_ldscore} ${chrom}

  echo "Computing LD scores with the annot file ${prefix_annot}.${chrom}.annot.gz"

  python $HOME/softwares/ldsc/ldsc.py \
  --l2 \
  --bfile ${dir_LDSC}/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
  --print-snps ${dir_LDSC}/LDSCORE/listHM3.txt \
  --ld-wind-cm 1 \
  --annot ${dir_ldscore}/${prefix_annot}/${prefix_annot}.${chrom}.annot.gz \
  --out ${dir_ldscore}/${prefix_annot}/${prefix_annot}.${chrom}

done
