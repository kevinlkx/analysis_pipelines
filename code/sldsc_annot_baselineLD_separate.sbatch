#!/bin/bash

#SBATCH --job-name=sldsc
#SBATCH --output=sldsc_%J.out
#SBATCH --error=sldsc_%J.err
#SBATCH --partition=broadwl
#SBATCH --mem=10G

dir_GWAS=$1
trait=$2
prefix_annot=$3
dir_sLDSC_output=$4

dir_LDSC=/project2/xinhe/kevinluo/ldsc
dir_ldsc_annot=/project2/xinhe/kevinluo/ldsc/annot/ldscores
dir_baselineLD=/project2/xinhe/kevinluo/ldsc/LDSCORE/1000G_Phase3_baselineLD_v1.1_ldscores

source activate ldsc

echo "GWAS trait: ${trait}"

dir_out=${dir_sLDSC_output}/${trait}/baselineLDv1.1
mkdir -p ${dir_out}

python $HOME/softwares/ldsc/ldsc.py \
--h2 ${dir_GWAS}/${trait}.sumstats.gz \
--ref-ld-chr ${dir_baselineLD}/baselineLD.,${dir_ldsc_annot}/${prefix_annot}/${prefix_annot}. \
--frqfile-chr ${dir_LDSC}/LDSCORE/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${dir_LDSC}/LDSCORE/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out ${dir_out}/${trait}_${prefix_annot}_baselineLDv1.1

