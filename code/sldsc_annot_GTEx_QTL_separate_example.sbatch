#!/bin/bash

#SBATCH --job-name=sldsc_QTL
#SBATCH --output=sldsc_QTL_%J.out
#SBATCH --error=sldsc_QTL_%J.err
#SBATCH --partition=broadwl
#SBATCH --mem=10G

trait=$1
prefix_annot=$2

dir_GWAS=/project2/xinhe/kevinluo/GWAS/GWAS_summary_stats/GWAS_summary_alkesgroup/independent_sumstats/

source activate ldsc

dir_annot=/project2/xinhe/kevinluo/ldsc/results/sLDSC_Hormozdiari_NG2018/LDSC_QTL/Annot/

echo "GWAS trait: ${trait}"
echo "dir_annot: ${dir_annot}"
echo "Annot: ${prefix_annot}"

cd ${dir_GWAS}

## Adjusting for baselineLD_v1.1
dir_out=/project2/xinhe/kevinluo/ldsc/results/sLDSC_Hormozdiari_NG2018/LDSC_QTL/results_sLDSC/${trait}/baselineLDv1.1
mkdir -p ${dir_out}

python $HOME/softwares/ldsc/ldsc.py \
--h2 ${dir_GWAS}/${trait}.sumstats.gz \
--ref-ld-chr /project2/xinhe/kevinluo/ldsc/LDSCORE/1000G_Phase3_baselineLD_v1.1_ldscores/baselineLD.,${dir_annot}/${prefix_annot}/${prefix_annot}. \
--frqfile-chr /project2/xinhe/kevinluo/ldsc/LDSCORE/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr /project2/xinhe/kevinluo/ldsc/LDSCORE/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out ${dir_out}/${trait}_${prefix_annot}_baselineLDv1.1

echo "Results saved at: ${dir_out}"
