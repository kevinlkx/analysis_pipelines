#!/bin/bash

#SBATCH --job-name=sldsc
#SBATCH --output=sldsc_%J.out
#SBATCH --error=sldsc_%J.err
#SBATCH --partition=broadwl
#SBATCH --mem=10G

dir_GWAS=$1
trait=$2
dir_sLDSC_output=$3

echo "GWAS trait: ${trait}"

dir_code=~/projects/analysis_pipelines/code/

sh ${dir_code}/sldsc_annot_baselineLD_separate.sbatch ${dir_GWAS} ${trait} CN_all_peaks.narrowPeak.cleaned.hg19.merged ${dir_sLDSC_output}

sh ${dir_code}/sldsc_annot_baselineLD_separate.sbatch ${dir_GWAS} ${trait} DN_all_peaks.narrowPeak.cleaned.hg19.merged ${dir_sLDSC_output}

sh ${dir_code}/sldsc_annot_baselineLD_separate.sbatch ${dir_GWAS} ${trait} GA_all_peaks.narrowPeak.cleaned.hg19.merged ${dir_sLDSC_output}

sh ${dir_code}/sldsc_annot_baselineLD_separate.sbatch ${dir_GWAS} ${trait} ips_all_peaks.narrowPeak.cleaned.hg19.merged ${dir_sLDSC_output}

sh ${dir_code}/sldsc_annot_baselineLD_separate.sbatch ${dir_GWAS} ${trait} NSC_all_peaks.narrowPeak.cleaned.hg19.merged ${dir_sLDSC_output}

sh ${dir_code}/sldsc_annot_baselineLD_separate.sbatch ${dir_GWAS} ${trait} ASoC_glut_anno_hg19 ${dir_sLDSC_output}

sh ${dir_code}/sldsc_annot_baselineLD_separate.sbatch ${dir_GWAS} ${trait} ASoC_npc_anno_hg19 ${dir_sLDSC_output}


