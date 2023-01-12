## selected generic annotations for LDSC baselineLDv2.2.
args <- commandArgs(trailingOnly = TRUE)
name_annot <- args[1]
dir_annot <- args[2]
chrom <- args[3]

# name_annot <- "baseline_gene_MAF_LD_histone_v2.2"
# dir_annot <- "/project2/xinhe/kevinluo/ldsc/annot/ldscores/"
# chrom <- 10

library(data.table)
library(dplyr)

dir_baselineLD <- "/project2/xinhe/kevinluo/ldsc/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores"
annot <- fread(paste0(dir_baselineLD, "/baselineLD.", chrom, ".annot.gz"), header = T)

if(name_annot == "baseline_gene_MAF_LD_v2.2"){
  annot_list <- c("CHR", "BP", "SNP", "CM", "base",
                  "Coding_UCSC", "Coding_UCSC.flanking.500",
                  "UTR_5_UCSC", "UTR_5_UCSC.flanking.500",
                  "UTR_3_UCSC", "UTR_3_UCSC.flanking.500",
                  "Intron_UCSC", "Intron_UCSC.flanking.500",
                  "Promoter_UCSC", "Promoter_UCSC.flanking.500", "PromoterFlanking_Hoffman", "PromoterFlanking_Hoffman.flanking.500",
                  "Human_Promoter_Villar", "Human_Promoter_Villar.flanking.500", "Human_Promoter_Villar_ExAC", "Human_Promoter_Villar_ExAC.flanking.500",
                  paste0("MAFbin", 1:10),
                  "MAF_Adj_Predicted_Allele_Age","MAF_Adj_LLD_AFR", "MAF_Adj_ASMC",
                  "Recomb_Rate_10kb", "Nucleotide_Diversity_10kb", "Backgrd_Selection_Stat", "CpG_Content_50kb")
}else if(name_annot == "baseline_gene_MAF_LD_histone_v2.2"){
  annot_list <- c("CHR", "BP", "SNP", "CM", "base",
                  "Coding_UCSC", "Coding_UCSC.flanking.500",
                  "UTR_5_UCSC", "UTR_5_UCSC.flanking.500",
                  "UTR_3_UCSC", "UTR_3_UCSC.flanking.500",
                  "Intron_UCSC", "Intron_UCSC.flanking.500",
                  "Promoter_UCSC", "Promoter_UCSC.flanking.500", "PromoterFlanking_Hoffman", "PromoterFlanking_Hoffman.flanking.500",
                  "Human_Promoter_Villar", "Human_Promoter_Villar.flanking.500", "Human_Promoter_Villar_ExAC", "Human_Promoter_Villar_ExAC.flanking.500",
                  paste0("MAFbin", 1:10),
                  "MAF_Adj_Predicted_Allele_Age","MAF_Adj_LLD_AFR", "MAF_Adj_ASMC",
                  "Recomb_Rate_10kb", "Nucleotide_Diversity_10kb", "Backgrd_Selection_Stat", "CpG_Content_50kb",
                  "H3K27ac_Hnisz", "H3K27ac_Hnisz.flanking.500", "H3K27ac_PGC2", "H3K27ac_PGC2.flanking.500",
                  "H3K4me1_peaks_Trynka", "H3K4me1_Trynka", "H3K4me1_Trynka.flanking.500",
                  "H3K4me3_peaks_Trynka", "H3K4me3_Trynka", "H3K4me3_Trynka.flanking.500",
                  "H3K9ac_peaks_Trynka", "H3K9ac_Trynka", "H3K9ac_Trynka.flanking.500")
}else if(name_annot == "baseline_gene_v2.2"){
  annot_list <- c("CHR", "BP", "SNP", "CM", "base",
                  "Coding_UCSC", "Coding_UCSC.flanking.500",
                  "UTR_5_UCSC", "UTR_5_UCSC.flanking.500",
                  "UTR_3_UCSC", "UTR_3_UCSC.flanking.500",
                  "Intron_UCSC", "Intron_UCSC.flanking.500",
                  "Promoter_UCSC", "Promoter_UCSC.flanking.500", "PromoterFlanking_Hoffman", "PromoterFlanking_Hoffman.flanking.500",
                  "Human_Promoter_Villar", "Human_Promoter_Villar.flanking.500", "Human_Promoter_Villar_ExAC", "Human_Promoter_Villar_ExAC.flanking.500")
}else{
  stop("Please specify the correct name_annot! \n")
}

stopifnot(all(annot_list %in% colnames(annot)))

cat("Extract annotations from : ", paste0(dir_baselineLD, "/baselineLD.", chrom, ".annot.gz"), "\n")
annot_included <- annot %>% dplyr::select(all_of(annot_list))
print(colnames(annot_included))

dir.create(paste0(dir_annot, "/", name_annot), showWarnings = F, recursive = T)
fwrite(annot_included, paste0(dir_annot, "/", name_annot, "/", name_annot, ".", chrom, ".annot.gz"), sep = "\t", col.names = T, row.names = F)

