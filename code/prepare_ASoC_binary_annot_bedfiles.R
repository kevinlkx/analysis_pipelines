## prepare ASoC annotations in BED format for LDSC analysis

dir_annot_bed <- "/project2/xinhe/kevinluo/ldsc/annot/annot_bed/"

### ASoC_glut_anno_hg19
annot_name <- "ASoC_glut_anno_hg19"
annot_filename <- "ASoC_glut_anno_hg19.txt"
ASoC_annot <- read.table(paste0(dir_annot_bed, "/", annot_filename), header = F, stringsAsFactors = F)
colnames(ASoC_annot) <- c("chr", "SNP_POS", "annot")

ASoC_sig <- ASoC_annot[ASoC_annot$annot == 1, ]
ASoC_sig.bed <- data.frame(chr = ASoC_sig$chr, start = ASoC_sig$SNP_POS - 1, end = ASoC_sig$SNP_POS)
ASoC_sig.bed$chr <- factor(ASoC_sig.bed$chr, levels = paste0("chr", 1:22))
ASoC_sig.bed <- ASoC_sig.bed[order(ASoC_sig.bed$chr, ASoC_sig.bed$start), ]
ASoC_sig.bed <- unique(ASoC_sig.bed)
cat(nrow(ASoC_sig.bed), "SNPs with annot \n")
write.table(ASoC_sig.bed, paste0(dir_annot_bed, "/", annot_name, ".bed"), sep = "\t", col.names = F, row.names = F, quote = F)

### ASoC_npc_anno_hg19
annot_name <- "ASoC_npc_anno_hg19"
annot_filename <- "ASoC_npc_anno_hg19.txt"
ASoC_annot <- read.table(paste0(dir_annot_bed, "/", annot_filename), header = F, stringsAsFactors = F)
colnames(ASoC_annot) <- c("chr", "SNP_POS", "annot")

ASoC_sig <- ASoC_annot[ASoC_annot$annot == 1, ]
ASoC_sig.bed <- data.frame(chr = ASoC_sig$chr, start = ASoC_sig$SNP_POS - 1, end = ASoC_sig$SNP_POS)
ASoC_sig.bed$chr <- factor(ASoC_sig.bed$chr, levels = paste0("chr", 1:22))
ASoC_sig.bed <- ASoC_sig.bed[order(ASoC_sig.bed$chr, ASoC_sig.bed$start), ]
ASoC_sig.bed <- unique(ASoC_sig.bed)
cat(nrow(ASoC_sig.bed), "SNPs with annot \n")
write.table(ASoC_sig.bed, paste0(dir_annot_bed, "/", annot_name, ".bed"), sep = "\t", col.names = F, row.names = F, quote = F)

cat("annotation bed files saved at:", dir_annot_bed, "\n")
