## make continuous annotation file from bim file and annotation bed file
args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]
bim_file <- args[2]
annot_file <- args[3]
thin_annot <- args[4]

if(is.na(thin_annot) | is.null(thin_annot)){
  thin_annot <- "full"
}

suppressPackageStartupMessages(library(GenomicRanges))

annot_bed <- read.table(bed_file, header = F, stringsAsFactors = F)
colnames(annot_bed) <- c("chr", "start", "end", "value")
annot_bed$chr <- gsub("chr", "", annot_bed$chr)
annot_bed$start <- annot_bed$start + 1

dir_annot_ldsc <- dirname(annot_file)
dir.create(dir_annot_ldsc, showWarnings = F, recursive = T)

# cat("Make LDSC continuous annot for: ", bed_file, "\n")

bim_data <- read.table(bim_file, header = F, stringsAsFactors = F)
colnames(bim_data) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

annot_ldsc <- bim_data[, c("CHR", "BP", "SNP", "CM")]
annot_ldsc.gr <- makeGRangesFromDataFrame(annot_ldsc, seqnames.field = "CHR", start.field = "BP", end.field = "BP", keep.extra.columns = T)

annot_bed.gr <- makeGRangesFromDataFrame(annot_bed, keep.extra.columns = T)

idx_overlaps <- data.frame(findOverlaps(query = annot_ldsc.gr, subject = annot_bed.gr))

annot_ldsc[, "ANNOT"] <- 0
annot_ldsc[idx_overlaps$queryHits, "ANNOT"] <- annot_bed$value[idx_overlaps$subjectHits]

if(thin_annot == "thin-annot"){
  write.table(data.frame(ANNOT = annot_ldsc[, "ANNOT"]), gzfile(annot_file), col.names = T, row.names = F, quote = F, sep = "\t")
}else{
  write.table(annot_ldsc, gzfile(annot_file), col.names = T, row.names = F, quote = F, sep = "\t")
}

cat("LDSC annotation file saved at: ", annot_file, "\n")

