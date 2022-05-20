#!/usr/bin/env Rscript
library(mapgen)
library(doParallel)

args <- commandArgs(trailingOnly=TRUE)

EUR_LD_1KG <- '/project2/xinhe/1kg/bigsnpr/EUR_variable_1kg.rds'
bigSNP <- bigsnpr::snp_attach(EUR_LD_1KG)

data.dir <- '/project2/gca/Heart_Atlas/reorganized_data/example_data'
sumstats.sigloci <- readRDS(paste0(data.dir, '/finemapping/sumstats.sigloci.rds'))
locus_list <- unique(sumstats.sigloci$locus)
z <- locus_list[as.integer(args[1])]

outfile <- paste0(data.dir, '/finemapping/AF_finemapping_result_torusprior_', z,'.rds')
if(!file.exists(outfile)){
  sumstats.locus <- sumstats.sigloci[sumstats.sigloci$locus == z, ]
  finemap_res <- run_finemapping(sumstats.locus, bigSNP, priortype = 'torus', L = 1)[[1]]
  saveRDS(finemap_res, outfile)
}

# finemap_sumstats <- merge_susie_sumstats(susie_finemap_L1, sumstats.sigloci)
# saveRDS(finemap_sumstats, paste0(data.dir, '/finemapping/AF_finemapping_sumstats_torusprior.rds'))
