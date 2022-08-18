
## Functions for making track plots
finemapping_annot_trackplot <- function(finemapstats,
                                        region,
                                        gene.annots,
                                        bigSNP,
                                        genome = "hg19",
                                        countsdata,
                                        peaks,
                                        HiC_loops,
                                        txdb=NULL,
                                        data_colors = seq_along(countsdata),
                                        highlight_snps = NULL,
                                        highlight_colors = "pink") {

  cat("Making trackplots ...\n")

  # Prepare GWAS summary stats
  if( min(finemapstats$pval) >=0 && max(finemapstats$pval) <= 1 ){
    cat("Convert GWAS p-value to -log10(pvalue). \n")
    finemapstats$pval <- -log10(finemapstats$pval)
  }
  seqlevelsStyle(finemapstats) <- "UCSC"

  # Limit to genomic region to visualize
  seqlevelsStyle(region) <- "UCSC"
  curr_finemapstats <- finemapstats[(as.character(seqnames(finemapstats)) == as.character(seqnames(region))) &
                                      (finemapstats$pos >= start(region)) &
                                      (finemapstats$pos <= end(region)),  ]
  curr_finemapstats <- as.data.frame(curr_finemapstats)

  # p-value track
  # Add LD information
  if(!missing(bigSNP)){
    curr_finemapstats <- add_LD_bigSNP(curr_finemapstats, bigSNP)
    cat(nrow(curr_finemapstats), "snps included.\n")

    pval.df <- curr_finemapstats %>% dplyr::select(chr, pos, pval, r2) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pval.df <- pval.df %>% tidyr::pivot_wider(names_from = r2, names_sort = TRUE, values_from = "pval")
    pval.gr <- makeGRangesFromDataFrame(pval.df, keep.extra.columns = T)
    seqlevelsStyle(pval.gr) <- "UCSC"
    avail_ld_groups <- names(mcols(pval.gr))

    ld_colors <- c("black","blue","green","orange","red")
    names(ld_colors) <- c("0-0.1","0.1-0.25","0.25-0.75","0.75-0.9","0.9-1")
    pval.track <- DataTrack(range = pval.gr,
                            genome = genome,
                            groups = avail_ld_groups,
                            col = ld_colors[avail_ld_groups],
                            name = "-log10 P",
                            legend = TRUE,
                            box.legend = TRUE)
  }else{
    cat(nrow(curr_finemapstats), "snps included.\n")

    pval.df <- curr_finemapstats %>% dplyr::select(chr, pos, pval) %>%
      dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
    pval.gr <- makeGRangesFromDataFrame(pval.df, keep.extra.columns = T)
    seqlevelsStyle(pval.gr) <- "UCSC"

    pval.track <- DataTrack(range = pval.gr,
                            genome = genome,
                            name = "-log10 P",
                            col = "black",
                            legend = FALSE)
  }

  dpars.pval <- list(col.title = "black",
                     col.axis = "black",
                     col.border.title = "lightgray",
                     col.frame = "lightgray",
                     frame = TRUE,
                     cex.axis = 0.6)
  displayPars(pval.track) <- dpars.pval

  # PIP track
  pip.df <- curr_finemapstats %>% dplyr::select(chr, pos, pip, locus) %>%
    dplyr::mutate(start = pos, end = pos) %>% dplyr::select(-pos)
  pip.df <- pip.df %>% tidyr::pivot_wider(names_from = locus, names_sort = TRUE, values_from = "pip")
  pip.gr <- makeGRangesFromDataFrame(pip.df, keep.extra.columns = T)
  seqlevelsStyle(pip.gr) <- "UCSC"
  avail_locus_groups <- names(mcols(pip.gr))

  pip.track <- DataTrack(range = pip.gr,
                         genome = genome,
                         groups = avail_locus_groups,
                         name = "PIP",
                         legend = FALSE)

  dpars.pip <- list(col.title = "black",
                    col.axis = "black",
                    col.border.title = "lightgray",
                    col.frame = "lightgray",
                    frame = TRUE,
                    cex.axis = 0.6)
  displayPars(pip.track) <- dpars.pip

  # Data tracks
  if(!missing(countsdata) && (length(countsdata) > 0)){
    dpars.data <- list(col.title = "black",
                       col.axis = "black",
                       col.border.title = "lightgray",
                       col.frame = "lightgray",
                       # rotation.title = 0,
                       frame = FALSE,
                       cex.axis = 0.2)

    countsdata.tracks <- lapply(names(countsdata), function(x){
      cat("Adding", x, "track...\n")
      countsdata.gr <- countsdata[[x]]
      i <- which(names(countsdata) == x)
      seqlevelsStyle(countsdata.gr) <- "UCSC"
      seqlevels(countsdata.gr, pruning.mode = "coarse") <- paste0("chr",1:22)
      countsdata.track <- DataTrack(range = countsdata.gr,
                                    type = 'h',
                                    genome = genome,
                                    col = data_colors[i],
                                    name = x,
                                    showAxis=FALSE,
                                    ylim = c(0,0.8))

      displayPars(countsdata.track) <- dpars.data
      countsdata.track
    })
  }else{
    countsdata.tracks <- NULL
  }

  # Peak annotation tracks
  if(!missing(peaks) && (length(peaks) > 0)){
    dpars.peaks <- list(col.title = "black",
                        col.axis = "black",
                        col.border.title = "lightgray",
                        col.frame = "lightgray",
                        # rotation.title = 0,
                        frame = FALSE,
                        cex.axis = 0.2)

    peaks.tracks <- lapply(names(peaks), function(x){
      cat("Adding", x, "track... \n")
      peaks.gr <- peaks[[x]]
      seqlevelsStyle(peaks.gr) <- "UCSC"
      seqlevels(peaks.gr, pruning.mode = "coarse") <- paste0("chr",1:22)
      peaks.gr$score <- 1
      peaks.gr <- peaks.gr[countOverlaps(peaks.gr, region) > 0]
      peaks.track <- DataTrack(range = peaks.gr,
                               type = "h",
                               genome = genome,
                               name = x,
                               col = "navy", showAxis=F,
                               ylim = c(0,1))
      displayPars(peaks.track) <- dpars.peaks
      peaks.track
    })
  }else{
    peaks.tracks <- NULL
  }

  # HiC loops tracks
  if(!missing(HiC_loops) && (length(HiC_loops) > 0)){
    dpars.HiC <- list(col.interactions = "black",
                      col.anchors.fill = "blue",
                      col.anchors.line = "black",
                      interaction.dimension = "width",
                      interaction.measure = "score",
                      plot.trans = FALSE,
                      plot.outside = FALSE,
                      col.outside="lightblue",
                      anchor.height = 0.1,
                      # rotation.title = 0,
                      col.title = "black",
                      col.axis = "black",
                      col.border.title = "lightgray",
                      col.frame = "lightgray",
                      frame = FALSE,
                      cex.axis = 0.2)

    gene.annots$chr <- as.character(seqnames(gene.annots))
    gene.annots$tss <- start(resize(gene.annots, width = 1))

    HiC_loops.tracks <- lapply(names(HiC_loops), function(x){
      cat("Adding", x, "track...\n")
      HiC_loops.gr <- HiC_loops[[x]]
      HiC_loops.gr <- HiC_loops.gr[HiC_loops.gr$gene_name %in% gene.annots$gene_name,]
      m <- match(HiC_loops.gr$gene_name, gene.annots$gene_name)
      HiC_loops.gr$promoter_chr <- gene.annots$chr[m]
      HiC_loops.gr$promoter_tss <- gene.annots$tss[m]

      HiC_promoters.gr <- GRanges(seqnames = HiC_loops.gr$promoter_chr,
                                  ranges = IRanges(start = HiC_loops.gr$promoter_tss, end = HiC_loops.gr$promoter_tss),
                                  score = HiC_loops.gr$score,
                                  gene = HiC_loops.gr$gene_name)
      HiC_enhancers.gr <- GRanges(seqnames = seqnames(HiC_loops.gr),
                                  ranges = IRanges(start = start(HiC_loops.gr), end = end(HiC_loops.gr)),
                                  score = HiC_loops.gr$score)
      HiC_loops.obj <- GenomicInteractions::GenomicInteractions(anchor1 = HiC_promoters.gr, anchor2 = HiC_enhancers.gr)
      HiC_loops.obj$score <- HiC_loops.obj$anchor1.score
      HiC_loops.track <- InteractionTrack(HiC_loops.obj, name = x)
      displayPars(HiC_loops.track) <- dpars.HiC
      HiC_loops.track
    })
  }else{
    HiC_loops.tracks <- NULL
  }

  # gene track
  gene.track <- knownGeneObject(region, genome = genome, track_name = "Genes")
  dpars.genes <- list(col.title = "black",
                      col.axis = "black",
                      col.border.title = "lightgray",
                      col.frame = "lightgray",
                      # rotation.title = 0,
                      frame = FALSE,
                      cex.axis = 0.2)
  displayPars(gene.track) <- dpars.genes

  # # restrict to protein coding genes
  # gene.track <- gene.track[symbol(gene.track) %in% gene.annots$gene_name]

  # axis track
  axisTrack <- GenomeAxisTrack()

  # List all tracks
  list.of.tracks <- c(pval.track,
                      pip.track,
                      countsdata.tracks,
                      peaks.tracks,
                      HiC_loops.tracks,
                      gene.track,
                      axisTrack)
  track.sizes <- c(1,
                   0.6,
                   rep(0.4, length(countsdata.tracks)),
                   rep(0.3, length(peaks.tracks)),
                   rep(0.6, length(HiC_loops.tracks)),
                   0.5,
                   0.4)

  # Highlight SNPs
  if( !missing(highlight_snps) && (length(highlight_snps) > 0)){

    if("topSNP" %in% highlight_snps){
      topsnp_idx <- which.max(curr_finemapstats$pip)
      highlight_snps[which(highlight_snps == "topSNP")] <- curr_finemapstats$snp[topsnp_idx]
    }
    highlight_snps <- intersect(highlight_snps, curr_finemapstats$snp)
    highlight_pos <- curr_finemapstats$pos[match(highlight_snps, curr_finemapstats$snp)]
    cat("Highlight SNPs:", highlight_snps[which(highlight_snps %in% curr_finemapstats$snp)], "\n")
    if(length(highlight_pos) == 1){
      list.of.tracks <- HighlightTrack(trackList = list.of.tracks,
                                       start = c(highlight_pos-500), width = 1000,
                                       chromosome = as.character(seqnames(region)),
                                       col = highlight_colors)
    }else if(length(highlight_pos) > 1){
      cat("Highlight SNP positions:", highlight_pos, "\n")
      if(length(highlight_colors) > 1){
        highlight_colors <- highlight_colors[which(highlight_snps %in% curr_finemapstats$snp)]
      }
      list.of.tracks <- HighlightTrack(trackList = list.of.tracks,
                                       start = highlight_pos, width = 1,
                                       chromosome = as.character(seqnames(region)),
                                       col = highlight_colors)
    }
  }

  plotTracks(list.of.tracks,
             chromosome = as.character(seqnames(region)),
             transcriptAnnotation = "symbol",
             collapseTranscripts= "meta",
             from = start(region),
             to = end(region),
             sizes = track.sizes,
             just.group = "right")

}


knownGeneObject <- function(region, genome, txdb, keytype = "ENTREZID", track_name = "Genes"){

  if(missing(txdb)){
    cat("Making gene track...\n")
    if(genome == "hg19"){
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    }else if(genome == "hg38"){
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }else{
      stop("Genome needs to be hg19 or hg38, or please specify txdb. ")
    }
  }else{
    cat("Load txdb...\n")
  }

  grtrack <- Gviz::GeneRegionTrack(range = txdb,
                                   genome = genome,
                                   chromosome = as.character(seqnames(region)),
                                   start = start(region),
                                   end = end(region),
                                   name = track_name)

  symbol(grtrack) <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                           keys=sub("\\.\\d+$", "", gene(grtrack)),
                                           keytype=keytype, column="SYMBOL")

  symbol(grtrack) <- ifelse(is.na(symbol(grtrack)), gene(grtrack), symbol(grtrack))

  return(grtrack)
}

# Add LD information
add_LD_bigSNP <- function(sumstats, bigSNP,
                          r2_breaks = c(0, 0.1, 0.25, 0.75, 0.9, 1),
                          r2_labels = c("0-0.1","0.1-0.25","0.25-0.75","0.75-0.9","0.9-1")) {

  # only include SNPs in bigSNP markers
  sumstats <- sumstats[sumstats$snp %in% bigSNP$map$marker.ID, ]
  sumstats$bigSNP_idx <- match(sumstats$snp, bigSNP$map$marker.ID)

  locus_list <- unique(sumstats$locus)

  sumstats.r2.df <- data.frame()
  for(locus in locus_list){
    curr_sumstats <- sumstats[sumstats$locus == locus, ]
    top_snp_idx <- curr_sumstats$bigSNP_idx[which.max(curr_sumstats$pval)]
    top_snp_genotype <- bigSNP$genotypes[,top_snp_idx]
    genotype.mat <- bigSNP$genotypes[,curr_sumstats$bigSNP_idx]

    r2.vals <- as.vector(cor(top_snp_genotype, genotype.mat))^2
    r2.brackets <- cut(r2.vals, breaks = r2_breaks, labels = r2_labels)
    curr_sumstats$r2 <- r2.brackets
    sumstats.r2.df <- rbind(sumstats.r2.df, curr_sumstats)
  }

  return(sumstats.r2.df)
}

# get gene region
get_gene_region <- function(gene.mapping.res, genes.of.interest, ext = 10000, select.region = c("all", "locus")){

  select.region <- match.arg(select.region)

  high.conf.snp.df <- gene.mapping.res %>% dplyr::filter(pip > 0.2) %>%
    dplyr::group_by(snp) %>% dplyr::arrange(-gene_pip) %>% dplyr::slice(1)
  gene.gr <- gene.annots[match(high.conf.snp.df$gene_name, gene.annots$gene_name),]
  gene.gr$tss <- start(resize(gene.gr, width = 1))
  gene.gr <- gene.gr[,c("gene_name","tss")]
  high.conf.snp.df$tss <- gene.gr$tss

  gene.snp.tss <- high.conf.snp.df %>%
    dplyr::filter(gene_name %in% genes.of.interest) %>%
    dplyr::group_by(locus) %>%
    dplyr::arrange(-pip) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(distToTSS = pos-tss) %>%
    dplyr::select(gene_name, locus, chr, pos, tss, distToTSS)

  if(select.region == "all"){
    chr <- gene.snp.tss$chr[1]
    locus <- paste(gene.snp.tss$locus, collapse = ",")
    region_start <- min(c(gene.snp.tss$pos, gene.snp.tss$tss)) - ext
    region_end <- max(c(gene.snp.tss$pos, gene.snp.tss$tss)) + ext
    region <- GRanges(seqnames = chr,
                      IRanges(start = region_start, end = region_end),
                      locus = locus)
    seqlevelsStyle(region) <- "UCSC"
  }else if(select.region == "locus"){
    region <- lapply(gene.snp.tss$locus, function(l){
      locus <- l
      locus.snp.tss <- gene.snp.tss[gene.snp.tss$locus == locus, ]
      chr <- locus.snp.tss$chr
      if(distToTSS < 0){ #snp is upstream
        region_start <- locus.snp.tss$pos - ext
        region_end <- locus.snp.tss$tss + ext
      } else{ # snp is downstream
        region_start <- locus.snp.tss$tss - ext
        region_end <- locus.snp.tss$pos + ext
      }
      gene_locus_region.gr <- GRanges(seqnames = chr,
                                      IRanges(start = region_start, end = region_end),
                                      locus = locus)
      seqlevelsStyle(gene_locus_region.gr) <- "UCSC"
      gene_locus_region.gr
    })
  }
  return(region)
}
