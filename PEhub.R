# ==================================================================================================
# PEhub v1.0 — CLI I/O Contract (comment-only; implement argument parsing in a thin wrapper script)
#
# Goal
#   Run the complete pipeline non-interactively via a single Rscript command where ALL inputs/outputs are supplied
#   from the command line (no hard-coded file paths in the wrapper). This file contains the core functions and
#   example invocations; production use should call the wrappers with CLI-provided arguments.
#
# Required CLI inputs (typical)
#   --gtf                Path to GTF (or an equivalent annotation source used to construct promoter GRanges).
#   --loop_sig           Path to significant interactions BEDPE (HicDCPlus significant_interactions).
#   --loop_all           Path to all interactions BEDPE(.gz) for null pool construction.
#   --outdir             Output directory for all intermediate and final deliverables.
#   --sample             Sample identifier used in output file naming (e.g., 'Caudate1').
#
# Optional CLI inputs (recommended to expose)
#   --mode               One of: test | brain | custom (controls which sample set is iterated).
#   --promoter_window    Promoter window around TSS (bp).
#   --method             Co-membership normalization: log_minmax | log_zscore | log_maxnorm.
#   --weight_mode        Weight model (see compute_weights() for enumerated modes).
#   --k_min              Minimum enhancers per promoter to consider hub discovery.
#   --resolution         Leiden resolution parameter.
#   --quantile_cutoff    Edge sparsification quantile for co-membership matrices.
#   --B_stability        Bootstrap iterations for stability (calculate_hub_stability).
#   --B_null             Monte Carlo iterations for global p-values (build.null.pvalue.calculate).
#   --null_mode          Null sampling mode: distance_matched | hist_matched | global.
#   --stat               Hub score statistic for null: density | sum.
#   --workers            Parallel workers (future::plan).
#   --python             Python interpreter path for leiden via reticulate.
#
# Outputs (all written under --outdir, names parameterized by --sample/--method/--weight_mode)
#   1) multiple_result.<sample>.hub.all.preprocess.RData
#   2) multiple_result.<sample>.hub.<method>.<weight_mode>.preprocess.RData
#   3) multiple_result.<sample>.hub.<method>.<weight_mode>.RData
#   4) multiple_result.<sample>.hub.<method>.<weight_mode>.bed / .txt / .bedpe
#   5) multiple_result.<sample>.pairwise.<method>.<weight_mode>.bed / .txt / .bedpe
#
# Notes
#   - This script does NOT implement argument parsing; your CLI wrapper should parse args (optparse/argparse) and
#     then call build.run.preprocess(), build.run.all(), build.cutoff(), and build.postprocess() with those args.
#   - The example invocations at the bottom are illustrative; for true CLI execution, drive them via args.
# ==================================================================================================


# ==============================================================================
# [0] Dependencies and runtime configuration
#   - Install R/Bioconductor packages once if missing.
#   - Configure parallel backend for future/furrr.
# ==============================================================================
# BiocManager::install(c("GenomicRanges", "dplyr", "purrr", "igraph", "parallel"))
library(optparse)
library(GenomicRanges)
library(dplyr)
library(purrr)
library(igraph)
library(parallel)
library(Matrix)
library(tidyr)
library(stringr)
library(furrr)
library(data.table)
library(progress) 

# ==============================================================================
# [0.1] Inputs (planned to be provided via command-line args)
# ==============================================================================
option_list <- list(
  make_option("--tss_gtf", type="character", default=NULL,
              help="GTF file for gene annotation (transcript features used for TSS/promoters)."),
  make_option("--loop_sig", type="character", default=NULL,
              help="Significant interactions BEDPE (used for observed hubs)."),
  make_option("--loop_all", type="character", default=NULL,
              help="All interactions BEDPE(.gz) (used for global null pool)."),
  make_option("--outdir", type="character", default=NULL,
              help="Output directory."),
  make_option("--sample", type="character", default=NULL,
              help="Output filename prefix / sample name."),
  make_option("--python_path", type="character", default=NULL,
              help="Python path."),

  make_option("--method", type="character", default="log_minmax",
              help="Co-membership normalization method: log_minmax | log_zscore | log_maxnorm"),
  make_option("--weight_method", type="character", default="bin_log_ratio_sig",
              help="Edge weight mode in compute_weights()."),

  make_option("--promoter_window", type="integer", default=0,
              help="Promoter window around TSS (bp)."),
  make_option("--k_min", type="integer", default=3,
              help="Minimum enhancers per promoter to consider hub."),
  make_option("--resolution", type="double", default=1.0,
              help="Leiden resolution parameter."),
  make_option("--quantile_cutoff", type="double", default=0.2,
              help="Quantile cutoff for sparsifying comembership matrix."),
  make_option("--pvalue_cutoff", type="double", default=0.05,
              help="p-value cutoff for each hub."),
  make_option("--stability_cutoff", type="double", default=0.5,
              help="Hub stability cutoff for each hub."),

  make_option("--B_stability", type="integer", default=10,
              help="Bootstrap iterations for stability."),
  make_option("--B_pvalue", type="integer", default=1000,
              help="Null iterations for p-value."),
  make_option("--null_mode", type="character", default="hist_matched",
              help="Null mode: distance_matched | hist_matched | global"),

  make_option("--workers", type="integer", default=10,
              help="Number of parallel workers for future/furrr.")
)

opt <- parse_args(OptionParser(option_list = option_list))
required <- c("tss_gtf","loop_sig","loop_all","outdir","sample","python_path")
tss_gtf <- opt$tss_gtf
loop_file <- opt$loop_sig
loop_file_all <- opt$loop_all
outdir <- opt$outdir
filename <- opt$sample
python_path <- opt$python_path

method <- opt$method
weight_method <- opt$weight_method
promoter_window <- opt$promoter_window
k_min <- opt$k_min
resolution <- opt$resolution
quantile_cutoff <- opt$quantile_cutoff
pvalue_cutoff <- opt$pvalue_cutoff
stability_cutoff <- opt$stability_cutoff
B_stability <- opt$B_stability
B_pvalue <- opt$B_pvalue
null_mode <- opt$null_mode
workers <- opt$workers
# Parallel workers (planned: args$workers). Tune based on cluster resources.
plan(multisession, workers = workers)

# ==============================================================================
# [0.2] Leiden clustering backend
#   - Uses the 'leiden' R package (and optionally python via reticulate).
#   - Python path should be configurable via args$python_bin (planned).
# ==============================================================================
library(reticulate)
use_python(python_path, required = TRUE)
py_config() 
library(leiden)


# ==============================================================================
# [1] Promoter/TSS annotation
#
# Goal
#   Build promoter/TSS annotations for mapping loop anchors to promoters.
#
# Supported annotation modes (documentation only)
#   (A) TxDb: derive promoters from TxDb.Hsapiens.UCSC.hg38.knownGene
#   (B) OrgDb: derive gene coordinates from org.Hs.eg.db
#   (C) GTF: parse transcript features and use transcript start sites (TSS)
#
# Current implementation
#   Uses GTF transcript records. This will be parameterized by args$tss_gtf.
# ==============================================================================
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# promoters_annotation <- promoters(
#   TxDb.Hsapiens.UCSC.hg38.knownGene,
#   upstream = 0,
#   columns  = "gene_id"
# )


# ## orgdb 模式：从 org.Hs.eg.db 的基因坐标取 promoter
# library(org.Hs.eg.db)
# genes_gr <- genes(org.Hs.eg.db)
# promoters_annotation <- promoters(
#   genes_gr,
#   upstream = 0,
#   columns  = "gene_id"
# )


#' @title Extract Transcript TSS from GTF File
#' @description This function processes a GTF file to extract the transcription start sites (TSS) 
#' of transcripts. It is designed to work with genomic data in GTF format.
#' @details The function reads the GTF file, identifies transcript features, and extracts their 
#' transcription start sites (TSS). The extracted TSS data is then saved to the specified output file 
#' in a format suitable for downstream analysis.
library(rtracklayer)
gtf <- rtracklayer::import(tss_gtf)
promoters_annotation <- gtf %>%
  subset(type == "transcript") %>%
  unique()



# ==============================================================================
# [2] Preprocess HiChIP loops (BEDPE) and annotate anchors
#   - Standardize loop anchors and label P/E types.
#   - Map anchors to promoter/enhancer peak IDs.
#   - Optionally annotate H3K27ac peak overlap.
# Outputs:
#   list(loops, hubs, peaks)
# ==============================================================================

# [F01] loops_convert()
# Purpose : Normalize and augment a HiChIP/Hi-C loop table into a unified 'hubs' table representation.
# Inputs  : dt (data.table/data.frame) containing a loop BEDPE-like schema plus PEhub metadata columns.
# Outputs : data.frame with standardized columns (chr/start/end/...) and one synthetic PP anchor row per hubid.
# Notes   : This function is pure in-memory transformation; all file I/O must be handled by the CLI wrapper.
loops_convert <- function(dt) {

  setDT(dt)

  rest <- dt[, -(1:3)]
  data.table::setnames(rest, old = names(rest)[1:3], new = c("chr", "start", "end"))

  anchor_cols <- c(
    "chr1","start1","end1","type1","type2","pair_type",
    "promoter_id","promoter_gene_id","promoter_transcript_id",
    "promoter_gene_name","hubid"
  )
  anchor_cols <- intersect(anchor_cols, names(dt))

  anchor <- dt[, ..anchor_cols]
  data.table::setnames(anchor, old = c("chr1","start1","end1"), new = c("chr","start","end"), skip_absent = TRUE)

  anchor[, pair_type := "PP"]

  if ("hubid" %in% names(anchor)) {
    data.table::setkey(anchor, hubid)
    anchor <- anchor[!duplicated(hubid)]
  }

  miss <- setdiff(names(rest), names(anchor))
  if (length(miss)) {
    anchor[, (miss) := NA]
  }

  data.table::setcolorder(anchor, names(rest))

  final <- data.table::rbindlist(list(rest, anchor), use.names = TRUE, fill = TRUE)

  if (all(c("enhancer_id", "promoter_id") %in% names(final))) {
    final[is.na(enhancer_id), enhancer_id := promoter_id]
  }

  if ("hubid" %in% names(final)) {
    data.table::setorder(final, hubid)
  }

  return(as.data.frame(final))
}


# [F02] make_peak_id()
# Purpose : Generate a stable peak identifier from genomic coordinates using a fixed binning scale.
# Inputs  : chr, start, end; optional scale (default 1000L).
# Outputs : character vector of peak IDs (e.g., peak_chr1_123_124).
# Notes   : Coordinate binning is an implementation detail; keep consistent across runs for reproducibility.
make_peak_id <- function(chr, start, end, scale = 1000L) {
  paste0("peak_",chr, "_", start %/% scale, "_", end %/% scale)
}

# [F03] preprocess_hichip()
# Purpose : Ingest a BEDPE(-like) interaction file, annotate anchors as Promoter/Enhancer, and build PEhub-ready tables.
# Inputs  : loop_file (path via CLI), tss_input (GRanges with gene_id/gene_name/transcript_id), mode ('all' or 'promoter'),
#          promoter_window, optional h3k27ac_peak (path via CLI), min_dist, max_dist.
# Outputs : list(loops=annotated loops table, hubs=loops converted to hub-style table, peaks=GRanges of unique anchors).
# Notes   : All file paths must be passed as CLI args and bound to the corresponding variables before calling.
preprocess_hichip <- function(
    loop_file,
    tss_input = NULL,                  # GRanges object，include gene_id column
    mode = c("all"),       # all = EP+PP, promoter = 只EP
    promoter_window = 0,
    h3k27ac_peak = NULL,
    min_dist = 10000,
    max_dist = 2000000
) {
  message("Running in ", mode, " mode")

  # load promoter annotation
  if (!"gene_id" %in% names(mcols(tss_input))) {
    stop("Your TSS object must have 'gene_id' metadata column")
  }
  promoters_genome <- promoters(tss_input, upstream = promoter_window, downstream = promoter_window)
  message("Got ", length(tss_input), " promoters from annotation")
  
  # read loops
  loops_ori <- data.table::fread(loop_file, header=TRUE)
  colnames(loops_ori) = c("chr1", "start1", "end1", "chr2", "start2", "end2",colnames(loops_ori)[7:ncol(loops_ori)])
  #loops_ori = loops_ori[loops_ori$counts>0, ]
  message("Initial loop number: ", nrow(loops_ori))
  
  # take each side of the loop as peaks and annotate gene_id / transcript_id
  gr1 <- GRanges(loops_ori$chr1, IRanges(loops_ori$start1, loops_ori$end1))    # bed是0-based
  gr2 <- GRanges(loops_ori$chr2, IRanges(loops_ori$start2, loops_ori$end2))
  peaks <- sort(unique(c(gr1, gr2)))
  peaks$name <- make_peak_id(as.character(seqnames(peaks)), start(peaks), end(peaks))

  get_overlap_gene_id <- function(query_gr, promoter_gr, id_col="gene_id") {
    hits <- findOverlaps(query_gr, promoter_gr)
    # build a data.frame mapping query index → gene_id
    out <- data.frame(
      query_index = queryHits(hits),
      gene_id = mcols(promoter_gr)[,id_col][subjectHits(hits)],
      stringsAsFactors = FALSE
    )
    # If multiple promoters overlap a single query region → collapse to comma-separated
    out2 <- out %>%
      group_by(query_index) %>%
      summarise(gene_id = paste(unique(gene_id), collapse = ",")) %>%
      right_join(
        data.frame(query_index = seq_along(query_gr)),
        by = "query_index"
      ) %>% arrange(query_index)
    return(out2$gene_id)
  }

  gene_anno = get_overlap_gene_id(peaks, promoters_genome,"gene_id")
  genename_anno = get_overlap_gene_id(peaks, promoters_genome,"gene_name")
  transcript_anno = get_overlap_gene_id(peaks, promoters_genome,"transcript_id")
  peaks$gene_id = gene_anno
  peaks$transcript_id = transcript_anno
  peaks$gene_name = genename_anno
  promoters = peaks[unique(queryHits(findOverlaps(peaks, promoters_genome)))]
  enhancers = peaks[-unique(queryHits(findOverlaps(peaks, promoters_genome)))]

  # 1. label P / E 
  hits1 <- findOverlaps(gr1, promoters, type = "equal")
  type1 <- ifelse(seq_along(gr1) %in% unique(queryHits(hits1)), "P", "E")
  hits2 <- findOverlaps(gr2, promoters, type = "equal")
  type2 <- ifelse(seq_along(gr2) %in% unique(queryHits(hits2)), "P", "E")
  loops <- loops_ori %>% mutate(type1 = type1, type2 = type2,
           pair_type = paste0(type1, type2))
  message("Initial loop type number: ", nrow(loops), 
          "  PE = ", sum(loops$pair_type=="PE"), 
          "  EP = ", sum(loops$pair_type=="EP"), 
          ", EE = ", sum(loops$pair_type=="EE"),
          ", PP = ", sum(loops$pair_type=="PP"), ")")
  
  # 2. Swap EP to PE（promoter in the front）
  setDT(loops)

  idx <- loops[type1 == "E" & type2 == "P", which = TRUE]

  if (length(idx) > 0) {
    # swap coords in-place
    tmp <- loops[idx, .(chr1, start1, end1, type1)]
    loops[idx, `:=`(
      chr1   = chr2,
      start1 = start2,
      end1   = end2,
      chr2   = tmp$chr1,
      start2 = tmp$start1,
      end2   = tmp$end1,
      type1  = type2,
      type2  = tmp$type1
    )]
  }

  loops[, pair_type := paste0(type1, type2)]
  table(loops$pair_type)
  message("Initial loop type swap number: ", nrow(loops), 
          "  PE = ", sum(loops$pair_type=="PE"), 
          ", EE = ", sum(loops$pair_type=="EE"),
          ", PP = ", sum(loops$pair_type=="PP"), ")")

  # 3. get id for ehancer / promoter
  get_id <- function(gr, range_set) {
    hits <- findOverlaps(gr, range_set, type = "equal")
    return(list(name=range_set[subjectHits(hits)]$name,
                gene_id=range_set[subjectHits(hits)]$gene_id,
                transcript_id=range_set[subjectHits(hits)]$transcript_id,
                gene_name=range_set[subjectHits(hits)]$gene_name
                ))
  }
  
  node1_id <- get_id(GRanges(loops$chr1, IRanges(loops$start1, loops$end1)), peaks)
  node2_id <- get_id(GRanges(loops$chr2, IRanges(loops$start2, loops$end2)), peaks)  
  loops <- loops %>% mutate(promoter_id = node1_id$name, enhancer_id=node2_id$name
                         , promoter_gene_id = node1_id$gene_id, enhancer_gene_id = node2_id$gene_id
                         , promoter_transcript_id = node1_id$transcript_id, enhancer_transcript_id = node2_id$transcript_id
                          , promoter_gene_name = node1_id$gene_name
                         )
  message("Before filtering loop number: ", nrow(loops))

  # keep only PE or not
  if (mode == "promoter") {
    loops <- loops %>% filter(pair_type %in% c("PE"))
    message("Only keep enhancer promoter interactions, loop number: ", nrow(loops))
  } else {
    loops <- loops %>% filter(pair_type %in% c("PE","EE"))
    message("Keep all PE and EE interaction types, loop number: ", nrow(loops))
  }
    
  # Filtering using distance
  loops <- loops %>%
    mutate(distance = abs((start1+end1)/2 - (start2+end2)/2)) %>%
    filter(distance >= min_dist, distance <= max_dist)
  message("After distance filtering, loop number: ", nrow(loops))
  
  # 6. if H3K27ac peaks overlap 
  if(!is.null(h3k27ac_peak)){
    h3k27ac_peaks <- rtracklayer::import(h3k27ac_peak)
    
    gr_left  <- GRanges(loops$chr1, IRanges(loops$start1 + 1, loops$end1))
    gr_right <- GRanges(loops$chr2, IRanges(loops$start2 + 1, loops$end2))

    over_left  <- overlapsAny(gr_left,  h3k27ac_peaks + 5000)
    over_right <- overlapsAny(gr_right, h3k27ac_peaks + 5000)

    # loops <- loops[over_left & over_right, ]
    loops$peakoverlap1 <- over_left
    loops$peakoverlap2 <- over_right
    message("After H3K27ac peaks filtering, loop number: ", nrow(loops))
  }

  setDT(loops)
  loops[, hubid := paste0("hub_", frankv(list(chr1, start1, end1), ties.method = "dense"))]
  
  message("Final loop total number: ", nrow(loops), 
          "  PE = ", sum(loops$pair_type=="PE"), 
          ", EE = ", sum(loops$pair_type=="EE"),
          ", PP = ", sum(loops$pair_type=="PP"), ")")
  
  hubs <- loops_convert(loops)

  return(list(loops=as.data.frame(loops), hubs=hubs, peaks=peaks))

}

# [F04] compute_weights()
# Purpose : Compute per-EP interaction weights under multiple scoring models (counts, significance, distance, bin-corrected, etc.).
# Inputs  : df (loops table), weight_mode, dist_col, q_col, breaks, alpha, sig_cap, sig_beta.
# Outputs : list(loops=loops with weight_raw/weight_percentage, hubs=loops converted to hub-style table).
# Notes   : Weight definitions must be reported in Methods; record weight_mode and parameters in output metadata externally.
compute_weights <- function(df,
                            weight_mode = c("count_only",
                                            "sig_only",
                                            "distance_only",
                                            "count_sig",
                                            "zscore_residual",
                                            "count_sig_plus_dist_linear",
                                            "bin_percentile_plus_sig",
                                            "bin_log_ratio_sig",
                                            "bin_diff_global",
                                            "bin_diff_binmax",
                                            "bin_percentile"),
                            dist_col = c("D"),
                            q_col = "qvalue",
                            breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6),
                            alpha = 0.8,
                            sig_cap = 10,
                            sig_beta = 0.5) {

  weight_mode <- match.arg(weight_mode)
  dist_col <- dist_col[dist_col %in% colnames(df)][1]
  if (is.na(dist_col)) stop("No distance column found. Provide 'distance' or 'D' in df, or set dist_col.")

  # --- base quantities (always computed) ---
  df <- df %>%
    mutate(
      w_count = log1p(counts),
      w_sig_raw = if (q_col %in% colnames(.)) -log10(.data[[q_col]] + 1e-10) else NA_real_,
      w_sig_cap = ifelse(is.na(w_sig_raw), NA_real_, pmin(w_sig_raw, sig_cap)),
      # two variants:
      w_sig_mult = ifelse(is.na(w_sig_cap), 1, w_sig_cap),                              # multiplicative version (stronger)
      reliability = ifelse(is.na(w_sig_cap), 1, 1 + sig_beta * (w_sig_cap / sig_cap)),  # 1 ~ 1+beta (recommended)
      dist = .data[[dist_col]]
    )

  # helper: global normalize to [0,1]
  norm01 <- function(x) {
    mx <- max(x, na.rm = TRUE)
    if (!is.finite(mx) || mx <= 0) return(rep(0, length(x)))
    x / mx
  }

  # =========================
  # choose weight_mode
  # =========================
  if (weight_mode == "count_only") {
    # baseline: only counts
    df <- df %>%
      mutate(
        weight_raw = w_count,
        weight_percentage = norm01(weight_raw)
      )

  } else if (weight_mode == "sig_only") {
    # baseline: only significance (capped), treat NA as 0 signal
    df <- df %>%
      mutate(
        weight_raw = ifelse(is.na(w_sig_cap), 0, w_sig_cap),
        weight_percentage = norm01(weight_raw)
      )

  } else if (weight_mode == "distance_only") {
    # negative control: only distance (nearer gets larger weight)
    # use inverse log-distance; add 1 to avoid div-by-zero when dist=0
    df <- df %>%
      mutate(
        weight_raw = 1 / log1p(pmax(dist, 0)),
        weight_percentage = norm01(weight_raw)
      )

  } else if (weight_mode == "count_sig") {
    # weight_raw = log1p(counts) * capped(-log10(q))
    df <- df %>%
      mutate(
        weight_raw = w_count * w_sig_mult,
        weight_percentage = norm01(weight_raw)
      )

  } else if (weight_mode == "count_sig_plus_dist_linear") {
    # your "alpha*signal + (1-alpha)*dist_norm" (global dist norm)
    df <- df %>%
      mutate(
        w_dist_norm = {
          x <- log1p(dist)
          m <- max(x, na.rm = TRUE)
          if (!is.finite(m) || m == 0) rep(0, length(x)) else x / m
        },
        weight_raw_signal = w_count * w_sig_mult,
        weight_raw = alpha * norm01(weight_raw_signal) + (1 - alpha) * w_dist_norm,
        weight_percentage = norm01(weight_raw)
      )

  } else if (weight_mode == "bin_diff_global") {
    # w_dist = pmax(w_count - bin_median, 0); then global normalize
    df <- df %>%
      mutate(
        bin = cut(dist, breaks = breaks, labels = FALSE, include.lowest = TRUE)
      ) %>%
      group_by(bin) %>%
      mutate(
        bin_med = median(w_count, na.rm = TRUE),
        w_dist_raw = pmax(w_count - bin_med, 0)
      ) %>%
      ungroup() %>%
      mutate(
        weight_raw = w_dist_raw * reliability,
        weight_percentage = norm01(weight_raw)
      )

  } else if (weight_mode == "bin_diff_binmax") {
    # w_dist_raw = pmax(w_count - bin_median, 0); then BIN-wise normalize; then global normalize
    df <- df %>%
      mutate(
        bin = cut(dist, breaks = breaks, labels = FALSE, include.lowest = TRUE)
      ) %>%
      group_by(bin) %>%
      mutate(
        bin_med = median(w_count, na.rm = TRUE),
        w_dist_raw = pmax(w_count - bin_med, 0),
        w_dist_bin01 = {
          mx <- max(w_dist_raw, na.rm = TRUE)
          if (!is.finite(mx) || mx <= 0) rep(0, length(w_dist_raw)) else w_dist_raw / mx
        }
      ) %>%
      ungroup() %>%
      mutate(
        weight_raw = w_dist_bin01 * reliability,
        weight_percentage = norm01(weight_raw)
      )

  } else if (weight_mode == "bin_percentile") {
    # weight is purely bin-wise percentile of w_count (no qvalue use)
    df <- df %>%
      mutate(
        bin = cut(dist, breaks = breaks, labels = FALSE, include.lowest = TRUE)
      ) %>%
      group_by(bin) %>%
      mutate(
        w_distcorr = dplyr::percent_rank(w_count)
      ) %>%
      ungroup() %>%
      mutate(
        weight_raw = pmax(w_distcorr, 0),
        weight_percentage = norm01(weight_raw)
      )

  } else if (weight_mode == "bin_percentile_plus_sig") {
    # bin percentile + mild reliability from qvalue
    df <- df %>%
      mutate(
        bin = cut(dist, breaks = breaks, labels = FALSE, include.lowest = TRUE)
      ) %>%
      group_by(bin) %>%
      mutate(
        w_distcorr = dplyr::percent_rank(w_count)
      ) %>%
      ungroup() %>%
      mutate(
        weight_raw = pmax(w_distcorr, 0) * reliability,
        weight_percentage = norm01(weight_raw)
      )

  } else if (weight_mode == "bin_log_ratio_sig") {
    # log(Observed / Expected) within distance bins (Expected = bin median counts)
    df <- df %>%
      mutate(
        bin = cut(dist, breaks = breaks, labels = FALSE, include.lowest = TRUE)
      ) %>%
      group_by(bin) %>%
      mutate(
        bin_med_count = median(pmax(counts, 1), na.rm = TRUE),
        w_dist_ratio = log(pmax(counts, 1) / bin_med_count)
      ) %>%
      ungroup() %>%
      mutate(
        w_dist_raw = pmax(w_dist_ratio, 0),
        weight_raw = w_dist_raw * reliability,
        weight_percentage = norm01(weight_raw)
      )

  } else if (weight_mode == "zscore_residual") {
    # uses HicDCPlus model outputs: mu, sdev (must exist)
    if (!all(c("mu", "sdev") %in% colnames(df))) {
      stop("zscore_residual requires columns 'mu' and 'sdev' in df.")
    }
    df <- df %>%
      mutate(
        z = (counts - mu) / pmax(sdev, 1e-10),
        weight_raw = pmax(z, 0) * reliability,
        weight_percentage = norm01(weight_raw)
      )
  }

  hubs=loops_convert(df)
  return(list(loops=as.data.frame(df), hubs=as.data.frame(hubs)))
}


# [F05] fast_comembership()
# Purpose : Build an enhancer-by-enhancer weighted co-membership matrix for a single promoter using a synergy-like model.
# Inputs  : enhancers (character vector), weights (numeric vector aligned to enhancers), quantile_cutoff, method.
# Outputs : dgCMatrix sparse symmetric matrix with enhancer names as dimnames.
# Notes   : quantile_cutoff sparsifies edges; choose consistently for observed/null computations to avoid bias.
fast_comembership <- function(enhancers, weights, quantile_cutoff = 0, method = c("log_minmax")) { #"log_minmax", "log_zscore", "log_maxnorm"
  # 1. data
  w <- as.numeric(weights)
  w[!is.finite(w)] <- 0
  w <- pmax(w, 0)
  n <- length(w)
  
  if (n < 2) {
    mat <- Matrix(0, n, n, sparse = TRUE)
    dimnames(mat) <- list(enhancers, enhancers)
    return(mat)
  }
  
  # 2. Log-space calculation (preserves your Synergy formula logic perfectly)
  # Original formula: Score_ij ≈ w_i * w_j * [ T_all / ((1+w_i)*(1+w_j)) ]
  # We ignore "-1" since when T_all is large, "-1" has almost no effect on relative ranking
  
  # Precompute log values to prevent overflow
  # Use log1p(x) for more accurate computation of log(1+x)
  log_w <- log(w + 1e-10)    # 防止 log(0)
  log_1_plus_w <- log1p(w)   # log(1+w)
  sum_log_T <- sum(log_1_plus_w) # 这是原来的 log(T_all)
  
  # Co-membership Log(Score)
  # Log(S_ij) = log(w_i) + log(w_j) + log(T_all) - log(1+w_i) - log(1+w_j)
  
  # For better performance, use outer (vector outer product)
  term_wi_wj <- outer(log_w, log_w, "+")         # log(w_i) + log(w_j)
  term_denom <- outer(log_1_plus_w, log_1_plus_w, "+") # log(1+w_i) + log(1+w_j)
  
  log_S_matrix <- term_wi_wj + sum_log_T - term_denom
  
  # 3. Normalization (Key Step)
  # Find the maximum log value in the matrix and shift all values so that the maximum becomes 0 (i.e., original value becomes 1).
  # This ensures that after applying exp(), the values are in the range [0, 1] while maintaining their relative proportions.
  diag(log_S_matrix) <- NA
  max_log_val <- max(log_S_matrix, na.rm = TRUE)
  min_log_val <- min(log_S_matrix, na.rm = TRUE)
  if (method == "log_minmax") {
    den <- (max_log_val - min_log_val)
    if (!is.finite(den) || den <= 0) {
      S_normalized <- matrix(0, n, n)
    } else {
      S_normalized <- (log_S_matrix - min_log_val) / den
    }    
  } else if (method == "log_zscore") {
    mu <- mean(log_S_matrix, na.rm = TRUE)
    sigma <- sd(log_S_matrix, na.rm = TRUE)
    
    if (!is.finite(sigma) || sigma <= 0) {
      S_normalized <- matrix(0, n, n)
    } else {
      S_normalized <- (log_S_matrix - mu) / sigma
      S_normalized <- pmax(S_normalized, 0)
    }
  } else if (method == "log_maxnorm") {
    S_normalized <- exp(log_S_matrix - max_log_val)
  }
  
  # 4. Sparsification - Reduce overly large hubs
  # Keep only the strongest synergy edges
  S_normalized[is.na(S_normalized)] <- 0
  S_normalized[!is.finite(S_normalized)] <- 0
  diag(S_normalized) <- 0 # 去掉自环

  # Prevent corner cases where quantile might return NA
  vals <- S_normalized[upper.tri(S_normalized)]
  if (length(vals) == 0 || all(!is.finite(vals)) || all(vals == 0, na.rm = TRUE)) {
      return(as(S_normalized, "dgCMatrix"))
  }
  
  # Check if all values are 0 or NA
  if (all(vals == 0)) {
    S_normalized[,] <- 0
  } else {
  # Calculate the threshold to retain only the top edges (e.g., top 25% strongest edges)
  # This step is crucial for breaking large hubs into smaller hubs
    threshold <- quantile(vals, probs = quantile_cutoff, na.rm = TRUE)
    S_normalized[S_normalized < threshold] <- 0
  }

  # 5. Output the sparse matrix
  dimnames(S_normalized) <- list(enhancers, enhancers)
  return(as(S_normalized, "dgCMatrix"))
}


# [F06] run_leiden_for_promoter()
# Purpose : Cluster an enhancer co-membership graph for one promoter using Leiden and return ranked hub partitions.
# Inputs  : promoter_id, comembership (dgCMatrix), resolution, seed.
# Outputs : tibble(promoter_id, hub_id, enhancers(list), dominance).
# Notes   : dominance is defined as within-graph node strength sum per cluster; document this choice for interpretation.
run_leiden_for_promoter <- function(promoter_id,
                                    comembership,
                                    resolution = 0.5,
                                    seed = 42) {

  g <- graph_from_adjacency_matrix(
    comembership,
    mode     = "undirected",
    weighted = TRUE
  )

  enh_names <- V(g)$name

  # If there are too few enhancers or no edges: treat the entire set as one hub
  if (length(enh_names) == 0 || vcount(g) < 3 || ecount(g) == 0) {
    return(tibble(
      promoter_id = promoter_id,
      hub_id      = paste0(promoter_id, "_hub1"),
      enhancers   = list(enh_names),
      dominance   = 0
    ))
  }

  set.seed(seed)
  cl <- leiden(g, resolution_parameter = resolution)

  cl_df <- tibble(
    enhancer_id = V(g)$name,
    cluster     = as.integer(cl)
  )

  ## ===== NEW: dominance = cluster total node strength (internal + external) =====
  # Weighted degree (strength) for each node
  node_strength <- igraph::strength(g, vids = V(g), weights = E(g)$weight)

  # Sum node strength within each cluster
  hub_strength <- cl_df %>%
    mutate(node_strength = node_strength[enhancer_id]) %>%  # relies on vertex names
    group_by(cluster) %>%
    summarise(
      dominance = sum(node_strength, na.rm = TRUE),
      .groups   = "drop"
    )

  ## ===== Original hubs construction =====
  hubs <- cl_df %>%
    group_by(cluster) %>%
    summarise(
      promoter_id = promoter_id,
      enhancers   = list(enhancer_id),
      .groups     = "drop"
    ) %>%
    left_join(hub_strength, by = "cluster") %>%
    mutate(
      dominance = ifelse(is.na(dominance), 0, dominance)
    ) %>%
    arrange(desc(dominance)) %>%  # ranks main hub first under this dominance definition
    mutate(
      hub_rank = row_number(),
      hub_id   = paste0(promoter_id, "_hub", hub_rank)
    ) %>%
    select(promoter_id, hub_id, enhancers, dominance)

  hubs
}

# [F07] build_hubs_from_EP()
# Purpose : End-to-end hub construction from EP interactions: candidate promoter filtering, co-membership, Leiden partitioning, and metrics.
# Inputs  : EP (loops table with weight_raw/weight_percentage), k_min, resolution, seed, quantile_cutoff, method, use_leiden.
# Outputs : list(observed_hubs=hub summary tibble, loops=EP rows restricted to hub members, hubs=hub-style table).
# Notes   : This step is computationally heavy; parallelism is configured globally (future plan/workers) via CLI wrapper.
build_hubs_from_EP <- function(EP,
                               k_min      = 3,
                               resolution = 0.5,
                               seed       = 42,
                               quantile_cutoff = 0,
                               method = "log_zscore",
                               use_leiden = 'on') { #"log_minmax", "log_zscore", "log_maxnorm"
    # Step 1: 每个 promoter 至少连 k_min 个 enhancer 才算候选 hub
    message("Building observed weighted co-membership networks...")
    candidate_hubs <- EP %>%
      group_by(promoter_id) %>%
      filter(n() >= k_min) %>%
      summarise(
        enhancer_id = list(enhancer_id),
        weight_percentage = list(weight_percentage),
        weight_raw  = list(weight_raw),
        n_enh = n(),
        .groups = "drop"
      )
    cat("Step 1: 候选 promoter 数 =", nrow(candidate_hubs), "\n")

    if (nrow(candidate_hubs) == 0) {
      return(tibble(
        promoter_id    = character(0),
        hub_id         = character(0),
        enhancers      = list(),
        size           = integer(0),
        internal_score = numeric(0),
        density        = numeric(0)
      ))
    }

    # Step 2: co-membership 矩阵快速计算函数,hypergeometric-like weighting
    candidate_hubs_comembership <- candidate_hubs %>%
      mutate(comembership = map2(enhancer_id, weight_raw,
                            ~ fast_comembership(.x, .y, quantile_cutoff, method = method)) ##"log_minmax", "log_zscore", "log_maxnorm"
      )

    cat("Fast co-membership", dim(candidate_hubs_comembership))

    #Step 3. 对单个 promoter 的 comembership 矩阵跑 Leiden，得到若干 enhancer hub
    #first promoter test leiden clustering 
    # test_row <- candidate_hubs_comembership %>% slice(1)

    # test_promoter_id <- test_row$promoter_id[[1]]
    # test_mat         <- test_row$comembership[[1]]

    # test_result <- run_leiden_for_promoter(
    #   promoter_id = test_promoter_id,
    #   comembership = test_mat,
    #   resolution   = resolution,
    #   seed         = seed
    # )
    #   cat("single promoter finished.\n")
    #   print(test_result)

    # Run Leiden clustering for all individual promoters in batch
    if (use_leiden=='on') {
      print("Running Leiden clustering for each promoter...")
      candidate_hubs_comembership_subhubs <- candidate_hubs_comembership %>%
        mutate(
          hubs = map2(promoter_id, comembership,
                      ~ run_leiden_for_promoter(.x, .y,
                                                resolution = resolution,
                                                seed       = seed))
        ) %>%
        tidyr::unnest(hubs, names_sep = "_", keep_empty = TRUE)

      cat("Total", nrow(candidate_hubs_comembership_subhubs), "promoter-hub \n")
      print(candidate_hubs_comembership_subhubs, width = Inf)
    } else {
      candidate_hubs_comembership_subhubs <- candidate_hubs_comembership
      candidate_hubs_comembership_subhubs$hubs_hub_id <- candidate_hubs_comembership_subhubs$promoter_id
      candidate_hubs_comembership_subhubs$hubs_enhancers <- candidate_hubs_comembership_subhubs$enhancer_id
      candidate_hubs_comembership_subhubs$hubs_dominance <- candidate_hubs_comembership_subhubs$weight_raw
    }


  print("Running observed counts for each promoter...")
  observed_hubs <- candidate_hubs_comembership_subhubs %>%
    transmute(
      promoter_id,
      hub_id = hubs_hub_id,
      hub_index = regmatches(hubs_hub_id, regexpr("hub[0-9]+", hubs_hub_id)),
      enhancers = hubs_enhancers,
      n_enh = n_enh,
      weight_percentage = weight_percentage,
      weight_raw = weight_raw,
      hubs_dominance = hubs_dominance,
      comembership
    ) %>%
    rowwise() %>%
    mutate(
      size = length(enhancers),

      cm_ok = !is.null(dimnames(comembership)) &&
              !is.null(rownames(comembership)) &&
              !is.null(colnames(comembership)),

      # 1. Total edge weight (Internal Score)
      internal_score = if (cm_ok && size >= 2) {
        cm <- comembership[enhancers, enhancers, drop = FALSE]
        sum(cm) / 2
      } else 0,

      # 2. Classic weighted density
      weighted_density = if (cm_ok && size >= 2) {
        max_edges = size * (size - 1) / 2
        internal_score / max_edges
      } else 0,

      # 3. Average weight of non-zero edges
      avg_non_zero_weight = if (cm_ok && size >= 2) {
        cm <- comembership[enhancers, enhancers, drop = FALSE]
        valid_weights <- cm[upper.tri(cm) & cm > 0]
        if (length(valid_weights) > 0) mean(valid_weights) else 0
      } else 0,

      # 4. Connection ratio (Graph Density)
      graph_density = if (cm_ok && size >= 2) {
        cm <- comembership[enhancers, enhancers, drop = FALSE]
        num_edges = sum(cm[upper.tri(cm)] > 0)
        max_possible_edges = size * (size - 1) / 2
        num_edges / max_possible_edges
      } else 0
    ) %>%
    ungroup() 

    # Step 4: Build full_hubs_merge, which contains all enhancer-promoter relationships
    full_hubs <- observed_hubs %>%
      unnest(enhancers)

    full_hubs_merge <- EP %>%
      right_join(
        full_hubs,
        by = c(
          "promoter_id" = "promoter_id",
          "enhancer_id" = "enhancers"
        )
      )

    hubs=loops_convert(full_hubs_merge)
    hubs[is.na(hubs$hub_index),"hub_index"] <- "hub1"
    hubs[is.na(hubs$hub_id),"hub_id"] <- paste0(hubs[is.na(hubs$hub_id),"promoter_id"],"_hub1")

    cat("Step 8: 最终检测到", nrow(observed_hubs), "个 MEI hubs\n")
    print(head(observed_hubs))                               
    return(list(observed_hubs = observed_hubs, loops = full_hubs_merge, hubs = hubs) )
  }
  

# [F08] PEhub_full()
# Purpose : Convenience orchestrator for the core MEI discovery pipeline (weights -> observed hubs).
# Inputs  : raw_data (preprocessed loops), weight_mode, dist_col, q_col, breaks, alpha, sig_cap, sig_beta, k_min, resolution, seed, quantile_cutoff, method.
# Outputs : list from build_hubs_from_EP() for the observed data.
# Notes   : Permutation/null testing and post-processing are implemented in downstream steps; keep this function side-effect free.
PEhub_full <- function(raw_data, weight_mode = c("count_only"), dist_col = c("D"), q_col = "qvalue", breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6), alpha = 0.8, sig_cap = 10, sig_beta = 0.5, n_perm = 1000, k_min = 3, resolution = 0.1, seed = 42, quantile_cutoff = 0, method = "log_zscore") {
  # raw_data = prep_hichip_weights$loops
  # n_perm = 10
  # k_min = 2
  # resolution = 1
  # seed = 42
  # quantile_cutoff = 0
  # method = "log_zscore" #"log_minmax", "log_zscore", "log_maxnorm"

  print("Starting MEI discovery pipeline:")
  print(weight_mode)
  # 1. 计算 weights
  EP <- compute_weights(raw_data,weight_mode = weight_mode, dist_col = dist_col, q_col = q_col, breaks = breaks, alpha = alpha, sig_cap = sig_cap, sig_beta = sig_beta)

  set.seed(seed)
  
  ## ====== 在 observed EP 上跑一次，得到 observed_hubs_per_promoter ======

  observed_hubs_per_promoter <- build_hubs_from_EP(
    EP         = EP$loops,
    k_min      = k_min,
    resolution = resolution,
    seed       = seed,
    quantile_cutoff = quantile_cutoff,
    method = method,
    use_leiden = 'on'
  )

  cat("Observed hubs:", nrow(observed_hubs_per_promoter), "\n")
  return(observed_hubs_per_promoter)
  }


# [F09] jaccard_threshold_by_size()
# Purpose : Size-aware Jaccard threshold heuristic for defining 'reproducible' hubs under bootstrap.
# Inputs  : size (integer hub size).
# Outputs : numeric threshold in [0,1].
# Notes   : This is a policy choice; adjust thresholds based on empirical stability calibration and report explicitly.
jaccard_threshold_by_size <- function(size) {
  if (size <= 3) return(0.7)  # Small hubs require higher consistency
  if (size <= 5) return(0.6)
  return(0.5)                 # Larger hubs allow some flexibility
}

# [F10] calculate_jaccard()
# Purpose : Compute Jaccard similarity between two enhancer sets.
# Inputs  : set1, set2 (vectors of enhancer IDs).
# Outputs : numeric Jaccard index in [0,1].
# Notes   : Empty-set behavior is defined explicitly (returns 0 when both empty).
calculate_jaccard <- function(set1, set2) {
  # If both sets are empty (shouldn't happen in theory), define as 0 or 1 based on logic; here we use 0
  if (length(set1) == 0 && length(set2) == 0) return(0)
  
  inter <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  
  return(inter / union)
}

# [F11] calculate_hub_stability()
# Purpose : Bootstrap-based stability assessment for observed hubs using Poisson-resampled counts and re-discovery.
# Inputs  : raw_data, observed_hubs, weight_mode, dist_col, q_col, breaks, alpha, sig_cap, sig_beta, B, subsample_frac, k_min, resolution, quantile_cutoff, method, seed.
# Outputs : observed_hubs table augmented with jaccard_stability_score, reproducibility_rate, existence_rate.
# Notes   : B, subsample_frac, and thresholds (tau) materially affect results; expose via CLI and log in a run manifest.
calculate_hub_stability <- function(raw_data, observed_hubs, weight_mode = c("count_only"), dist_col = c("D"), q_col = "qvalue", breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6), alpha = 0.8, sig_cap = 10, sig_beta = 0.5, 
                                    B = 10, subsample_frac = 0.8,
                                    k_min = 3, resolution = 0.5,
                                    quantile_cutoff = 0, method = "log_zscore",
                                    seed = 42) {
  # raw_data = prep_hichip$loops
  # dist_col = c("D")
  # q_col = "qvalue"
  # breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6)
  # alpha = 0.8
  # sig_cap = 10
  # sig_beta = 0.5
  # B = 10
  # subsample_frac = 0.8
  # n_perm = 10
  # k_min = 2
  # resolution = 1
  # seed = 42
  # quantile_cutoff = 0
  # method = "log_minmax" #"log_minmax", "log_zscore", "log_maxnorm"
  # weight_mode =weight_method
  
  # Bootstrap runs
  message(paste0("Running ", B, " bootstrap iterations..."))

  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent pts: :elapsedfull ETA: :eta",
    total = B, clear = FALSE, width = 60)

  stability_runs <- map(1:B, function(i) {
    print(paste0("Bootstrap iteration ", i, "/", B))
    pb$tick() 
    
    set.seed(seed + i)
    EP_bs <- raw_data %>% 
      mutate(counts = rpois(n(), lambda = pmax(0, subsample_frac * counts))) 
    
    EP_bs <- compute_weights(EP_bs,weight_mode = weight_mode, dist_col = dist_col, q_col = q_col, breaks = breaks, alpha = alpha, sig_cap = sig_cap, sig_beta = sig_beta)

    run_res <- build_hubs_from_EP(EP_bs$loops, k_min=k_min, resolution=resolution, seed=seed + i, quantile_cutoff=quantile_cutoff, method=method, use_leiden = 'on')

    rm(EP_bs)

    if (i %% 5 == 0) gc() 

    if (is.null(run_res)) return(NULL)

    valid_hubs <- run_res$observed_hubs %>% 
      filter(!is.null(enhancers), lengths(enhancers) > 0)
    
    return(split(valid_hubs$enhancers, valid_hubs$promoter_id))
  })


  valid_runs <- compact(stability_runs)
  n_valid <- length(valid_runs)
  
  # 3. Calculate stability scores
  scored_results <- observed_hubs %>%
    mutate(
      stability_metrics = map2(promoter_id, enhancers, function(p_id, ref_set) {
        # Get the performance of this promoter across all valid iterations
        tau <- 0.5 # Jaccard threshold for reproducibility

        js <- map_dbl(valid_runs, function(m) {
          cand_sets <- m[[p_id]] # Get all candidate hubs for this promoter in the current iteration
          if (is.null(cand_sets)) return(0)
          # Find the most similar hub
          max(vapply(cand_sets, function(x) calculate_jaccard(x, ref_set), numeric(1)))
        })
        ov_best <- map_int(valid_runs, function(m) {
          cand_sets <- m[[p_id]]
          if (is.null(cand_sets) || length(cand_sets) == 0) return(0L)

          # Find the best matching hub (highest Jaccard similarity)
          jvals <- vapply(cand_sets, function(x) calculate_jaccard(x, ref_set), numeric(1))
          k <- which.max(jvals)

          # Count the overlap size of the best matching hub
          length(intersect(cand_sets[[k]], ref_set))
        })

        list(
          jaccard_stability_score = mean(js, trim = 0.1),
          reproducibility_rate = mean(js >= tau),
          existence_rate = mean(ov_best >= 2)
        )
      })) %>%
    unnest_wider(stability_metrics)
  print(paste("Hub stability number:", nrow(scored_results[scored_results$reproducibility_rate >= 0.5, ])))
  return(scored_results)
}



# [F12] build.null.pvalue()
# Purpose : Prepare distance-stratified weight pools and per-hub node/bin lists for global p-value estimation.
# Inputs  : hub_tbl (hub summaries), hub_tbl_ep (EP edges for hubs), hub_tbl_all_ep (genome-wide EP pool), method, weight_method, B, quantile_cutoff, stat.
# Outputs : list(hubs=hub table with w_list/bin_list attached, global_bins=distance-bin pooled weights).
# Notes   : This function prepares inputs only; heavy Monte Carlo is executed by build.null.pvalue.calculate().
build.null.pvalue <- function(hub_tbl, hub_tbl_ep, hub_tbl_all_ep, method = "log_minmax", weight_method = "bin_log_ratio_sig", B = 1000, quantile_cutoff = 0, stat = "density") {
  print("Building null model for hub p-value calculation...")
  # B <- 100
  # method <- "log_minmax"      
  # quantile_cutoff <- 0        # Should match the fast_comembership setting
  # stat <- "density"           # Recommended; use "sum" if you prefer total strength

  ## =========================================================
  ## 1) Distance binning (global pool bins)
  ## =========================================================
  # Define distance bins
  breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6, Inf)

  ## Note: The EP table may use either "D" or "distance" as the column name.
  ## Here, we assume "D". If your column is named "distance", replace "D" with "distance".

  ## =========================================================
  ## 2) Select the weight column for null and hub score calculation
  ## =========================================================
  ## Use "weight_raw" for comembership calculations. The column name might vary
  ## (e.g., "weight_raw.x" or "weight_raw.y" after a join).
  ## The function below ensures robust selection of the correct column.
  pick_weight_col <- function(df) {
    if ("weight_raw.x" %in% names(df)) return("weight_raw.x")
    if ("weight_raw"   %in% names(df)) return("weight_raw")
    if ("weight_raw.y" %in% names(df)) return("weight_raw.y")
    stop("Cannot find weight_raw column in the object (tried weight_raw.x / weight_raw / weight_raw.y).")
  }
  wcol_all <- pick_weight_col(hub_tbl_all_ep)
  wcol <- pick_weight_col(hub_tbl_ep)

  ## =========================================================
  ## 3) Build a global distance-stratified weight pool (node level)
  ## =========================================================
  global_bins <- hub_tbl_all_ep %>%
    filter(pair_type == "PE") %>%                 # Use only EP interactions
    filter(is.finite(D), D >= 0) %>%
    mutate(dist_bin = cut(D, breaks = breaks, include.lowest = TRUE, right = TRUE)) %>%
    group_by(dist_bin) %>%
    summarise(weight_pool = list(.data[[wcol_all]]), .groups = "drop")
  ## Safety check: Ensure each bin has enough values
  if (any(lengths(global_bins$weight_pool) < 50)) {
    message("Warning: some distance bins have small pool size (<50). Consider merging bins or using fewer breaks.")
  }

  ## =========================================================
  ## 4) Prepare dist_bins for each hub (extract enhancer distances from EP table)
  ## =========================================================
  ## Extract hub_id + enhancer_id + D + weight_raw from hub_tbl_ep
  ## Deduplicate and group by hub to align with hub_tbl enhancers
  hub_nodes <- hub_tbl_ep %>%
    select(promoter_id, enhancer_id, D, w = all_of(wcol)) %>%
    distinct(promoter_id, enhancer_id, .keep_all = TRUE) %>%
    mutate(dist_bin = cut(D, breaks = breaks, include.lowest = TRUE, right = TRUE)) %>%
    group_by(promoter_id) %>%
    summarise(
      enh_list = list(enhancer_id),
      w_list   = list(w),
      bin_list = list(dist_bin),
      .groups  = "drop"
    )

  ## Merge hub_nodes back into hub_tbl
  hub_tbl2 <- hub_tbl %>%
    select(promoter_id, hub_id, enhancers, weight_raw, size, cm_ok, weighted_density, internal_score) %>%
    left_join(hub_nodes, by = c("promoter_id"))

  ## =========================================================
  ## 6) Compute global p-values (recommended stat: density, aligns with weighted_density)
  ## =========================================================
  ## Step 1: Flags
  hub_tbl3 <- hub_tbl2 %>%
    mutate(
      has_w   = map_lgl(w_list,   ~ !is.null(.x) && length(.x) >= 2),
      has_bin = map_lgl(bin_list, ~ !is.null(.x) && length(.x) >= 2),
      has_inputs = has_w & has_bin
    )

  return(list(
    hubs = hub_tbl3,
    global_bins = global_bins
  ))
}

# [F13] hub_score_from_weights()
# Purpose : Compute a hub score from node weights by constructing a co-membership matrix and summarizing its edge weights.
# Inputs  : w (numeric), method, quantile_cutoff, stat ('density' or 'sum').
# Outputs : numeric hub score (density-normalized or raw sum).
# Notes   : Must match the scoring used in observed discovery; otherwise p-values are not interpretable.
hub_score_from_weights <- function(w,
                                  method = "log_minmax",
                                  quantile_cutoff = 0,
                                  stat = c("density", "sum")) {
  stat <- match.arg(stat)
  w <- as.numeric(w)
  w[!is.finite(w)] <- 0
  # w <- pmax(w, 0)
  n <- length(w)
  if (n < 2) return(0)

  enh <- paste0("e", seq_len(n))

  cm <- fast_comembership(
    enhancers = enh,
    weights = w,
    quantile_cutoff = quantile_cutoff,
    method = method
  )

  s <- as.numeric(sum(cm) / 2)  
  if (stat == "sum") return(s)
  return(s / (n * (n - 1) / 2)) # density
}

# [F14] calculate_global_p()
# Purpose : Monte Carlo p-value for a hub score under a chosen null sampling scheme (distance-matched / histogram-matched / global).
# Inputs  : observed_weights, dist_bins, pool_map, B, method, quantile_cutoff, stat, null_mode.
# Outputs : list(p, obs, null_median).
# Notes   : Uses a smoothed p-value to avoid 0; ensure pool_map sizes are adequate per distance bin.
calculate_global_p <- function(observed_weights,
                 dist_bins,
                 pool_map,
                 B = 10000,
                 method = "log_minmax",
                 quantile_cutoff = 0,
                 stat = c("density", "sum"),
                 null_mode = c("distance_matched", "hist_matched", "global")) {
  stat <- match.arg(stat)
  null_mode <- match.arg(null_mode)

  # Calculate the observed hub score
  obs_score <- hub_score_from_weights(
  observed_weights,
  method = method,
  quantile_cutoff = quantile_cutoff,
  stat = stat
  )

  # Map each node/bin to its corresponding pool (dist_bins might be a factor)
  key <- as.character(dist_bins)

  if (null_mode == "distance_matched") {
  # Strict mode: each node gets its own pool (length == number of nodes)
  pools <- pool_map[key]
  if (any(vapply(pools, is.null, logical(1)))) {
    return(list(p = NA_real_, obs = obs_score, null_median = NA_real_))
  }

  # Generate null scores by sampling one weight per node from its pool
  null_scores <- replicate(B, {
    w_star <- vapply(pools, function(pool) sample(pool, 1), numeric(1))
    hub_score_from_weights(
    w_star,
    method = method,
    quantile_cutoff = quantile_cutoff,
    stat = stat
    )
  })

  } else if (null_mode == "hist_matched") {
  # Relaxed mode: match the hub's distance-bin histogram
  bin_counts <- table(key)
  bins <- names(bin_counts)

  pools_by_bin <- pool_map[bins]
  if (any(vapply(pools_by_bin, is.null, logical(1)))) {
    return(list(p = NA_real_, obs = obs_score, null_median = NA_real_))
  }
  counts <- as.integer(bin_counts)

  # Generate null scores by sampling weights based on the bin histogram
  null_scores <- replicate(B, {
    w_star <- unlist(
    Map(function(pool, k) sample(pool, size = k, replace = TRUE),
      pools_by_bin, counts),
    use.names = FALSE
    )

    hub_score_from_weights(
    w_star,
    method = method,
    quantile_cutoff = quantile_cutoff,
    stat = stat
    )
  })
  } else if (null_mode == "global") {
  # Global mode: sample weights from the entire global pool
  global_pool <- unlist(pool_map[-(1:2)], use.names = FALSE)
  
  # Determine the number of nodes in the current hub
  n_nodes <- length(observed_weights)
  
  # Generate null scores by sampling weights from the global pool
  null_scores <- replicate(B, {
    w_star <- sample(global_pool, size = n_nodes, replace = TRUE)
    
    hub_score_from_weights(
    w_star,
    method = method,
    quantile_cutoff = quantile_cutoff,
    stat = stat
    )
  })
  }

  # Calculate the smoothed p-value to avoid zero
  p_val <- (1 + sum(null_scores >= obs_score, na.rm = TRUE)) / (1 + length(null_scores))

  list(
  p = p_val,
  obs = obs_score,
  null_median = median(null_scores, na.rm = TRUE)
  )
}



# [F15] build.null.pvalue.calculate()
# Purpose : Parallelized computation of global p-values for all hubs using the prepared pools and per-hub inputs.
# Inputs  : hub_tbl3 (from build.null.pvalue()), global_bins, B, method, quantile_cutoff, stat, null_mode.
# Outputs : hub table with hub_p_value_global, hub_p_adj_global (BH), observed/null scores, and OE_ratio_global.
# Notes   : The parallel backend and worker count must be configured via CLI; monitor memory due to large pools.
build.null.pvalue.calculate <- function(hub_tbl3, global_bins, B = 1000, method = "log_minmax", quantile_cutoff = 0, stat = "sum", null_mode = "distance_matched") {
  ## Step 2: heavy compute (res_list)
  print("Start the heavy step, calculating global p-values for each hub...")
  pool_map <- global_bins$weight_pool
  names(pool_map) <- as.character(global_bins$dist_bin)
  shrink_pool_map <- function(pool_map, max_per_bin = 1000000, seed = 1) {
    set.seed(seed)
    lapply(pool_map, function(v) {
      v <- as.numeric(v)
      v <- v[is.finite(v)]
      if (length(v) <= max_per_bin) return(v)
      sample(v, size = max_per_bin, replace = FALSE)
    })
  }

  pool_map_small <- shrink_pool_map(pool_map, max_per_bin = 1000000, seed = 1)

  res_list <- furrr::future_map2(
    hub_tbl3$w_list,
    hub_tbl3$bin_list,
    function(w_vec, bin_vec) {
      ok <- !is.null(w_vec) && length(w_vec) >= 2 &&
            !is.null(bin_vec) && length(bin_vec) >= 2
      if (!ok) return(list(p = NA_real_, obs = NA_real_, null_median = NA_real_))

      calculate_global_p(
        observed_weights = w_vec,
        dist_bins        = bin_vec,
        # global_bins      = global_bins,
        pool_map         = pool_map_small,
        B                = B,
        method           = method,
        quantile_cutoff  = quantile_cutoff,
        stat             = stat,
        null_mode = null_mode
      )
    },
    .options = furrr::furrr_options(
      seed = TRUE,
      packages = c("Matrix"),
        globals = list(
          calculate_global_p = calculate_global_p,
          pool_map            = pool_map_small,
          B                   = B,
          method              = method,
          quantile_cutoff     = quantile_cutoff,
          stat                = stat,
          null_mode = null_mode,
          hub_score_from_weights = hub_score_from_weights,
          fast_comembership   = fast_comembership
        )
    )
  )
  print("Finished calculating global p-values.")

  hub_tbl_p <- hub_tbl3 %>%
    mutate(
      .res = res_list,
      hub_p_value_global    = map_dbl(.res, "p"),
      hub_score_obs_global  = map_dbl(.res, "obs"),
      hub_score_null_median = map_dbl(.res, "null_median"),
      OE_ratio_global = ifelse(
        is.finite(hub_score_null_median) & hub_score_null_median > 0,
        hub_score_obs_global / hub_score_null_median,
        NA_real_
      ),
      hub_p_adj_global = p.adjust(hub_p_value_global, method = "BH", n = n())
    ) %>%
    select(-.res, -has_w, -has_bin, -has_inputs)

  ## 你也可以看看显著 hub 的数量
  print(paste0("Pvalue num FDR<0.05 hubs: ", sum(hub_tbl_p$hub_p_adj_global < 0.05, na.rm = TRUE)))
  
  return(hub_tbl_p)
}

# [F16] build.run.preprocess()
# Purpose : CLI-facing wrapper to preprocess inputs for a sample and persist intermediate RData for downstream steps.
# Inputs  : loop_file (significant interactions; CLI), loop_file_all (all interactions; CLI), outdir (CLI), tss_input, filename (sample ID), promoter_window.
# Outputs : Writes '<outdir>/multiple_result.<filename>.hub.all.preprocess.RData'.
# Notes   : This function performs file I/O; ensure outdir exists and is writable before invoking via Rscript.
build.run.preprocess <- function(loop_file, loop_file_all, outdir, tss_input, filename, promoter_window=0) {
  prep_hichip_all <- preprocess_hichip(loop_file = loop_file, tss_input = tss_input, mode = "all", promoter_window = promoter_window)
  prep_hichip <- preprocess_hichip(loop_file = loop_file, tss_input = tss_input, mode = "promoter", promoter_window = promoter_window)
  all_ep_data = preprocess_hichip(loop_file = loop_file_all, tss_input = tss_input, mode = "promoter", promoter_window = promoter_window)

  save(prep_hichip_all, prep_hichip, all_ep_data, file = file.path(outdir, paste("multiple_result",filename,"hub","all.preprocess.RData",sep=".")))
  print("Preprocessing done.")
}

# [F17] build.run.all()
# Purpose : CLI-facing wrapper to run observed hub discovery for one (method, weight_method) setting and persist preprocess objects for cutoff.
# Inputs  : outdir (CLI), method, weight_method, filename, promoter_window.
# Outputs : Writes '<outdir>/multiple_result.<filename>.hub.<method>.<weight_method>.preprocess.RData'.
# Notes   : Requires that build.run.preprocess() has been executed for the same filename/outdir.
build.run.all <- function(outdir, method, weight_method, filename, promoter_window=0) {
  load(file.path(outdir, paste("multiple_result",filename,"hub","all.preprocess.RData",sep=".")))
  observed_hubs_per_promoter <- PEhub_full(raw_data = prep_hichip$loops, weight_mode = weight_method, dist_col = c("D"), q_col = "qvalue", breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6), alpha = 0.8, sig_cap = 10, sig_beta = 0.5, n_perm = 10, resolution = resolution, k_min = k_min, seed = 42, quantile_cutoff = quantile_cutoff, method = method)
  observed_hubs_per_promoter_sub <- observed_hubs_per_promoter$loops[!is.na(observed_hubs_per_promoter$loops$hub_index) & observed_hubs_per_promoter$loops$hub_index=="hub1",]
  observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter$hubs[!is.na(observed_hubs_per_promoter$hubs$hub_index) & observed_hubs_per_promoter$hubs$hub_index=="hub1",]
  observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_sub_hub[order(observed_hubs_per_promoter_sub_hub$hubid),]

  #pvalue prepare
  all_ep_data_weight <- compute_weights(all_ep_data$loops, weight_mode = weight_method, dist_col = c("D"), q_col = "qvalue", breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6), alpha = 0.8, sig_cap = 10, sig_beta = 0.5)
  observed_hubs_per_promoter_sub_observed <- observed_hubs_per_promoter$observed_hubs[observed_hubs_per_promoter$observed_hubs$promoter_id %in% observed_hubs_per_promoter_sub$promoter_id,]
  observed_hubs_per_promoter_pvalue_prepare <- build.null.pvalue(hub_tbl=observed_hubs_per_promoter_sub_observed, hub_tbl_ep=observed_hubs_per_promoter_sub, hub_tbl_all_ep=all_ep_data_weight$loops, method = method, B = B_pvalue, quantile_cutoff = quantile_cutoff, stat = "density")

  save(observed_hubs_per_promoter, observed_hubs_per_promoter_sub, observed_hubs_per_promoter_sub_hub, observed_hubs_per_promoter_pvalue_prepare, file = file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"preprocess.RData",sep=".")))
  print("All runs done.")
}

# [F18] build.cutoff()
# Purpose : CLI-facing wrapper to compute hub stability metrics and global p-values, then persist merged results.
# Inputs  : outdir (CLI), method, weight_method, filename.
# Outputs : Writes '<outdir>/multiple_result.<filename>.hub.<method>.<weight_method>.RData'.
# Notes   : Uses both bootstrap (B) and Monte Carlo (B) parameters hardcoded in calls; expose via CLI in your wrapper.
build.cutoff <- function(outdir, method, weight_method, filename) {
  print(paste("Starting post-cutoff:", filename))
  load(file.path(outdir, paste("multiple_result",filename,"hub","all.preprocess.RData",sep=".")))
  load(file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"preprocess.RData",sep=".")))
  
  # #stability 
  observed_hubs_per_promoter_stability_all <- calculate_hub_stability(raw_data = prep_hichip$loops, observed_hubs = observed_hubs_per_promoter$observed_hubs[!is.null(observed_hubs_per_promoter$observed_hubs$enhancers) & observed_hubs_per_promoter$observed_hubs$hub_index=="hub1",], weight_mode = weight_method, dist_col = c("D"), q_col = "qvalue", breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6), alpha = 0.8, sig_cap = 10, sig_beta = 0.5, 
    B = B_stability, subsample_frac = 0.8,
    k_min = k_min, resolution = resolution,
    quantile_cutoff = quantile_cutoff, method = method,
    seed = 42)
  observed_hubs_per_promoter_stability <- observed_hubs_per_promoter_stability_all %>% select(promoter_id, hub_id, jaccard_stability_score, reproducibility_rate, existence_rate) # %>% filter(!is.na(hub_p_adj_global), hub_p_adj_global <= 0.05)
  observed_hubs_per_promoter_sub <- observed_hubs_per_promoter_sub %>% left_join( observed_hubs_per_promoter_stability, by = c("promoter_id", "hub_id") )
  observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_sub_hub %>% left_join( observed_hubs_per_promoter_stability, by = c("promoter_id", "hub_id") )
  #save(prep_hichip_all, observed_hubs_per_promoter, observed_hubs_per_promoter_sub, observed_hubs_per_promoter_sub_hub, observed_hubs_per_promoter_stability, file = file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"stability.RData",sep=".")))
  #observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_stability %>%  filter(reproducibility_rate >= 0.5 & existence_rate >= 0.5)

  #pvalue
  observed_hubs_per_promoter_pvalue_all <- build.null.pvalue.calculate(hub_tbl3=observed_hubs_per_promoter_pvalue_prepare$hubs, global_bins=observed_hubs_per_promoter_pvalue_prepare$global_bins, method = method, B = B_pvalue, quantile_cutoff = quantile_cutoff, stat = "sum", null_mode = null_mode)
  observed_hubs_per_promoter_pvalue <- observed_hubs_per_promoter_pvalue_all %>% select(promoter_id, hub_id, hub_p_value_global, hub_p_adj_global, hub_score_obs_global, hub_score_null_median, OE_ratio_global) # %>% filter(!is.na(hub_p_adj_global), hub_p_adj_global <= 0.05)
  observed_hubs_per_promoter_sub <- observed_hubs_per_promoter_sub %>% left_join( observed_hubs_per_promoter_pvalue, by = c("promoter_id", "hub_id") )
  observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_sub_hub %>% left_join( observed_hubs_per_promoter_pvalue, by = c("promoter_id", "hub_id") )
  save(prep_hichip_all, prep_hichip, observed_hubs_per_promoter, observed_hubs_per_promoter_sub, observed_hubs_per_promoter_sub_hub, observed_hubs_per_promoter_stability, observed_hubs_per_promoter_pvalue, file = file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"RData",sep=".")))
}

# [F19] build.postprocess()
# Purpose : Export final significant hubs and remaining pairwise interactions to BED/BEDPE/TXT deliverables.
# Inputs  : outdir (CLI), method, weight_method, filename.
# Outputs : Writes hub and pairwise exports under '<outdir>/' with prefixes 'multiple_result' and 'pairwise'.
# Notes   : Selection thresholds (FDR, reproducibility_rate, etc.) are policy choices; keep them configurable via CLI wrapper.
build.postprocess <- function(outdir, method, weight_method, filename) {
  print(paste("Starting post-processing:", filename))
  load(file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"RData",sep=".")))
  
  #export significant hubs
  observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_sub_hub %>% filter(!is.na(hub_p_adj_global), hub_p_adj_global <= pvalue_cutoff, reproducibility_rate >= stability_cutoff) 
  print(paste0("Significant hubs number: ", nrow(observed_hubs_per_promoter_sub_hub)))
  observed_hubs_per_promoter_sub_hub$name <- with(observed_hubs_per_promoter_sub_hub[order(observed_hubs_per_promoter_sub_hub$hubid),], { h <- gsub("_","", observed_hubs_per_promoter_sub_hub$hubid); paste0(observed_hubs_per_promoter_sub_hub$promoter_gene_name, "_", h, "_", ave(h, h, FUN = seq_along)) })
  gr <- makeGRangesFromDataFrame(observed_hubs_per_promoter_sub_hub,keep.extra.columns = TRUE)
  export(gr, con=file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"bed",sep=".")))
  write.table(as.data.frame(observed_hubs_per_promoter_sub_hub)[,c(1:21,ncol(observed_hubs_per_promoter_sub_hub))], file=file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"txt",sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
  colnames(observed_hubs_per_promoter_sub)[1] <- paste0("#",colnames(observed_hubs_per_promoter_sub)[1])
  write.table(observed_hubs_per_promoter_sub[,1:12], file=file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"bedpe",sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
  
  #pairwise
  observed_hubs_per_promoter_sub_hub_pairs <- prep_hichip$hubs[!prep_hichip$hubs$hubid %in% observed_hubs_per_promoter_sub_hub$hubid,]
  observed_hubs_per_promoter_sub_pairs <- prep_hichip$loops[!prep_hichip$loops$hubid %in% observed_hubs_per_promoter_sub$hubid,]
  gr_pairs <- makeGRangesFromDataFrame(observed_hubs_per_promoter_sub_hub_pairs,keep.extra.columns = TRUE)
  gr_pairs$name <- gr_pairs$hubid
  export(gr_pairs, con=file.path(outdir, paste("multiple_result",filename,"pairwise",method,weight_method,"bed",sep=".")))
  write.table(as.data.frame(observed_hubs_per_promoter_sub_hub_pairs), file=file.path(outdir, paste("multiple_result",filename,"pairwise",method,weight_method,"txt",sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
  colnames(observed_hubs_per_promoter_sub_pairs)[1] <- paste0("#",colnames(observed_hubs_per_promoter_sub_pairs)[1])
  write.table(observed_hubs_per_promoter_sub_pairs[,1:12], file=file.path(outdir, paste("multiple_result",filename,"pairwise",method,weight_method,"bedpe",sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
  
  print("Postprocessing done.")
}


##########################
##Run end-to-end workflow
##########################
build.run.preprocess(loop_file = loop_file, loop_file_all = loop_file_all, 
              outdir = outdir, tss_input = promoters_annotation, filename = filename, promoter_window = promoter_window)
build.run.all(outdir = outdir, 
              method = method, weight_method = weight_method, filename = filename, promoter_window = promoter_window)
build.cutoff(outdir = outdir,  
              method = method, weight_method = weight_method, filename = filename)
build.postprocess(outdir = outdir,  
              method = method, weight_method = weight_method, filename = filename)
