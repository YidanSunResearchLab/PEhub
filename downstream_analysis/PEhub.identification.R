# ==================== PEhub v1.0 ====================
# 安装依赖（一次性）
# BiocManager::install(c("GenomicRanges", "dplyr", "purrr", "igraph", "parallel"))

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
plan(multisession, workers = 10)

##leiden clustering
library(reticulate)
use_python("/mnt/citadel2/research/syidan/miniconda3/envs/r_porc/bin/python", required = TRUE)
py_config() 
library(leiden)

##Input
tss_input = "~/syidan/Genomes/GRCh38/release-47-index/annotation/genes.gtf"


##Prepare gene annotation: promoter regions ±2kb around TSS
# # 支持四种输入模式：TxDb / orgdb / gtf / bed
# ## TxDb 模式：直接从 TxDb.Hsapiens.UCSC.hg38.knownGene 取 promoter
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
gtf <- rtracklayer::import(tss_input)
promoters_annotation <- gtf %>%
  subset(type == "transcript") %>%
  unique()



# ------------------- Step 0: preprocess the hichip bedpe input -------------------

loops_convert <- function(dt) {

  setDT(dt)  # 不会深拷贝到每一列

  # rest_df: 等价于 loops[, -(1:3)]，并把前3列改名 chr/start/end
  rest <- dt[, -(1:3)]
  data.table::setnames(rest, old = names(rest)[1:3], new = c("chr", "start", "end"))

  # anchor_df: 只取需要的列
  anchor_cols <- c(
    "chr1","start1","end1","type1","type2","pair_type",
    "promoter_id","promoter_gene_id","promoter_transcript_id",
    "promoter_gene_name","hubid"
  )
  # 只取存在的列（避免缺列报错）
  anchor_cols <- intersect(anchor_cols, names(dt))

  anchor <- dt[, ..anchor_cols]
  data.table::setnames(anchor, old = c("chr1","start1","end1"), new = c("chr","start","end"), skip_absent = TRUE)

  # 强制 anchor 的 pair_type 为 PP（不存在就新建）
  anchor[, pair_type := "PP"]

  # 每个 hubid 只保留第一条（比 duplicated 更快/更省）
  if ("hubid" %in% names(anchor)) {
    data.table::setkey(anchor, hubid)
    anchor <- anchor[!duplicated(hubid)]
  }

  # 给 anchor 补齐 rest 中缺失的列（一次性，不要循环逐列赋值）
  miss <- setdiff(names(rest), names(anchor))
  if (length(miss)) {
    anchor[, (miss) := NA]
  }

  # 让 anchor 列顺序与 rest 一致
  data.table::setcolorder(anchor, names(rest))

  # rbind（data.table 的 rbindlist 比 base rbind 快）
  final <- data.table::rbindlist(list(rest, anchor), use.names = TRUE, fill = TRUE)

  # enhancer_id NA -> promoter_id（用 data.table 的 fcoalesce 更快）
  if (all(c("enhancer_id", "promoter_id") %in% names(final))) {
    final[is.na(enhancer_id), enhancer_id := promoter_id]
  }

  # 如果你真的需要按 hubid 排序（很耗时，能省则省）
  if ("hubid" %in% names(final)) {
    data.table::setorder(final, hubid)
  }

  return(as.data.frame(final))
}


make_peak_id <- function(chr, start, end, scale = 1000L) {
  paste0("peak_",chr, "_", start %/% scale, "_", end %/% scale)
}

preprocess_hichip <- function(
    loop_file,
    tss_input = NULL,                  # GRanges object，include gene_id column
    mode = c("all"),       # all = EP+PP, promoter = 只EP
    promoter_window = 0,
    h3k27ac_peak = NULL,
    # h3k27ac_bw = NULL,
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
  
  # 0. read loops
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

  # 4. 关键开关：选择保留哪些类型
  if (mode == "promoter") {
    loops <- loops %>% filter(pair_type %in% c("PE"))
    message("Only keep enhancer promoter interactions, loop number: ", nrow(loops))
  } else {
    loops <- loops %>% filter(pair_type %in% c("PE","EE"))
    message("Keep all PE and EE interaction types, loop number: ", nrow(loops))
  }
    
  # 5. 计算距离 + 过滤
  loops <- loops %>%
    mutate(distance = abs((start1+end1)/2 - (start2+end2)/2)) %>%
    filter(distance >= min_dist, distance <= max_dist)
  message("After distance filtering, loop number: ", nrow(loops))
  
  # 6. 是否与 H3K27ac peaks overlap 过滤
  if(!is.null(h3k27ac_peak)){
    ## 1. 读 H3K27ac peaks
    h3k27ac_peaks <- rtracklayer::import(h3k27ac_peak)
    
    ## 2. 判断是否与 h3k27ac_peaks overlap
    gr_left  <- GRanges(loops$chr1, IRanges(loops$start1 + 1, loops$end1))
    gr_right <- GRanges(loops$chr2, IRanges(loops$start2 + 1, loops$end2))

    over_left  <- overlapsAny(gr_left,  h3k27ac_peaks + 5000)
    over_right <- overlapsAny(gr_right, h3k27ac_peaks + 5000)

    ## 3. 至少一端 overlap 就保留
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
  
  # 7. make hub object
  hubs <- loops_convert(loops)

  return(list(loops=as.data.frame(loops), hubs=hubs, peaks=peaks))

}

# ------------------- Step 1: Compute weights for each EP interaction in 11 different ways -------------------
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


# ------------------- Step 2: Fast co-membership matrix computation -------------------
fast_comembership <- function(enhancers, weights, quantile_cutoff = 0, method = c("log_minmax")) { #"log_minmax", "log_zscore", "log_maxnorm"
  # 1. 准备数据
  w <- as.numeric(weights)
  w[!is.finite(w)] <- 0
  w <- pmax(w, 0)
  n <- length(w)
  
  if (n < 2) {
    mat <- Matrix(0, n, n, sparse = TRUE)
    dimnames(mat) <- list(enhancers, enhancers)
    return(mat)
  }
  
  # 2. Log-Space 计算 (完美保留你的 Synergy 公式逻辑)
  # 原公式: Score_ij ≈ w_i * w_j * [ T_all / ((1+w_i)*(1+w_j)) ]
  # 我们忽略 "-1"，因为当 T_all 很大时，-1 几乎不影响相对排序
  
  # 预计算 Log 值，防止溢出
  # 使用 log1p(x) 计算 log(1+x) 更精确
  log_w <- log(w + 1e-10)    # 防止 log(0)
  log_1_plus_w <- log1p(w)   # log(1+w)
  sum_log_T <- sum(log_1_plus_w) # 这是原来的 log(T_all)
  
  # 矩阵化计算 Log(Score)
  # Log(S_ij) = log(w_i) + log(w_j) + log(T_all) - log(1+w_i) - log(1+w_j)
  
  # 为了速度，使用 outer (向量外积)
  term_wi_wj <- outer(log_w, log_w, "+")         # log(w_i) + log(w_j)
  term_denom <- outer(log_1_plus_w, log_1_plus_w, "+") # log(1+w_i) + log(1+w_j)
  
  log_S_matrix <- term_wi_wj + sum_log_T - term_denom
  
  # 3. 归一化 (关键步骤)
  # 找出矩阵中最大的 Log 值，将所有值平移，使最大值为 0 (即原值为 1)
  # 这样保证了 exp() 之后数值在 [0, 1] 之间，且相对比例不变
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
  
  # 4. 阈值化 (Sparsification) - 解决 Hub 过于臃肿的问题
  # 只有原本就很强的 Synergy 边才保留
  S_normalized[is.na(S_normalized)] <- 0
  S_normalized[!is.finite(S_normalized)] <- 0
  diag(S_normalized) <- 0 # 去掉自环

  #防止 quantile 出 NA 的 corner case
  vals <- S_normalized[upper.tri(S_normalized)]
  # 1. 检查是否为空 (防止 n=1 或 n=0 的情况)
  if (length(vals) == 0 || all(!is.finite(vals)) || all(vals == 0, na.rm = TRUE)) {
      # 矩阵太小，没有非对角线元素，直接返回
      return(as(S_normalized, "dgCMatrix"))
  }
  # 2. 检查是否全为 0 或全为 NA
  if (all(vals == 0)) {
    # 所有边都很弱，直接返回全零矩阵
    S_normalized[,] <- 0
  } else {
  # 计算保留的阈值 (例如只保留前 25% 强的边)
  # 注意：这步对于把大 Hub 打散成小 Hub 至关重要
    threshold <- quantile(vals, probs = quantile_cutoff, na.rm = TRUE)
  # 低于阈值的置为 0
    S_normalized[S_normalized < threshold] <- 0
  }

  # 5. 输出稀疏矩阵
  dimnames(S_normalized) <- list(enhancers, enhancers)
  return(as(S_normalized, "dgCMatrix"))
}


# ------------------- Step 3: Leiden clustering for each promoter's co-membership matrix -------------------
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

  # enhancer 太少或没有边：整体一个 hub
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

# ------------------- Step 4: Build enhancer hubs from EP interactions using weights (including comembership and leiden cluster) -------------------
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

    cat("超快版完成！矩阵", dim(candidate_hubs_comembership), "，用时 <10 秒\n")

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
    #   cat("测试单个 promoter 完成，结果：\n")
    #   print(test_result)

    # 对所有单个 promoter 批量跑 Leiden
    if (use_leiden=='on') {
      print("Running Leiden clustering for each promoter...")
      candidate_hubs_comembership_subhubs <- candidate_hubs_comembership %>%
        mutate(
          hubs = map2(promoter_id, comembership,
                      ~ run_leiden_for_promoter(.x, .y,
                                                resolution = resolution,
                                                seed       = seed))
        ) %>%
        # select(promoter_id, hubs) %>%
        tidyr::unnest(hubs, names_sep = "_", keep_empty = TRUE)

      cat("检测到总共", nrow(candidate_hubs_comembership_subhubs), "个 promoter-hub 组合\n")
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

      # 新增：comembership 是否有有效的行列名（最小化改动的核心）
      cm_ok = !is.null(dimnames(comembership)) &&
              !is.null(rownames(comembership)) &&
              !is.null(colnames(comembership)),

      # 1. 总边权 (Internal Score)
      internal_score = if (cm_ok && size >= 2) {
        cm <- comembership[enhancers, enhancers, drop = FALSE]
        sum(cm) / 2
      } else 0,

      # 2. 经典加权密度
      weighted_density = if (cm_ok && size >= 2) {
        max_edges = size * (size - 1) / 2
        internal_score / max_edges
      } else 0,

      # 3. 平均非零边权
      avg_non_zero_weight = if (cm_ok && size >= 2) {
        cm <- comembership[enhancers, enhancers, drop = FALSE]
        valid_weights <- cm[upper.tri(cm) & cm > 0]
        if (length(valid_weights) > 0) mean(valid_weights) else 0
      } else 0,

      # 4. 连通比例 (Graph Density)
      graph_density = if (cm_ok && size >= 2) {
        cm <- comembership[enhancers, enhancers, drop = FALSE]
        num_edges = sum(cm[upper.tri(cm)] > 0)
        max_possible_edges = size * (size - 1) / 2
        num_edges / max_possible_edges
      } else 0
    ) %>%
    ungroup() 

    # Step 4: 构建 full_hubs_merge，包含所有 enhancer - promoter 对应关系
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
  

# ------------------- Step 7: Full MEI discovery pipeline -------------------
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


# ------------------- Step 8: Hub stability assessment via bootstrap -------------------
# size-aware threshold（你可以按数据再调）
# 示例：根据集合大小动态设定阈值
jaccard_threshold_by_size <- function(size) {
  if (size <= 3) return(0.7)  # 小 Hub 要求更高的一致性
  if (size <= 5) return(0.6)
  return(0.5)                 # 大 Hub 允许一定的漂移
}

# 定义 Jaccard 相似度计算函数
calculate_jaccard <- function(set1, set2) {
  # 如果两个都是空的（理论上不应发生），定义为 0 或 1 取决于逻辑，这里给 0
  if (length(set1) == 0 && length(set2) == 0) return(0)
  
  inter <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  
  return(inter / union)
}

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

  # 1. 预分配进度条
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent pts: :elapsedfull ETA: :eta",
    total = B, clear = FALSE, width = 60)

  # 2. 使用 map 进行迭代
  stability_runs <- map(1:B, function(i) {
    print(paste0("Bootstrap iteration ", i, "/", B))
    pb$tick() # 更新进度条
    
    set.seed(seed + i)
    # 执行计算逻辑
    EP_bs <- raw_data %>% 
      mutate(counts = rpois(n(), lambda = pmax(0, subsample_frac * counts))) 
    
    EP_bs <- compute_weights(EP_bs,weight_mode = weight_mode, dist_col = dist_col, q_col = q_col, breaks = breaks, alpha = alpha, sig_cap = sig_cap, sig_beta = sig_beta)

    run_res <- build_hubs_from_EP(EP_bs$loops, k_min=k_min, resolution=resolution, seed=seed + i, quantile_cutoff=quantile_cutoff, method=method, use_leiden = 'on')

    # 关键点：手动清理内存，防止随着循环进行越来越慢
    rm(EP_bs)

    if (i %% 5 == 0) gc() 

    if (is.null(run_res)) return(NULL)

    # 提取有效信息
    valid_hubs <- run_res$observed_hubs %>% 
      filter(!is.null(enhancers), lengths(enhancers) > 0)
    
    return(split(valid_hubs$enhancers, valid_hubs$promoter_id))
  })


  valid_runs <- compact(stability_runs)
  n_valid <- length(valid_runs)
  
  # 3. 计算得分
  scored_results <- observed_hubs %>%
    mutate(
      stability_metrics = map2(promoter_id, enhancers, function(p_id, ref_set) {
        # 获取该 promoter 在所有有效迭代中的表现
        tau <- 0.5 #jaccard_threshold_by_size(length(ref_set))
        
        js <- map_dbl(valid_runs, function(m) {
          cand_sets <- m[[p_id]] # 获取该次迭代中该 promoter 下的所有候选 hubs
          if (is.null(cand_sets)) return(0)
          # 找到最像的那一个
          max(vapply(cand_sets, function(x) calculate_jaccard(x, ref_set), numeric(1)))
        })
        ov_best <- map_int(valid_runs, function(m) {
          cand_sets <- m[[p_id]]
          if (is.null(cand_sets) || length(cand_sets) == 0) return(0L)

          # 先找出“最像的那个 hub”（按 Jaccard 最大）
          jvals <- vapply(cand_sets, function(x) calculate_jaccard(x, ref_set), numeric(1))
          k <- which.max(jvals)

          # 计算该最佳匹配 hub 的 overlap 数
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



# ------------------- Step 8: Build null model for hub p-value calculation -------------------
build.null.pvalue <- function(hub_tbl, hub_tbl_ep, hub_tbl_all_ep, method = "log_minmax", weight_method = "bin_log_ratio_sig", B = 1000, quantile_cutoff = 0, stat = "density") {
  print("Building null model for hub p-value calculation...")
  # B <- 100
  # method <- "log_minmax"      
  # quantile_cutoff <- 0        # 与你 fast_comembership 设置一致
  # stat <- "density"           # 推荐；若你坚持“总强度”，用 "sum"

  ## =========================================================
  ## 1) 距离分箱（全局 pool 的 bin）
  ## =========================================================
  # breaks <- c(0, 1e4, 5e4, 1e5, 2.5e5, 5e5, 1e6, 1.5e6, 2e6, Inf)
  breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6, Inf)

  ## 注意：你的 EP 表里同时出现过 D 和 distance，这里统一使用 D。
  ## 如果你的列名叫 distance，就把 D 改成 distance。

  ## =========================================================
  ## 2) 选择用哪个“节点权重”做 null 与 hub score
  ## =========================================================
  ## 你说 weight_raw 是用于算 comembership 的——这里就用 weight_raw。
  ## 但你的 EP 明细表里列名可能是 weight_raw.x 或 weight_raw.y（join 后）。
  ## 下面做一个鲁棒选择：优先 weight_raw.x，其次 weight_raw，再其次 weight_raw.y
  pick_weight_col <- function(df) {
    if ("weight_raw.x" %in% names(df)) return("weight_raw.x")
    if ("weight_raw"   %in% names(df)) return("weight_raw")
    if ("weight_raw.y" %in% names(df)) return("weight_raw.y")
    stop("Cannot find weight_raw column in the object (tried weight_raw.x / weight_raw / weight_raw.y).")
  }
  wcol_all <- pick_weight_col(hub_tbl_all_ep)
  wcol <- pick_weight_col(hub_tbl_ep)

  ## =========================================================
  ## 3) 构建全局 distance-stratified 权重池（节点层）
  ## =========================================================
  global_bins <- hub_tbl_all_ep %>%
    filter(pair_type == "PE") %>%                 # 只用 EP
    filter(is.finite(D), D >= 0) %>%
    mutate(dist_bin = cut(D, breaks = breaks, include.lowest = TRUE, right = TRUE)) %>%
    group_by(dist_bin) %>%
    summarise(weight_pool = list(.data[[wcol_all]]), .groups = "drop")
  ## 一个安全检查：每个 bin 至少要有一些值
  if (any(lengths(global_bins$weight_pool) < 50)) {
    message("Warning: some distance bins have small pool size (<50). Consider merging bins or using fewer breaks.")
  }

  ## =========================================================
  ## 4) 为每个 hub 准备 dist_bins（关键：从 EP 明细表提取每个 enhancer 的 D）
  ## =========================================================
  ## 我们从 hub_tbl_ep 中提取：hub_id + enhancer_id + D + weight_raw
  ## 去重后按 hub 汇总成 list（与 hub_tbl 的 enhancers 对齐）
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

  ## 将 hub_nodes 合并回 hub_tbl
  hub_tbl2 <- hub_tbl %>%
    select(promoter_id, hub_id, enhancers, weight_raw, size, cm_ok, weighted_density, internal_score) %>%
    left_join(hub_nodes, by = c("promoter_id"))

  ## =========================================================
  ## 5.5) 关键：避免闭包捕获巨大环境（并行 globals 爆炸）
  ## =========================================================
  # environment(hub_score_from_weights) <- .GlobalEnv
  # environment(calculate_global_p)     <- .GlobalEnv

  # assign("hub_score_from_weights", hub_score_from_weights, envir = .GlobalEnv)
  # assign("calculate_global_p",     calculate_global_p,     envir = .GlobalEnv)

  ## =========================================================
  ## 6) 计算全局 p-value（建议 stat 用 density，与 weighted_density 同类）
  ## =========================================================
  ## Step 1: flags
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

## =========================================================
## 5) hub score：用你已有的 fast_comembership() 逻辑计算
## =========================================================
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

  s <- as.numeric(sum(cm) / 2)  # 上三角总和（无向边）
  if (stat == "sum") return(s)
  return(s / (n * (n - 1) / 2)) # density
}

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

  # observed
  obs_score <- hub_score_from_weights(
    observed_weights,
    method = method,
    quantile_cutoff = quantile_cutoff,
    stat = stat
  )

  # 为每个节点/bin 找到对应 pool（注意 dist_bins 可能是 factor）
  key <- as.character(dist_bins)

  if (null_mode == "distance_matched") {
    # strict: one pool per node (length == n)
    pools <- pool_map[key]
    if (any(vapply(pools, is.null, logical(1)))) {
      return(list(p = NA_real_, obs = obs_score, null_median = NA_real_))
    }

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
    # relaxed but still distance-aware: match the hub's distance-bin histogram
    bin_counts <- table(key)
    bins <- names(bin_counts)

    pools_by_bin <- pool_map[bins]
    if (any(vapply(pools_by_bin, is.null, logical(1)))) {
      return(list(p = NA_real_, obs = obs_score, null_median = NA_real_))
    }
    counts <- as.integer(bin_counts)

    null_scores <- replicate(B, {
      # draw k samples from each bin pool, then concatenate
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
        global_pool <- unlist(pool_map[-(1:2)], use.names = FALSE)
        
        # 确定当前 Hub 需要抽取的节点数量 (n)
        n_nodes <- length(observed_weights)
        
        null_scores <- replicate(B, {
          # 直接从全基因组大池子里随机抽 n 个权重
          w_star <- sample(global_pool, size = n_nodes, replace = TRUE)
          
          hub_score_from_weights(
            w_star,
            method = method,
            quantile_cutoff = quantile_cutoff,
            stat = stat
          )
        })}

  # smoothed p-value to avoid 0
  p_val <- (1 + sum(null_scores >= obs_score, na.rm = TRUE)) / (1 + length(null_scores))

  list(
    p = p_val,
    obs = obs_score,
    null_median = median(null_scores, na.rm = TRUE)
  )
}



## calculate global p-values for each hub using parallelization
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

# ------------------- Step 10: End-to-end workflow -------------------
build.run.preprocess <- function(loop_file, loop_file_all, outdir, tss_input, filename, promoter_window=0) {
  prep_hichip_all <- preprocess_hichip(loop_file = loop_file, tss_input = tss_input, mode = "all", promoter_window = promoter_window)
  prep_hichip <- preprocess_hichip(loop_file = loop_file, tss_input = tss_input, mode = "promoter", promoter_window = promoter_window)
  all_ep_data = preprocess_hichip(loop_file = loop_file_all, tss_input = tss_input, mode = "promoter", promoter_window = promoter_window)

  save(prep_hichip_all, prep_hichip, all_ep_data, file = file.path(outdir, paste("multiple_result",filename,"hub","all.preprocess.RData",sep=".")))
  print("Preprocessing done.")
}

build.run.all <- function(outdir, method, weight_method, filename, promoter_window=0) {
  load(file.path(outdir, paste("multiple_result",filename,"hub","all.preprocess.RData",sep=".")))
  observed_hubs_per_promoter <- PEhub_full(raw_data = prep_hichip$loops, weight_mode = weight_method, dist_col = c("D"), q_col = "qvalue", breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6), alpha = 0.8, sig_cap = 10, sig_beta = 0.5, n_perm = 10, resolution = 1.0, k_min = 3, seed = 42, quantile_cutoff = 0.2, method = method)
  observed_hubs_per_promoter_sub <- observed_hubs_per_promoter$loops[!is.na(observed_hubs_per_promoter$loops$hub_index) & observed_hubs_per_promoter$loops$hub_index=="hub1",]
  observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter$hubs[!is.na(observed_hubs_per_promoter$hubs$hub_index) & observed_hubs_per_promoter$hubs$hub_index=="hub1",]
  observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_sub_hub[order(observed_hubs_per_promoter_sub_hub$hubid),]

  #pvalue prepare
  all_ep_data_weight <- compute_weights(all_ep_data$loops, weight_mode = weight_method, dist_col = c("D"), q_col = "qvalue", breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6), alpha = 0.8, sig_cap = 10, sig_beta = 0.5)
  observed_hubs_per_promoter_sub_observed <- observed_hubs_per_promoter$observed_hubs[observed_hubs_per_promoter$observed_hubs$promoter_id %in% observed_hubs_per_promoter_sub$promoter_id,]
  observed_hubs_per_promoter_pvalue_prepare <- build.null.pvalue(hub_tbl=observed_hubs_per_promoter_sub_observed, hub_tbl_ep=observed_hubs_per_promoter_sub, hub_tbl_all_ep=all_ep_data_weight$loops, method = method, B = 1000, quantile_cutoff = 0.2, stat = "density")

  save(observed_hubs_per_promoter, observed_hubs_per_promoter_sub, observed_hubs_per_promoter_sub_hub, observed_hubs_per_promoter_pvalue_prepare, file = file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"preprocess.RData",sep=".")))
  print("All runs done.")
}

build.cutoff <- function(outdir, method, weight_method, filename) {
  print(paste("Starting post-cutoff:", filename))
  load(file.path(outdir, paste("multiple_result",filename,"hub","all.preprocess.RData",sep=".")))
  load(file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"preprocess.RData",sep=".")))
  
  # #stability 
  observed_hubs_per_promoter_stability_all <- calculate_hub_stability(raw_data = prep_hichip$loops, observed_hubs = observed_hubs_per_promoter$observed_hubs[!is.null(observed_hubs_per_promoter$observed_hubs$enhancers) & observed_hubs_per_promoter$observed_hubs$hub_index=="hub1",], weight_mode = weight_method, dist_col = c("D"), q_col = "qvalue", breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6), alpha = 0.8, sig_cap = 10, sig_beta = 0.5, B = 10, resolution = 1.0,  k_min = 3, quantile_cutoff = 0.2, method = method)
  observed_hubs_per_promoter_stability <- observed_hubs_per_promoter_stability_all %>% select(promoter_id, hub_id, jaccard_stability_score, reproducibility_rate, existence_rate) # %>% filter(!is.na(hub_p_adj_global), hub_p_adj_global <= 0.05)
  observed_hubs_per_promoter_sub <- observed_hubs_per_promoter_sub %>% left_join( observed_hubs_per_promoter_stability, by = c("promoter_id", "hub_id") )
  observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_sub_hub %>% left_join( observed_hubs_per_promoter_stability, by = c("promoter_id", "hub_id") )
  #save(prep_hichip_all, observed_hubs_per_promoter, observed_hubs_per_promoter_sub, observed_hubs_per_promoter_sub_hub, observed_hubs_per_promoter_stability, file = file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"stability.RData",sep=".")))
  #observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_stability %>%  filter(reproducibility_rate >= 0.5 & existence_rate >= 0.5)

  #pvalue
  observed_hubs_per_promoter_pvalue_all <- build.null.pvalue.calculate(hub_tbl3=observed_hubs_per_promoter_pvalue_prepare$hubs, global_bins=observed_hubs_per_promoter_pvalue_prepare$global_bins, method = method, B = 1000, quantile_cutoff = 0.2, stat = "sum", null_mode = "hist_matched")
  observed_hubs_per_promoter_pvalue <- observed_hubs_per_promoter_pvalue_all %>% select(promoter_id, hub_id, hub_p_value_global, hub_p_adj_global, hub_score_obs_global, hub_score_null_median, OE_ratio_global) # %>% filter(!is.na(hub_p_adj_global), hub_p_adj_global <= 0.05)
  observed_hubs_per_promoter_sub <- observed_hubs_per_promoter_sub %>% left_join( observed_hubs_per_promoter_pvalue, by = c("promoter_id", "hub_id") )
  observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_sub_hub %>% left_join( observed_hubs_per_promoter_pvalue, by = c("promoter_id", "hub_id") )
  save(prep_hichip_all, prep_hichip, observed_hubs_per_promoter, observed_hubs_per_promoter_sub, observed_hubs_per_promoter_sub_hub, observed_hubs_per_promoter_stability, observed_hubs_per_promoter_pvalue, file = file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"RData",sep=".")))
}

build.postprocess <- function(outdir, method, weight_method, filename) {
  print(paste("Starting post-processing:", filename))
  load(file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"RData",sep=".")))
  
  #export significant hubs
  observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_sub_hub %>% filter(!is.na(hub_p_adj_global), hub_p_adj_global <= 0.05, reproducibility_rate >= 0.5) # | OE_ratio_global > 1.5
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
##For real data
##########################
##For GM12878 test data
build.run.preprocess(loop_file = "/home/syidan/syidan/Data/Processed/HiCHIP_GSE_test/significant_interactions/Hicdcplus.Human_GM12878_unknown_WT_unknown_HiCHIP_standard_mergedSRR.significant_interactions.bedpe", 
              loop_file_all = "~/syidan/Data/Processed/HiCHIP_GSE_test/hicdcplus/Human_GM12878_unknown_WT_unknown_HiCHIP_standard_mergedSRR/Hicdcplus.Human_GM12878_unknown_WT_unknown_HiCHIP_standard_mergedSRR.all_interactions.bedpe.gz", 
              outdir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer", tss_input = promoters_annotation, filename = "exampleGM12878")

for (method in c("log_minmax")) { #, "log_maxnorm", "log_zscore"
  # for (weight_method in c("bin_log_ratio_sig")) { #"bin_percentile_plus_sig", "bin_log_ratio_sig", "bin_diff_global", "bin_diff_binmax"
  for (weight_method in c("distance_only", "count_only", "sig_only", "count_sig", "zscore_residual", "count_sig_plus_dist_linear", "bin_percentile_plus_sig", "bin_diff_global", "bin_diff_binmax", "bin_percentile")) {
      build.run.all(outdir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer", 
                    method = method, weight_method = weight_method, filename = "exampleGM12878")
  }
}

for (method in c("log_minmax")) { #, "log_maxnorm", "log_zscore"
  # for (weight_method in c("bin_log_ratio_sig")) { #"bin_percentile_plus_sig", "bin_log_ratio_sig", "bin_diff_global", "bin_diff_binmax"
  for (weight_method in c("distance_only", "count_only", "sig_only", "count_sig", "zscore_residual", "count_sig_plus_dist_linear", "bin_percentile_plus_sig", "bin_diff_global", "bin_diff_binmax", "bin_percentile")) {
    build.cutoff(outdir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer",  
                  method = method, weight_method = weight_method, filename = "exampleGM12878")
  }
}

for (method in c("log_minmax")) { #, "log_maxnorm", "log_zscore"
  # for (weight_method in c("bin_log_ratio_sig")) { #"bin_percentile_plus_sig", "bin_log_ratio_sig", "bin_diff_global", "bin_diff_binmax"
  for (weight_method in c("distance_only", "count_only", "sig_only", "count_sig", "zscore_residual", "count_sig_plus_dist_linear", "bin_percentile_plus_sig", "bin_diff_global", "bin_diff_binmax", "bin_percentile")) {
    build.postprocess(outdir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer",  
                  method = method, weight_method = weight_method, filename = "exampleGM12878")
  }
}



##For brain data
brain.dir= "/home/syidan/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged/HiCHIP/significant_interactions/"
for (samplename in c("Caudate1", "Caudate2", "Hippocampus1", "Hippocampus2", "MiddleFrontalGyrus1", "MiddleFrontalGyrus2", "ParietalLobe1", "ParietalLobe2", "SubstantiaNigra1", "SubstantiaNigra2", "SuperiorTemporalGyri1", "SuperiorTemporalGyri2")) {
    # build.run.preprocess(loop_file = file.path(brain.dir, paste0("Hicdcplus.Human_",samplename,"_unknown_WT_unknown_HiChIP_standard_mergedSRR.significant_interactions.bedpe")),
    #               loop_file_all = file.path(gsub("significant_interactions","hicdcplus",brain.dir), paste0("Human_",samplename,"_unknown_WT_unknown_HiChIP_standard_mergedSRR"), paste0("Hicdcplus.Human_",samplename,"_unknown_WT_unknown_HiChIP_standard_mergedSRR.all_interactions.bedpe.gz")), 
    #               outdir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer", tss_input = promoters_annotation, filename = samplename)

    build.run.all(outdir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer", 
                  method = "log_minmax", weight_method = "bin_log_ratio_sig", filename = samplename)
}

for (samplename in c("Caudate1", "Caudate2", "Hippocampus1", "Hippocampus2", "MiddleFrontalGyrus1", "MiddleFrontalGyrus2", "ParietalLobe1", "ParietalLobe2", "SubstantiaNigra1", "SubstantiaNigra2", "SuperiorTemporalGyri1", "SuperiorTemporalGyri2")) {
    build.cutoff(outdir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer", 
                  method = "log_minmax", weight_method = "bin_log_ratio_sig", filename = samplename)
}

for (samplename in c("Caudate1", "Caudate2", "Hippocampus1", "Hippocampus2", "MiddleFrontalGyrus1", "MiddleFrontalGyrus2", "ParietalLobe1", "ParietalLobe2", "SubstantiaNigra1", "SubstantiaNigra2", "SuperiorTemporalGyri1", "SuperiorTemporalGyri2")) {
    build.postprocess(outdir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer", 
                  method = "log_minmax", weight_method = "bin_log_ratio_sig", filename = samplename)
}




# Annotate loop bedpe file
# h3k27ac_peak = "~/syidan/Data/Processed/HiCHIP_GSE/ATAC/merged.snakePipes.out/MACS2/Human_GM12878_unknown_WT_unknown_ATAC_standard_mergedSRR.filtered.short.BAM_summits.bed"
# h3k27ac_peak = "/home/syidan/syidan/Data/Processed/HiCHIP_GSE_test/peaks/Human_GM12878_unknown_WT_unknown_HiCHIP_standard_mergedSRR_peaks.narrowPeak")

outdir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer"
method = "log_minmax"
weight_method = "bin_log_ratio_sig"
filename = "MiddleFrontalGyrus1"










##using multiple sessions to do bootstrap
  # stability_runs <- future_map(1:B, function(i) {
  #   # 模拟噪声：Poisson 采样
  #   # 注意：这里直接覆盖 counts 列，确保 build_hubs 内部的权重计算使用新值
  #   EP_bs <- raw_data %>% 
  #     mutate(counts = rpois(n(), lambda = pmax(0, subsample_frac * counts))) 
    
  #   EP_bs <- compute_weights(EP_bs,weight_mode = weight_mode, dist_col = dist_col, q_col = q_col, breaks = breaks, alpha = alpha, sig_cap = sig_cap, sig_beta = sig_beta)

  #   run_res <- build_hubs_from_EP(EP_bs$loops, k_min=k_min, resolution=resolution, seed=42, quantile_cutoff=quantile_cutoff, method=method, use_leiden = 'on')
    
  #   # 建立快速查找表：Promoter_id -> List of Enhancer Sets
  #   iter_hubs <- run_res$observed_hubs
  #   valid_hubs <- run_res$observed_hubs %>% 
  #     filter(!is.null(enhancers)) %>%
  #     filter(lengths(enhancers) > 0) # 确保 list 里的 vector 长度 > 0
  #   split(valid_hubs$enhancers, valid_hubs$promoter_id)
    
  #   }, .options = furrr_options(seed = TRUE,
  #   # 显式告诉 future 需要导出哪些函数和变量
  #   globals = c("raw_data", "compute_weights", "build_hubs_from_EP", 
  #               "weight_mode", "dist_col", "q_col", "breaks", "alpha", 
  #               "sig_cap", "sig_beta", "k_min", "resolution", "quantile_cutoff", 
  #               "method", "subsample_frac", "loops_convert", "make_peak_id",
  #               "run_leiden_for_promoter", "fast_comembership"),
  #   packages = c("dplyr", "purrr", "data.table", "Matrix", "igraph", "leiden", "tidyr")
  #   )
  # ) 

 


# # ------------------- Step 10: End-to-end workflow -------------------
# build.run.preprocess <- function(loop_file, loop_file_all, outdir, tss_input, filename, promoter_window=0) {
#   prep_hichip_all <- preprocess_hichip(loop_file = loop_file, tss_input = tss_input, mode = "all", promoter_window = promoter_window)
#   prep_hichip <- preprocess_hichip(loop_file = loop_file, tss_input = tss_input, mode = "promoter", promoter_window = promoter_window)
#   all_ep_data = preprocess_hichip(loop_file = loop_file_all, tss_input = tss_input, mode = "promoter", promoter_window = promoter_window)

#   save(prep_hichip_all, prep_hichip, all_ep_data, file = file.path(outdir, paste("multiple_result",filename,"hub","all.preprocess.RData",sep=".")))
#   print("Preprocessing done.")
# }

# build.run.all <- function(outdir, method, weight_method, filename, promoter_window=0) {
#   load(file.path(outdir, paste("multiple_result",filename,"hub","all.preprocess.RData",sep=".")))
#   observed_hubs_per_promoter <- PEhub_full(raw_data = prep_hichip$loops, weight_mode = weight_method, dist_col = c("D"), q_col = "qvalue", breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6), alpha = 0.8, sig_cap = 10, sig_beta = 0.5, n_perm = 10, resolution = 1.0, k_min = 2, seed = 42, quantile_cutoff = 0.2, method = method)
#   observed_hubs_per_promoter_sub <- observed_hubs_per_promoter$loops[!is.na(observed_hubs_per_promoter$loops$hub_index) & observed_hubs_per_promoter$loops$hub_index=="hub1",]
#   observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter$hubs[!is.na(observed_hubs_per_promoter$hubs$hub_index) & observed_hubs_per_promoter$hubs$hub_index=="hub1",]
#   observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_sub_hub[order(observed_hubs_per_promoter_sub_hub$hubid),]
#   #pvalue prepare
#   all_ep_data_weight <- compute_weights(all_ep_data$loops, weight_mode = weight_method, dist_col = c("D"), q_col = "qvalue", breaks = c(0,1e4,2.5e4,5e4,1e5,2.5e5,5e5,1e6,2e6), alpha = 0.8, sig_cap = 10, sig_beta = 0.5)
#   observed_hubs_per_promoter_sub_observed <- observed_hubs_per_promoter$observed_hubs[observed_hubs_per_promoter$observed_hubs$promoter_id %in% observed_hubs_per_promoter_sub$promoter_id,]
#   observed_hubs_per_promoter_pvalue_prepare <- build.null.pvalue(hub_tbl=observed_hubs_per_promoter_sub_observed, hub_tbl_ep=observed_hubs_per_promoter_sub, hub_tbl_all_ep=all_ep_data_weight$loops, method = method, B = 1000, quantile_cutoff = 0.2, stat = "density")

#   save(observed_hubs_per_promoter, observed_hubs_per_promoter_sub, observed_hubs_per_promoter_sub_hub, observed_hubs_per_promoter_pvalue_prepare, file = file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"preprocess.RData",sep=".")))
#   print("All runs done.")
# }

# build.postprocess <- function(outdir, method, weight_method, filename) {
#   load(file.path(outdir, paste("multiple_result",filename,"hub","all.preprocess.RData",sep=".")))
#   load(file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"preprocess.RData",sep=".")))
#   observed_hubs_per_promoter_pvalue_all <- build.null.pvalue.calculate(hub_tbl3=observed_hubs_per_promoter_pvalue_prepare$hubs, global_bins=observed_hubs_per_promoter_pvalue_prepare$global_bins, method = method, B = 1000, quantile_cutoff = 0.2, stat = "sum", null_mode = "global")
#   observed_hubs_per_promoter_pvalue <- observed_hubs_per_promoter_pvalue_all %>% select(promoter_id, hub_id, hub_p_value_global, hub_p_adj_global, hub_score_obs_global, hub_score_null_median, OE_ratio_global) # %>% filter(!is.na(hub_p_adj_global), hub_p_adj_global <= 0.05)
#   observed_hubs_per_promoter_sub <- observed_hubs_per_promoter_sub %>% left_join( observed_hubs_per_promoter_pvalue, by = c("promoter_id", "hub_id") )
#   observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_sub_hub %>% left_join( observed_hubs_per_promoter_pvalue, by = c("promoter_id", "hub_id") )
#   save(prep_hichip_all, observed_hubs_per_promoter, observed_hubs_per_promoter_sub, observed_hubs_per_promoter_sub_hub, observed_hubs_per_promoter_pvalue, file = file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"RData",sep=".")))

#   #export significant hubs
#   observed_hubs_per_promoter_sub_hub <- observed_hubs_per_promoter_sub_hub %>% filter(!is.na(hub_p_adj_global), hub_p_adj_global <= 0.05 | OE_ratio_global > 1.5)
#   print(paste0("Significant hubs number: ", nrow(observed_hubs_per_promoter_sub_hub)))
#   observed_hubs_per_promoter_sub_hub$name <- with(observed_hubs_per_promoter_sub_hub[order(observed_hubs_per_promoter_sub_hub$hubid),], { h <- gsub("_","", observed_hubs_per_promoter_sub_hub$hubid); paste0(observed_hubs_per_promoter_sub_hub$promoter_gene_name, "_", h, "_", ave(h, h, FUN = seq_along)) })
#   gr <- makeGRangesFromDataFrame(observed_hubs_per_promoter_sub_hub,keep.extra.columns = TRUE)
#   export(gr, con=file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"bed",sep=".")))
#   write.table(as.data.frame(observed_hubs_per_promoter_sub_hub)[,c(1:21,ncol(observed_hubs_per_promoter_sub_hub))], file=file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"txt",sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
#   colnames(observed_hubs_per_promoter_sub)[1] <- paste0("#",colnames(observed_hubs_per_promoter_sub)[1])
#   write.table(observed_hubs_per_promoter_sub[,1:12], file=file.path(outdir, paste("multiple_result",filename,"hub",method,weight_method,"bedpe",sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
  
#   #pairwise
#   observed_hubs_per_promoter_sub_hub_pairs <- prep_hichip$hubs[!prep_hichip$hubs$hubid %in% observed_hubs_per_promoter_sub_hub$hubid,]
#   observed_hubs_per_promoter_sub_pairs <- prep_hichip$loops[!prep_hichip$loops$hubid %in% observed_hubs_per_promoter_sub$hubid,]
#   gr_pairs <- makeGRangesFromDataFrame(observed_hubs_per_promoter_sub_hub_pairs,keep.extra.columns = TRUE)
#   gr_pairs$name <- gr_pairs$hubid
#   export(gr_pairs, con=file.path(outdir, paste("multiple_result",filename,"pairwise",method,weight_method,"bed",sep=".")))
#   write.table(as.data.frame(observed_hubs_per_promoter_sub_hub_pairs), file=file.path(outdir, paste("multiple_result",filename,"pairwise",method,weight_method,"txt",sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
#   colnames(observed_hubs_per_promoter_sub_pairs)[1] <- paste0("#",colnames(observed_hubs_per_promoter_sub_pairs)[1])
#   write.table(observed_hubs_per_promoter_sub_pairs[,1:12], file=file.path(outdir, paste("multiple_result",filename,"pairwise",method,weight_method,"bedpe",sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
#   print("Postprocessing done.")
# }
