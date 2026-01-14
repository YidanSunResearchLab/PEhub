###############################################################################
### LIBRARIES
###############################################################################

library(GenomicRanges)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(readr)
indir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer"



# ------------------------- step 1: Build replicable hubs for one sample with 2 reps -------------------------

# -----------------------
# 1) parse sample name
# -----------------------
parse_sample <- function(sample) {
  tibble(
    sample = sample,
    rep    = as.integer(str_extract(sample, "\\d+$")),
    region = str_remove(sample, "\\d+$")
  )
}

# -----------------------
# 2) loop -> promoter-hub summary per replicate
#    (1 promoter = 1 hub, you已确认)
# -----------------------
hub_summary <- function(loop_tbl, sample) {
  meta <- parse_sample(sample)

  loop_tbl %>%
    mutate(sample = sample) %>%
    left_join(meta, by = "sample") %>%
    group_by(region, rep, promoter_id) %>%
    summarise(
      enhancers = list(sort(unique(enhancer_id))),
      size      = length(unique(enhancer_id)),
      fdr       = first(hub_p_adj_global),
      stability = first(reproducibility_rate),
      promoter_gene_id   = first(promoter_gene_id),
      promoter_gene_name = first(promoter_gene_name),
      .groups = "drop"
    ) %>%
    mutate(sig = !is.na(fdr) & fdr <= 0.05 & size >= 3)
}

# -----------------------
# 3) promoter-level replicable hubs
# -----------------------
find_replicable_hubs <- function(
  loops_list,
  min_overlap = 2,
  jaccard_thresh = 0.1
) {
  hubs <- imap_dfr(loops_list, hub_summary)

  hubs_wide <- hubs %>%
    mutate(rep = paste0("rep", rep)) %>%
    pivot_wider(
      names_from  = rep,
      values_from = c(enhancers, size, stability, fdr, sig)
    )

  hubs_wide %>%
    mutate(
      overlap = map2_int(enhancers_rep1, enhancers_rep2,
                         ~ length(intersect(.x, .y))),
      jaccard = map2_dbl(enhancers_rep1, enhancers_rep2, function(a, b) {
        u <- length(union(a, b))
        if (u == 0) return(0)
        length(intersect(a, b)) / u
      }),
      sig_any = sig_rep1 | sig_rep2,
      replicable = sig_any &
        overlap >= min_overlap &
        jaccard >= jaccard_thresh &
        stability_rep1 >= 0.5 & stability_rep2 >= 0.5
    ) %>%
    filter(replicable)
}

# -----------------------
# 4) enhancer-level counts (rep 内先 sum)
# -----------------------
enhancer_counts <- function(loops_list) {
  imap_dfr(loops_list, function(df, sample) {
    meta <- parse_sample(sample)

    df %>%
      mutate(sample = sample) %>%
      left_join(meta, by = "sample") %>%
      select(region, rep, promoter_id, enhancer_id,
             counts, promoter_gene_id, promoter_gene_name) %>%
      group_by(region, rep, promoter_id, enhancer_id,
               promoter_gene_id, promoter_gene_name) %>%
      summarise(counts = sum(counts, na.rm = TRUE), .groups = "drop")
  })
}

# -----------------------
# 5) final intersect / union tables (only replicable hubs)
# -----------------------
make_final_tables <- function(loops_list) {

  # promoter-level replicable hubs
  rep_hubs <- find_replicable_hubs(loops_list)
  print(paste0("Final replicable number: "))
  print(table(rep_hubs$region))

  # enhancer-level counts
  enh <- enhancer_counts(loops_list)

  enh_filt <- enh %>%
    semi_join(rep_hubs, by = c("region", "promoter_id"))

  enh_wide <- enh_filt %>%
    mutate(rep = paste0("rep", rep)) %>%
    pivot_wider(names_from = rep, values_from = counts)

  # INTERSECT
  enh_intersect <- enh_wide %>%
    filter(!is.na(rep1), !is.na(rep2)) %>%
    transmute(
      region, promoter_id, promoter_gene_id, promoter_gene_name,
      enhancer_id,
      counts_rep1 = rep1,
      counts_rep2 = rep2,
      counts_sum  = rep1 + rep2
    )

  # UNION
  enh_union <- enh_wide %>%
    transmute(
      region, promoter_id, promoter_gene_id, promoter_gene_name,
      enhancer_id,
      counts_rep1 = coalesce(rep1, 0),
      counts_rep2 = coalesce(rep2, 0),
      counts_sum  = counts_rep1 + counts_rep2
    )

  # filter loops_list by replicable hubs (region + promoter_id)
  # Build a lookup: region -> vector(promoter_id)
  keep_map <- rep_hubs %>%
    distinct(region, promoter_id) %>%
    group_by(region) %>%
    summarise(keep_promoters = list(unique(promoter_id)), .groups = "drop")

  # helper: infer region from sample name via your existing parse_sample()
  get_region_from_sample <- function(sample) {
    parse_sample(sample)$region[1]
  }

  new_loops_list <- purrr::imap(loops_list, function(loop_tbl, sample) {
    reg <- get_region_from_sample(sample)
    keep_p <- keep_map$keep_promoters[match(reg, keep_map$region)]
    if (is.na(keep_p)) {
      # no replicable hubs for this region -> return empty with same cols
      return(loop_tbl[0, , drop = FALSE])
    }
    loop_tbl %>% filter(promoter_id %in% keep_p[[1]])
  })

  list(
    replicable_hubs    = rep_hubs,
    enhancer_intersect = enh_intersect,
    enhancer_union     = enh_union,
    new_loops_list     = new_loops_list
  )
  }

###############################################################################
# ------------------------- step 2: Export into txt/bed/bedpe files -------------------------

# peak_chr10_103854_103860 -> chr10, start0, end0  (BED 0-based)
to_bed <- function(id, bin = 1000L) {
  m <- str_match(id, "^peak_(chr[^_]+)_(\\d+)_(\\d+)$")
  if (any(is.na(m))) stop("Bad peak id (example expected: peak_chr10_103854_103860).")
  s <- as.integer(m[,3]); e <- as.integer(m[,4])
  tibble(chr = m[,2],
         start = (s + 1L) * bin - 1L,
         end   = e * bin)
}

export_by_region <- function(df, outdir = ".", prefix = "replicable", bin = 1000L) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  df %>%
    dplyr::group_split(region, .keep = TRUE) %>%
    purrr::walk(function(dfr) {

      reg <- dfr$region[1]
      reg_safe <- stringr::str_replace_all(reg, "[^A-Za-z0-9_.-]+", "_")  # 文件名安全

      base <- file.path(outdir, paste("multiple_result", reg_safe, prefix, sep="."))

      # 1) txt
      readr::write_tsv(dfr, paste0(base, ".txt"))

      # promoter/enhancer coords
      prom <- dplyr::bind_cols(dfr, to_bed(dfr$promoter_id, bin)) %>%
        dplyr::rename(p_chr = chr, p_start = start, p_end = end)

      pe <- dplyr::bind_cols(prom, to_bed(dfr$enhancer_id, bin)) %>%
        dplyr::rename(e_chr = chr, e_start = start, e_end = end)

      # 2) bedpe (10列标准 bedpe)
      bedpe <- pe %>%
        dplyr::transmute(
          chrom1 = p_chr, start1 = p_start, end1 = p_end,
          chrom2 = e_chr, start2 = e_start, end2 = e_end,
          name   = paste(promoter_id, promoter_gene_id, promoter_gene_name, sep="|"),
          score  = counts_sum,
          strand1 = ".", strand2 = "."
        )
      readr::write_tsv(bedpe, paste0(base, ".bedpe"), col_names = FALSE)

      # 3) ONE combined BED (promoter + enhancer), sorted
      bed_all <- dplyr::bind_rows(
        # promoter: will repeat across rows in pe, so dedup within each promoter_id
        pe %>%
          dplyr::transmute(
            promoter_id,
            chrom = p_chr, start = p_start, end = p_end,
            promoter_gene_id, promoter_gene_name,
            feature = "promoter",
            enhancer_id = NA_character_
          ) %>%
          dplyr::distinct(promoter_id, chrom, start, end, .keep_all = TRUE),

        # enhancer: dedup ONLY within each promoter_id (do NOT dedup globally)
        pe %>%
          dplyr::transmute(
            promoter_id,
            chrom = e_chr, start = e_start, end = e_end,
            promoter_gene_id, promoter_gene_name,
            feature = "enhancer",
            enhancer_id
          ) %>%
          dplyr::distinct(promoter_id, enhancer_id, .keep_all = TRUE)
      ) %>%
        dplyr::group_by(promoter_id) %>%
        dplyr::arrange(
          dplyr::desc(feature == "promoter"),  # promoter first within hub
          start, end, enhancer_id,
          .by_group = TRUE
        ) %>%
        dplyr::mutate(idx = dplyr::row_number()) %>%
        dplyr::ungroup() %>%
        dplyr::transmute(
          chrom, start, end,
          name   = paste(promoter_id, idx, sep = "_"),
          score  = 0,
          strand = "."
        )

      readr::write_tsv(bed_all, paste0(base, ".bed"), col_names = FALSE)

      message("Wrote: ", reg, " -> ", base, ".[txt|bedpe|bed]")
    })
}

loops_convert <- function(dt) {

  data.table::setDT(dt)  # 不会深拷贝到每一列

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


save_loops_by_region <- function(new_loops_list, outdir = ".", prefix = "new_loops") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # 用你现有的 parse_sample() 解析 region（推荐，最稳）
  get_region <- function(sample) {
    parse_sample(sample)$region[1]
  }

  # 给每个元素加 region 信息
  meta <- tibble::tibble(
    sample = names(new_loops_list),
    region = vapply(names(new_loops_list), get_region, character(1))
  )

  # 按 region 拆分 sample 名
  split_samples <- split(meta$sample, meta$region)

  # 每个 region 保存一个 RData
  for (reg in names(split_samples)) {
    reg_safe <- stringr::str_replace_all(reg, "[^A-Za-z0-9_.-]+", "_")
    observed_hubs_per_promoter_sub <- new_loops_list[ split_samples[[reg]] ]
    observed_hubs_per_promoter_sub <- purrr::imap_dfr(
      observed_hubs_per_promoter_sub,
      ~ dplyr::mutate(.x, sample = .y)
    )
    observed_hubs_per_promoter_sub_hub <- loops_convert(observed_hubs_per_promoter_sub)
    # 文件名
    base <- file.path(outdir, paste("multiple_result", reg_safe, prefix, "RData",sep="."))

    # 保存对象名 observed_hubs_per_promoter_sub（你也可改成别的）
    save(observed_hubs_per_promoter_sub, observed_hubs_per_promoter_sub_hub, file = base)
    message("Saved: ", reg, " -> ", base)
  }
}


###############################################################################
##Pairwise file read and merge
read_one_pairwise <- function(f) {
  df <- readr::read_tsv(f, show_col_types = FALSE)

  # 只保留 promoter-enhancer（去掉 PP/EE 等）
  if (!"pair_type" %in% names(df)) stop("Missing column pair_type in: ", f)

  need <- c("promoter_id", "enhancer_id", "counts", "promoter_gene_id", "promoter_gene_name")
  miss <- setdiff(need, names(df))
  if (length(miss) > 0) stop("Missing columns in ", f, ": ", paste(miss, collapse = ","))

  # 从文件名解析 region + rep（假设 ...<region><1/2>... 结尾，如 SubstantiaNigra2）
  bn <- basename(f)
  m <- stringr::str_match(bn, "^multiple_result\\.(.+?)\\.pairwise\\..*\\.txt$")
  if (any(is.na(m))) stop("Bad filename format: ", bn)

  sample_tag <- m[, 2]                       # e.g. "SubstantiaNigra2"
  m2 <- stringr::str_match(sample_tag, "^(.*?)([12])$")
  if (any(is.na(m2))) stop("Cannot parse region/rep from: ", sample_tag)

  region <- m2[, 2]
  rep    <- as.integer(m2[, 3])

  df %>%
    dplyr::filter(pair_type == "PE") %>%     # 只要 PE
    dplyr::transmute(
      region = region,
      rep    = rep,
      promoter_id,
      enhancer_id,
      counts = as.integer(counts),
      promoter_gene_id,
      promoter_gene_name
    )
}

make_pairwise_union_intersect <- function(pairwise_df) {

  pw <- pairwise_df %>%
    dplyr::group_by(region, rep, promoter_id, enhancer_id, promoter_gene_id, promoter_gene_name) %>%
    dplyr::summarise(counts = sum(counts, na.rm = TRUE), .groups = "drop")

  pw_wide <- pw %>%
    dplyr::mutate(rep = paste0("rep", rep)) %>%
    tidyr::pivot_wider(names_from = rep, values_from = counts)

  pw_intersect <- pw_wide %>%
    dplyr::filter(!is.na(rep1), !is.na(rep2)) %>%
    dplyr::transmute(
      region, promoter_id, promoter_gene_id, promoter_gene_name,
      enhancer_id,
      counts_rep1 = rep1,
      counts_rep2 = rep2,
      counts_sum  = rep1 + rep2
    )

  pw_union <- pw_wide %>%
    dplyr::transmute(
      region, promoter_id, promoter_gene_id, promoter_gene_name,
      enhancer_id,
      counts_rep1 = dplyr::coalesce(rep1, 0L),
      counts_rep2 = dplyr::coalesce(rep2, 0L),
      counts_sum  = counts_rep1 + counts_rep2
    )

  list(pairwise_intersect = pw_intersect,
       pairwise_union     = pw_union)
}


##################################
## Read in each replicate hub results
hubs_list <- list()
loops_list <- list()
for (samplename in c("Caudate1", "Caudate2", "Hippocampus1", "Hippocampus2", "MiddleFrontalGyrus1", "MiddleFrontalGyrus2", "ParietalLobe1", "ParietalLobe2", "SubstantiaNigra1", "SubstantiaNigra2", "SuperiorTemporalGyri1", "SuperiorTemporalGyri2")) {
  load(file.path(indir, paste("multiple_result",samplename,"hub","log_minmax","bin_log_ratio_sig","RData",sep=".")))
  hubs_list[[samplename]] <- observed_hubs_per_promoter_sub_hub  #%>% filter(!is.na(hub_p_value_global), hub_p_value_global <= 0.05, reproducibility_rate >= 0.5) # | OE_ratio_global > 1.5
  loops_list[[samplename]] <- observed_hubs_per_promoter_sub  #%>% filter(!is.na(hub_p_value_global), hub_p_value_global <= 0.05, reproducibility_rate >= 0.5) # | OE_ratio_global > 1.5
}

print(paste0("hub number in each sample originally."))
print(sapply(hubs_list,function(hublist){length(unique(hublist$hubid))}))
save(hubs_list, loops_list, file = file.path(indir, "brain.original_hubs.all.RData"))


# ##################################
# ## merge replicates and select replicable hubs
load(file.path(indir, "brain.original_hubs.all.RData"))
replicable_hubs <- make_final_tables(loops_list)
save(replicable_hubs, file = file.path(indir, "brain.replicable_hubs.all.RData"))


##################################
# Generate txt, bed, bedpe files for replicable hubs in each sample
load(file.path(indir, "brain.replicable_hubs.all.RData"))
export_by_region(df=replicable_hubs$enhancer_intersect, outdir=indir, prefix="replicable_hubs_bin_intersect", bin = 1000L)
export_by_region(df=replicable_hubs$enhancer_union,     outdir=indir, prefix="replicable_hubs_bin_union",     bin = 1000L)
save_loops_by_region(new_loops_list=replicable_hubs$new_loops_list, outdir=indir, prefix="replicable_hubs_loops")


##################################
# pairwise replicate merge
pairwise_list = list()
for (samplename in c("Caudate1", "Caudate2", "Hippocampus1", "Hippocampus2", "MiddleFrontalGyrus1", "MiddleFrontalGyrus2", "ParietalLobe1", "ParietalLobe2", "SubstantiaNigra1", "SubstantiaNigra2", "SuperiorTemporalGyri1", "SuperiorTemporalGyri2")) {
  pairwise_list[[samplename]] <- read_one_pairwise(file.path(indir, paste("multiple_result",samplename,"pairwise","log_minmax","bin_log_ratio_sig","txt",sep=".")))
}

pairwise_hubs <- make_pairwise_union_intersect(dplyr::bind_rows(pairwise_list))
save(pairwise_hubs, file = file.path(indir, "brain.pairwise_hubs.all.RData"))

load(file.path(indir, "brain.pairwise_hubs.all.RData"))
export_by_region(df = pairwise_hubs$pairwise_intersect, outdir = indir,prefix = "pairwise_bin_intersect",bin = 1000L)
export_by_region(df = pairwise_hubs$pairwise_union, outdir = indir, prefix = "pairwise_bin_union", bin = 1000L)
