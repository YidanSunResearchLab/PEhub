###############################################################################
### LIBRARIES
###############################################################################

library(GenomicRanges)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(rtracklayer)
indir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer"
regions <- c("SubstantiaNigra","MiddleFrontalGyrus","Caudate","Hippocampus","ParietalLobe","SuperiorTemporalGyri")



# ------------------------- step 2: Build complete loss table -------------------------
build.complete.loss.table <- function(
  loops_list,
  outdir = ".",
  replicates = c("intersect", "union"),
  min_hub_size = 2,
  verbose = TRUE
) {
  stopifnot(is.list(loops_list), length(loops_list) > 0)
  replicates <- match.arg(replicates)

  # sample -> region: "Caudate1" -> "Caudate"
  sample2region <- function(x) sub("[0-9]+$", "", x)

  promoter_summary_one <- function(loops_df, sample_id, min_hub_size = 2) {
    loops_df %>%
      filter(pair_type == "PE") %>%
      group_by(promoter_id, promoter_gene_id, promoter_gene_name) %>%
      summarise(
        sample = sample_id,
        hub_size = n_distinct(enhancer_id),
        hub_strength = sum(counts, na.rm = TRUE),
        has_hub = hub_size >= min_hub_size,
        .groups = "drop"
      )
  }

  promoter_by_sample <- imap_dfr(
    loops_list,
    ~ promoter_summary_one(.x, .y, min_hub_size = min_hub_size)
  ) %>%
    mutate(region = sample2region(sample))

  promoter_by_region <- promoter_by_sample %>%
    group_by(region, promoter_id, promoter_gene_id, promoter_gene_name) %>%
    summarise(
      n_rep = n(),
      has_hub_any  = any(has_hub),
      has_hub_both = all(has_hub) & n_rep >= 2,
      hub_size_max = max(hub_size, na.rm = TRUE),
      hub_strength_sum = sum(hub_strength, na.rm = TRUE),
      .groups = "drop"
    )

  # replicate rule filter + choose correct presence column
  if (replicates == "intersect") {
    promoter_by_region_present <- promoter_by_region %>% filter(has_hub_both)
    presence_value_col <- "has_hub_both"
  } else {
    promoter_by_region_present <- promoter_by_region %>% filter(has_hub_any)
    presence_value_col <- "has_hub_any"
  }

  # Presence matrix (0/1)
  presence_mat <- promoter_by_region_present %>%
    select(promoter_id, promoter_gene_id, promoter_gene_name, region, !!sym(presence_value_col)) %>%
    rename(present = !!sym(presence_value_col)) %>%
    pivot_wider(
      names_from = region,
      values_from = present,
      values_fill = FALSE
    ) %>%
    mutate(across(-c(promoter_id, promoter_gene_id, promoter_gene_name), as.integer))

  region_cols <- setdiff(names(presence_mat), c("promoter_id", "promoter_gene_id", "promoter_gene_name"))
  mat <- as.matrix(presence_mat[, region_cols, drop = FALSE])

  complete_loss_table <- presence_mat %>%
    mutate(
      n_regions_present = rowSums(mat),
      present_regions = apply(mat, 1, \(r) paste(region_cols[r == 1], collapse = ",")),
      absent_regions  = apply(mat, 1, \(r) paste(region_cols[r == 0], collapse = ","))
    )

  # --------- Verbose summary prints ---------
  if (verbose) {
    # How many promoters are present under union vs intersect (per region and overall)
    union_counts <- promoter_by_region %>%
      summarise(n_union = sum(has_hub_any), .by = region)

    intersect_counts <- promoter_by_region %>%
      summarise(n_intersect = sum(has_hub_both), .by = region)

    rep_stats <- left_join(union_counts, intersect_counts, by = "region") %>%
      mutate(
        kept_frac = ifelse(n_union > 0, n_intersect / n_union, NA_real_),
        dropped_frac = ifelse(n_union > 0, 1 - kept_frac, NA_real_)
      ) %>%
      arrange(desc(dropped_frac))

    cat("\n[Replicate consistency summary]\n")
    print(rep_stats, n = Inf)

    total_union <- sum(union_counts$n_union, na.rm = TRUE)
    total_intersect <- sum(intersect_counts$n_intersect, na.rm = TRUE)

    cat("\n[Overall replicate filtering effect]\n")
    cat("Total promoter-region calls passing UNION (any rep): ", total_union, "\n", sep = "")
    cat("Total promoter-region calls passing INTERSECT (both reps): ", total_intersect, "\n", sep = "")
    if (total_union > 0) {
      cat("INTERSECT relative to UNION: ",
          round(100 * (total_intersect / total_union), 2), "%\n", sep = "")
    }

    cat("\n[Distribution of n_regions_present in the final table]\n")
    tab <- table(complete_loss_table$n_regions_present)
    print(tab)

    cat("\n[Interpretation]\n")
    cat("- n_regions_present = 1: promoter forms a hub in exactly one region (region-specific / 'complete loss' in others).\n")
    cat("- n_regions_present = 6: promoter forms hubs in all regions (shared/core hubs).\n")
    cat("- Intermediate values (2–5): partially shared; good candidates for rewiring / gain-loss analysis.\n")
    cat("- If INTERSECT drops a large fraction, many 'loss' calls may reflect replicate inconsistency or power differences.\n")
  }

  complete_loss_table_list = list(
    complete_loss_table = complete_loss_table,
    presence_mat = presence_mat,
    promoter_by_region = promoter_by_region,
    promoter_by_region_present = promoter_by_region_present
  )
  save(complete_loss_table_list, file = file.path(outdir, paste("brain.complete_loss_table",replicates,"RData",sep=".")))
  write.table(complete_loss_table, file = file.path(outdir, paste("brain.complete_loss_table",replicates,"txt",sep=".")),row.names=FALSE, quote=FALSE, sep="\t")
  # return(complete_loss_table_list)
}

#gene expression data compare between different brain regions for each specific gene set
build.gene.expression.out.plot = function(gene_expression, gene_sample, complete_loss_table, regions, target_region, expr_threshold = 1,outdir = ".", select_top_genes = FALSE) {
  print(target_region)
  ## =========================
  ## 0) 基础参数
  ## =========================
  # target_region <- "MiddleFrontalGyrus"
  # expr_threshold <- 1
  # regions=c("MiddleFrontalGyrus","SubstantiaNigra","Caudate","Hippocampus","ParietalLobe","SuperiorTemporalGyri")

  ## =========================
  ## 1) gene_expression: 宽表 -> 长表，并与 gene_sample 合并
  ##    - gene_expression: ensembl_gene_id, hgnc_symbol, + 很多 sample 列 (HA74, HC34, ...)
  ##    - gene_sample$name1: 对应这些 sample 列名 (HA3, HB3, ...)
  ## =========================
  gene_sample_clean <- gene_sample %>%
    filter(name3 %in% regions) %>%                    # 去掉 regions 里没有的
    mutate(
      name3 = factor(name3, levels = regions)          # 按指定顺序排序
    ) %>%
    arrange(name3) %>%  # 确保 sample1 按 regions 排序
    select(2:4)

  ge_long <- gene_expression %>%
    pivot_longer(
      cols = -c(ensembl_gene_id, hgnc_symbol),
      names_to = "sample",
      values_to = "expr"
    ) %>%
    mutate(
      ensembl_gene_id = as.character(ensembl_gene_id),
      sample = as.character(sample),
      expr = as.numeric(expr)
    ) %>%
    inner_join(
      gene_sample_clean %>%
        transmute(
          sample = as.character(name1),
          region = as.character(name3),
          name2 = as.character(name2)
        ),
      by = "sample"
    )

  ## =========================
  ## 2) complete_loss_table: promoter_gene_id 拆分成单基因并去版本号
  ##    例如 "ENSG...20" -> "ENSG..."; 多基因用逗号分隔
  ## =========================
  cl_genes <- complete_loss_table %>%
    mutate(row_id = row_number()) %>%
    transmute(
      row_id,
      promoter_id,
      present_regions,
      promoter_gene_id_raw = promoter_gene_id
    ) %>%
    mutate(promoter_gene_id_raw = str_replace_all(promoter_gene_id_raw, "\\s+", "")) %>%
    separate_rows(promoter_gene_id_raw, sep = ",") %>%
    mutate(
      ensembl_gene_id = str_replace(promoter_gene_id_raw, "\\..*$", "")  # 去掉版本号 .20
    ) %>%
    select(-promoter_gene_id_raw)

  ## =========================
  ## 3) 三表“合并大表”：complete_loss_table x gene_expression x gene_sample_clean
  ## =========================
  merged_all <- cl_genes %>%
    inner_join(ge_long, by = "ensembl_gene_id")
  # merged_all 现在至少包含：
  # row_id, promoter_id, present_regions, ensembl_gene_id, hgnc_symbol, sample, expr, region, GSM, name2

  if (target_region != "all") {
    print(paste0("Filtering for target region:", target_region, "\n"))
    ## =========================
    ## 4) 在合并层面做筛选：present_regions == "MiddleFrontalGyrus"
    ##    注意：你这里是“完全等于”，不会包含 "Caudate,MiddleFrontalGyrus" 那类
    ## =========================
    ## merged_all 已经是：cl_genes inner_join(ge_long) 的结果
    ## merged_mfg_rows 是：present_regions == target_region 的子集
    merged_mfg_rows <- merged_all %>%
      filter(present_regions == target_region)

    if(select_top_genes) {
      print("Selecting top genes based on expression in target region \n")
      ## =========================
      ## 5) 对每个 promoter：在 MFG RNA 样本中，选 sum(expr) 最大的 row_id
      ## =========================
      row_scores <- merged_mfg_rows %>%
        # filter(region == target_region) %>%                 # 只用 MFG 的 RNA 样本打分
        group_by(promoter_id, row_id) %>%
        summarise(total_expr_sum = sum(expr, na.rm = TRUE), .groups = "drop")

      best_rows_per_promoter <- row_scores %>%
        group_by(promoter_id) %>%
        slice_max(total_expr_sum, n = 1, with_ties = FALSE) %>%
        ungroup()

      ## 这些是“每个 promoter 选出来的最佳 row”，对应的基因集合
      selected_genes <- cl_genes %>%
        semi_join(best_rows_per_promoter, by = "row_id") %>%
        distinct(ensembl_gene_id) %>%
        pull(ensembl_gene_id)

      ## =========================
      ## 6) 过滤基因：去掉“在所有 sample 中 max(expr) < 1”的基因
      ##    更严格/更符合字面：用所有 RNA samples（所有 region）来算 max
      ## =========================
      genes_keep <- merged_all %>%                           # 注意这里用 merged_all 才是“所有sample”
        filter(ensembl_gene_id %in% selected_genes, !is.na(expr)) %>%
        group_by(ensembl_gene_id) %>%
        summarise(max_expr = max(expr,na.rm=TRUE), .groups = "drop") %>%
        filter(max_expr >= expr_threshold) %>%
        pull(ensembl_gene_id)

    } else {

      print("Not selecting top genes based on expression in target region \n")
      genes_keep <- merged_all %>%                           # 注意这里用 merged_all 才是“所有sample”
        group_by(ensembl_gene_id) %>%
        summarise(max_expr = max(expr,na.rm=TRUE), .groups = "drop") %>%
        filter(max_expr >= expr_threshold) %>%
        pull(ensembl_gene_id)

    }

    ## 如果你想把“所有 sample”限定为 MFG RNA 样本，把上面 merged_all 改成下面这行即可：
    ## genes_keep <- merged_all %>% filter(region == target_region, ...)

    ## =========================
    ## 7) 用筛出来的 genes_keep 直接在 merged_mfg_rows 里画图
    ##    每个 sample 一个 box，box 内是 genes_keep 在该 sample 的 expr 分布
    ## =========================
    plot_df <- merged_mfg_rows %>%
      filter(ensembl_gene_id %in% genes_keep)
  } else {
    print("No Filtering for target region \n")
    plot_df <- merged_all
  }
  
  sample_order <- gene_sample_clean %>%
    mutate(
      orig_order = row_number(),                     # 原始顺序
      prefix = str_extract(name1, "^[A-Z]+")         # HA / HB / HC / HD
    ) %>%
    mutate(
      prefix = factor(prefix, levels = c("HA", "HB", "HC", "HD"))
    ) %>%
    arrange(prefix, orig_order) %>%
    pull(name1)

  plot_df2 <- plot_df %>%
    mutate(
      group = str_extract(as.character(sample), "^[A-Z]+"),        # HA/HB/HC/HD
      group = factor(group, levels = c("HA","HB","HC","HD")),
      sample = factor(sample, levels = sample_order)
    )
  if (nrow(plot_df2) > 0) {
    print("Plotting selected genes expression \n")
    p <- ggplot(plot_df2, aes(x = sample, y = log1p(expr), fill = region, color=name2)) +
      geom_boxplot(outlier.alpha = 0.3) +
      facet_grid(. ~ group, scales = "free_x", space = "free_x") +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing.x = unit(0.8, "lines")
      ) +
      labs(
        title = target_region,
        x = "",
        y = "log1p(TPM)"
      )
    ggsave(plot=p, file = file.path(outdir, paste("brain.selected_genes_expression",target_region,"pdf",sep=".")), width = 12, height = 3)
  } else {
    print("No data to plot \n")
  }
}


build.circos.plot <- function(complete_loss_table, regions, outdir = ".") {
  # 1) Chord diagram
  library(circlize)
  # module_counts: present_regions -> n_promoters
  module_counts <- complete_loss_table %>%
    count(present_regions, name = "n_promoters") %>%
    mutate(
      region_list = str_split(present_regions, ","),
      region_list = map(region_list, ~str_trim(.x)),
      region_list = map(region_list, ~.x[.x %in% regions]),
      k = lengths(region_list)
    ) %>%
    filter(k >= 1)

  region_pair_share <- module_counts %>%
    filter(k >= 2) %>%
    mutate(region_pairs = map(region_list, ~{
      combn(.x, 2, simplify = FALSE) |>
        map_dfr(~tibble(from=.x[1], to=.x[2]))
    })) %>%
    unnest(region_pairs) %>%
    group_by(from, to) %>%
    summarise(value = sum(n_promoters), .groups = "drop") %>%
    filter(value > 0)

  region_singletons <- module_counts %>%
    filter(k == 1) %>%
    transmute(region = map_chr(region_list, 1), value = n_promoters) %>%
    group_by(region) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    right_join(tibble(region = regions), by = "region") %>%
    mutate(value = replace_na(value, 0))

  # 确保 region_singletons 包含所有 region（没的补 0）
  region_singletons <- tibble(region = regions) %>%
    left_join(region_singletons, by = "region") %>%
    mutate(value = replace_na(value, 0))

  max_single <- max(region_singletons$value, 1)
  ylim_top   <- max_single * 1.25

  pdf(file.path(outdir, "brain.complete_loss_table.circos_plot.pdf"), width = 6, height = 6)

  circos.clear()
  circos.par(start.degree = 90, gap.degree = 6)

  chordDiagram(
    x = region_pair_share,
    order = regions,
    transparency = 0.4,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.14)
  )

  circos.trackPlotRegion(
    track.index = 1,
    ylim = c(0, ylim_top),
    bg.border = NA,
    panel.fun = function(x, y) {
      reg  <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      v <- region_singletons$value[match(reg, region_singletons$region)]

      # 外圈 bar（只表示大小，不标数字）
      circos.rect(
        xleft  = xlim[1],
        ybottom= 0,
        xright = xlim[2],
        ytop   = v,
        border = NA
      )

      # brain region 名字（外侧）
      circos.text(
        x = mean(xlim),
        y = max_single * 1.18,
        labels = reg,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 0.7
      )
    }
  )

  dev.off()
}


# 2）heatmap
build.heatmap <- function(complete_loss_table, regions, outdir = ".") {
  # 1) gene × region 0/1（OR 聚合）
  gene_region_mat <- complete_loss_table %>%
    select(promoter_gene_id, all_of(regions)) %>%
    mutate(gene = str_split(promoter_gene_id, ",")) %>%
    unnest(gene) %>%
    mutate(gene = str_trim(gene)) %>%
    group_by(gene) %>%
    summarise(across(all_of(regions), ~ as.integer(any(.x == 1))), .groups = "drop") %>%
    mutate(n_regions_present = rowSums(across(all_of(regions))))

  # 2) 按 “1区→2区→…→6区” 排序
  #    同一层里，再按 presence pattern 排序（会更“成块”）
  gene_region_mat2 <- gene_region_mat %>%
    mutate(pattern = do.call(paste0, c(across(all_of(regions)), sep = ""))) %>%  # 例如 101010
    arrange(n_regions_present, pattern)

  # 3) 转 matrix（列固定顺序；不 cluster）
  mat <- gene_region_mat2 %>%
    select(gene, all_of(regions)) %>%
    column_to_rownames("gene") %>%
    as.matrix()

  # 可选：加一个行注释，显示每个 gene 属于 1/2/3/4/5/6 哪一层（不显示行名也能看到阶梯分段）
  ann_row <- gene_region_mat2 %>%
    select(gene, n_regions_present) %>%
    mutate(n_regions_present = factor(n_regions_present, levels = 1:6)) %>%
    column_to_rownames("gene")

  pdf(file.path(outdir, "brain.complete_loss_table.gene_region_heatmap.pdf"),
      width = 4, height = 8)

  pheatmap(
    mat[, rev(colnames(mat))],
    color = c("#F3F4F4", "#FDB5CE"),
    cluster_rows = FALSE,             # 关键：不聚类
    cluster_cols = FALSE,             # 关键：列不聚类（固定脑区顺序）
    show_rownames = FALSE,
    show_colnames = TRUE,
    fontsize_col = 10,
    border_color = NA,
    annotation_row = ann_row,         # 可选但很推荐：显示 1–6 阶梯分层
    main = ""
  )

  dev.off()
}

# jaccard similarity heatmap
build.jaccard <- function(complete_loss_table, regions, outdir = ".") {
  # 2.1 先把每个 region 的 hub genes 拿出来
  genes_by_region <- map(
    regions,
    ~ complete_loss_table %>%
        filter(.data[[.x]] == 1) %>%
        pull(promoter_gene_id) %>%
        unique()
  )

  names(genes_by_region) <- regions

  # 2.2 计算所有 region-pairs 的 Jaccard
  region_pairs <- combn(regions, 2, simplify = FALSE)

  jaccard_table <- map_dfr(region_pairs, function(pair) {
    A <- genes_by_region[[pair[1]]]
    B <- genes_by_region[[pair[2]]]

    tibble(
      region1 = pair[1],
      region2 = pair[2],
      n_A = length(A),
      n_B = length(B),
      intersection = length(intersect(A, B)),
      union = length(union(A, B)),
      jaccard = intersection / union
    )
  })

  jaccard_table %>% arrange(desc(jaccard))

  full_matrix_data <- bind_rows(
    jaccard_table %>% select(region1, region2, jaccard),
    jaccard_table %>% select(region1 = region2, region2 = region1, jaccard), # 交换顺序补全另一半
    tibble(region1 = unique(c(jaccard_table$region1, jaccard_table$region2)), 
          region2 = region1, jaccard = 1) # 补全对角线（自对比为1）
  )

  # 2. 转换为矩阵
  jaccard_matrix <- full_matrix_data %>%
    pivot_wider(names_from = region2, values_from = jaccard) %>%
    column_to_rownames("region1") %>%
    as.matrix()
  
  jaccard_matrix = jaccard_matrix[, rownames(jaccard_matrix)]

  pdf(file.path(outdir, "brain.complete_loss_table.gene_region_heatmap.jaccard_pairwise.pdf"),
      width = 6, height = 6)

    # 3. 绘制热图
    pheatmap(jaccard_matrix, 
            clustering_distance_rows = "euclidean",
            clustering_distance_cols = "euclidean",
            clustering_method = "complete",
            color = colorRampPalette(c("white", "#FEB05D", "#FF0087", "#FF0087", "#FF0087"))(100), # 漂亮的蓝色调
            display_numbers = TRUE, # 是否在格子里显示数字
            number_format = "%.2f",
            fontsize_number = 15,
            main = "Jaccard Similarity of Hubs across Brain Regions")
  dev.off()


}


# gene expression data compare between different gene set for each specific region
build.gene.expression.in.plot = function(gene_expression, gene_sample, hubs, pairs, regions, target_region, expr_threshold = 1, outdir = ".", select_top_genes = FALSE) {  
  print(target_region)
  ## =========================
  ## 参数
  ## =========================
  # target_region <- "MiddleFrontalGyrus"
  # expr_threshold <- 1          # 你要的 TPM 阈值（按 max TPM 过滤）
  set.seed(1000)

  ## =========================
  ## 1) gene_expression: 宽表 -> 长表，并与 gene_sample 合并
  ##    - gene_expression: ensembl_gene_id, hgnc_symbol, + 很多 sample 列 (HA74, HC34, ...)
  ##    - gene_sample$name1: 对应这些 sample 列名 (HA3, HB3, ...)
  ## =========================
  gene_sample_clean <- gene_sample %>%
    filter(name3 %in% regions) %>%                    # 去掉 regions 里没有的
    mutate(
      name3 = factor(name3, levels = regions)          # 按指定顺序排序
    ) %>%
    arrange(name3) %>%  # 确保 sample1 按 regions 排序
    select(2:4)

  ge_long_ori <- gene_expression %>%
    pivot_longer(
      cols = -c(ensembl_gene_id, hgnc_symbol),
      names_to = "sample",
      values_to = "expr"
    ) %>%
    mutate(
      ensembl_gene_id = as.character(ensembl_gene_id),
      sample = as.character(sample),
      expr = as.numeric(expr)
    ) %>%
    inner_join(
      gene_sample_clean %>%
        transmute(
          sample = as.character(name1),
          region = as.character(name3),
          name2 = as.character(name2)
        ),
      by = "sample"
    )

  ge_long <- ge_long_ori %>%
    mutate(
      sample_grp = str_extract(sample, "^HA|^HB|^HC|^HD")
    ) %>%
    filter(!is.na(sample_grp)) %>%
    group_by(ensembl_gene_id, hgnc_symbol, sample_grp, region) %>% #, name2
    summarise(
      expr = sum(expr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(sample = sample_grp)

  ge_long_mean <- ge_long_ori %>%
    mutate(
      sample_grp = str_extract(sample, "^HA|^HB|^HC|^HD")
    ) %>%
    filter(!is.na(sample_grp)) %>%
    group_by(ensembl_gene_id, hgnc_symbol, sample_grp, region) %>% #, name2
    summarise(
      expr = mean(expr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(sample = sample_grp)

  write.table(ge_long_mean, file=file.path(outdir, "gene_expression_merged.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


  ## =========================
  ## 2) complete_loss_table: promoter_gene_id 拆分成单基因并去版本号
  ##    例如 "ENSG...20" -> "ENSG..."; 多基因用逗号分隔
  ## =========================


  ## =========================
  ## 2.5) 定义：target-present vs target-absent（都限定在 complete_loss_table 的基因范围内）
  ## =========================
  # present_regions 可能包含多个区域；这里做“逗号边界”的包含匹配
  genes_target_present <- hubs

  genes_target_absent <- pairs

  cl_genes <- unique(c(hubs, pairs))

  ## =========================
  ## 2) gene_expression_rest：max(TPM) >= expr_threshold 且不与 complete_loss_table 重合
  ## =========================
  expr_gene_max <- ge_long %>%
    group_by(ensembl_gene_id) %>%
    summarise(
      max_tpm = max(expr, na.rm = TRUE),
      .groups = "drop"
    )

  genes_pass_threshold <- expr_gene_max %>%
    filter(is.finite(max_tpm), max_tpm >= expr_threshold) %>%
    pull(ensembl_gene_id)

  genes_rest <- setdiff(genes_pass_threshold, cl_genes)

  gene_expression_rest <- ge_long %>%
    filter(ensembl_gene_id %in% genes_rest)

  ## =========================
  ## 3) gene_expression_random：从 rest 随机抽 1500 个基因
  ## =========================
  n_random <- 1500
  if (length(genes_rest) < n_random) {
    warning(sprintf("genes_rest 只有 %d 个，少于 %d；将抽取全部 genes_rest。", length(genes_rest), n_random))
    genes_random <- genes_rest
  } else {
    genes_random <- sample(genes_rest, n_random, replace = FALSE)
  }

  gene_expression_random <- ge_long %>%
    filter(ensembl_gene_id %in% genes_random)

  ## =========================
  ## 4) 只在 target_region 样本中取表达，并打上分组标签
  ##    比较四组：target-present / target-absent / random-1500 / rest
  ## =========================
  expr_in_target_samples <- ge_long %>%
    filter(region == target_region) %>%
    mutate(
      group = case_when(
        ensembl_gene_id %in% genes_target_present ~ "Hubs",
        ensembl_gene_id %in% genes_target_absent  ~ "Pairs",
        ensembl_gene_id %in% genes_random         ~ "Random",
        ensembl_gene_id %in% genes_rest           ~ "All",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(group))

  # 为了展示更稳健，建议用 log1p(TPM)
  expr_in_target_samples <- expr_in_target_samples %>%
    mutate(expr_log = log1p(expr))

  p5  <- quantile(expr_in_target_samples$expr_log, 0.05, na.rm = TRUE)
  p95 <- quantile(expr_in_target_samples$expr_log, 0.95, na.rm = TRUE)

  expr_in_target_samples <- expr_in_target_samples %>%
    mutate(
      expr_log_w = pmin(pmax(expr_log, p5), p95)
    )
  ## =========================
  ## 5) 画图：四组在 target_region 样本里的表达分布
  ##    这里是把 “gene × sample” 的点作为分布（最贴近你原始 boxplot 逻辑）
  ## =========================
  expr_in_target_samples$group <- factor(
    expr_in_target_samples$group,
    levels = c(
      "Hubs",
      "Pairs",
      "Random",
      "All"
    )
  )

  ## =========================
  ## 6) （可选但强烈推荐）每个 sample 单独比较：facet by sample
  ## =========================
  grp_colors <- c("#FF0087","#FEB05D", "#B7B7B7", "#1B211A")
  expr_in_target_samples$sample <- gsub("H", "Human_", expr_in_target_samples$sample)
  ref <- levels(expr_in_target_samples$group)[1]

  comparisons_list <- lapply(
    levels(expr_in_target_samples$group)[-1],
    function(g) c(ref, g)
  )

  p2 <- ggplot(expr_in_target_samples, aes(x = group, y = expr_log_w, fill = group)) +
    geom_boxplot(alpha = 0.9, color = "black") +
    facet_wrap(~ sample, ncol = 6, scales = "free_y") +
    scale_fill_manual(values = grp_colors) +
    scale_color_manual(values = grp_colors) +
    theme_classic(base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.text = element_text(color = "black", face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5)
    ) +
    labs(x = "", y = "log1p(TPM)", title = paste0(target_region))+
    # coord_cartesian(ylim = c(quantile(expr_in_target_samples$expr_log, 0.05, na.rm = TRUE), quantile(expr_in_target_samples$expr_log, 0.95, na.rm = TRUE))) +
    stat_compare_means(
      comparisons = comparisons_list,
      method = "wilcox.test",
      label = "p.format",   # 或 "p.signif"
      size = 3,
      hide.ns = FALSE
    )

  ggsave(plot=p2, file = file.path(outdir, paste("brain.selected_region_expression",target_region,"pdf",sep=".")), width = 6.5, height = 3.5)
  write.table(gene_expression[gene_expression$ensembl_gene_id %in% genes_target_present & !is.na(gene_expression$hgnc_symbol) & gene_expression$hgnc_symbol != "","hgnc_symbol",drop=FALSE], 
    file = file.path(outdir, paste("brain.selected_region_expression",target_region,"txt",sep=".")),row.names=FALSE, quote=FALSE, col.names=FALSE, sep="\t")
}

## barplot
build.barplot <- function(complete_loss_table, regions, outdir = ".") {
  p1 <- complete_loss_table %>%
    count(n_regions_present) %>%
    ggplot(aes(x = factor(n_regions_present, levels=sort(unique(complete_loss_table$n_regions_present))), y = n)) +
    geom_col(fill = "#FEB05D", color = "white") +
    geom_text(aes(label = n), vjust = 0, size = 3) +
    labs(x = "Number of regions Hub presented", y = "Hub number") +
    theme_classic(base_family = "Helvetica") +
    theme(
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        strip.text = element_text(color = "black", face = "bold"),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5)
    )

  ggsave(plot=p1, file = file.path(outdir, paste("brain.complete_loss_table.barplot.overlap.pdf",sep=".")), width = 2.8, height = 2.5)

  region_counts <- complete_loss_table %>%
    summarise(across(all_of(regions), sum)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "region",
      values_to = "n_hubs"
    )

  p2 <- ggplot(region_counts, aes(x = reorder(region, -n_hubs), y = n_hubs)) +
    geom_col(fill = "#DE1A58", color = "white") +
    geom_text(aes(label = n_hubs), vjust = 0, size = 3) +
    labs(x = "", y = "Hub number") +
    theme_classic(base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        strip.text = element_text(color = "black", face = "bold"),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5)
    )

  ggsave(plot=p2, file = file.path(outdir, paste("brain.complete_loss_table.barplot.regions.pdf",sep=".")), width = 3, height = 3)


# 计算每个脑区的特异性百分比
spec_percentages <- sapply(regions, function(r) {
  total <- sum(complete_loss_table[[r]] == 1)
  spec  <- sum(complete_loss_table[[r]] == 1 & complete_loss_table$n_regions_present == 1)
  
  # 返回百分比（保留2位小数）
  round(spec / total * 100, 2)
})

# 查看结果
print(spec_percentages)


}



##Export gene lists based on hub classification
build.go.txt <- function(complete_loss_table, regions, outdir = ".") {
  ##This is only for class 1 and 3
  # 1. 执行分类逻辑
  classified_table <- complete_loss_table %>%
    mutate(hub_class = case_when(
      # 优先判定 Class 3
      n_regions_present >= 3 ~ "class3_common3",
      
      # 判定 3 个 Strong Pairs (Class 2)
      # 使用 grepl 确保即使顺序不同也能匹配，或者直接匹配你表中生成的标准字符串
      n_regions_present == 2 & (present_regions == "MiddleFrontalGyrus,SubstantiaNigra" | 
                                present_regions == "SubstantiaNigra,MiddleFrontalGyrus") ~ "class2_MidSub",
      
      n_regions_present == 2 & (present_regions == "Caudate,Hippocampus" | 
                                present_regions == "Hippocampus,Caudate") ~ "class2_CauHipp",
      
      n_regions_present == 2 & (present_regions == "ParietalLobe,SuperiorTemporalGyri" | 
                                present_regions == "SuperiorTemporalGyri,ParietalLobe") ~ "class2_PariSup",
      
      # 判定单脑区特异 (Class 1)
      n_regions_present == 1 ~ paste0("class1_", present_regions),
      
      # 其他不符合上述规则的 2 脑区组合
      TRUE ~ "class_others"
    ))

  classified_table %>%
    group_split(hub_class) %>%
    walk(~ {
      current_class <- unique(.x$hub_class)
      # 提取基因 ID，考虑到一行可能有多个基因（逗号分隔），先展开再取唯一值
      genes <- .x %>%
        pull(promoter_gene_name) %>%
        strsplit(",") %>% 
        unlist() %>%
        unique()
      
      # 写入文件
      write_lines(genes, file.path(outdir, paste("brain.selected_region_overlap",current_class,"txt",sep=".")))
      message(paste("Class 1 and 3 Exported:", current_class, "with", length(genes), "unique genes."))
    })

  classified_table %>%
    group_split(hub_class) %>%
    walk(~ {
      current_class <- unique(.x$hub_class)
      # 提取基因 ID，考虑到一行可能有多个基因（逗号分隔），先展开再取唯一值
      genes <- .x %>%
        pull(promoter_id) %>%
        unique()
      
      # 写入文件
      write_lines(genes, file.path(outdir, paste("brain.selected_region_overlap",current_class,"promoterid",sep=".")))
    })

  print(table(classified_table$hub_class))

  ##This is only for class 2
  classified_table <- complete_loss_table %>%
    mutate(hub_class = case_when(
      # 优先判定：核心对 (Class 2)
      ParietalLobe == 1 & SuperiorTemporalGyri == 1 ~ "class2only_PariSup",
      Caudate == 1 & Hippocampus == 1 ~ "class2only_CauHipp",    
      MiddleFrontalGyrus == 1 & SubstantiaNigra == 1 ~ "class2only_MidSub",
      
      # 其余情况（不包含核心对的多脑区组合）
      TRUE ~ "class_others"
    ))

  # 3. 按分类批量导出 txt 文件
  # 我们先根据 hub_class 分组，然后对每个组提取 promoter_gene_id 并去重
  classified_table %>%
    group_split(hub_class) %>%
    walk(~ {
      current_class <- unique(.x$hub_class)
      # 提取基因 ID，考虑到一行可能有多个基因（逗号分隔），先展开再取唯一值
      genes <- .x %>%
        pull(promoter_gene_name) %>%
        strsplit(",") %>% 
        unlist() %>%
        unique()
      
      # 写入文件
      write_lines(genes, file.path(outdir, paste("brain.selected_region_overlap",current_class,"txt",sep=".")))
      message(paste("Class 2 Exported:", current_class, "with", length(genes), "unique genes."))
    })
  print(table(classified_table$hub_class))

}

build.go.txt.to.bed <- function(bed_list, outdir = ".") {
  bed_list_grlist <- GRangesList(bed_list)
  bed_all_gr <- unlist(bed_list_grlist, use.names = FALSE)
  mcols(bed_all_gr)$clean_name <- sub("_\\d+$", "", mcols(bed_all_gr)$name)
  
  # 1. 定义所有需要处理的类别列表
  class_list <- c("class3_common3", "class2_MidSub", "class1_MiddleFrontalGyrus", "class1_SubstantiaNigra")

  # 3. 开始循环处理
  for (current_class in class_list) {
    keptid <- read.delim(file.path(outdir, paste("brain.selected_region_overlap",current_class,"promoterid",sep=".")))
    message("正在处理: ", current_class)
    
    # --- 步骤 A: 根据 ifelse 逻辑确定过滤条件 ---
    if (current_class == "class3_common3") {
      # 提取所有包含这些 ID 的区域 (不限脑区)
      target_gr <- bed_all_gr[mcols(bed_all_gr)$clean_name %in% keptid[,1] & 
                              mcols(bed_all_gr)$region %in% c("SubstantiaNigra", "MiddleFrontalGyrus")]
      
    } else if (current_class == "class2_MidSub") {
      # 提取在 Mid 或 Sub 中的区域
      target_gr <- bed_all_gr[mcols(bed_all_gr)$clean_name %in% keptid[,1] & 
                              mcols(bed_all_gr)$region %in% c("SubstantiaNigra", "MiddleFrontalGyrus")]
      
    } else if (current_class == "class1_MiddleFrontalGyrus") {
      # 仅提取 MiddleFrontalGyrus
      target_gr <- bed_all_gr[mcols(bed_all_gr)$clean_name %in% keptid[,1] & 
                              mcols(bed_all_gr)$region == "MiddleFrontalGyrus"]
      
    } else if (current_class == "class1_SubstantiaNigra") {
      # 仅提取 SubstantiaNigra
      target_gr <- bed_all_gr[mcols(bed_all_gr)$clean_name %in% keptid[,1] & 
                              mcols(bed_all_gr)$region == "SubstantiaNigra"]
    }

    sorted_gr <- sort(target_gr)

    if(current_class %in% c("class3_common3", "class2_MidSub")){
      export(sorted_gr[sorted_gr$region=="SubstantiaNigra"],con=file.path(outdir, paste("brain.selected_region_overlap",current_class,"SubstantiaNigra.bed",sep=".")))
      export(sorted_gr[sorted_gr$region=="MiddleFrontalGyrus"],con=file.path(outdir, paste("brain.selected_region_overlap",current_class,"MiddleFrontalGyrus.bed",sep=".")))
    }

    # 转换为 data.frame 进行复杂的“位置+名称”去重和重新编号
    final_df <- as.data.frame(sorted_gr) %>%
      arrange(clean_name, seqnames, start) %>%
      distinct(seqnames, start, end, clean_name) %>%
      group_by(clean_name) %>%
      mutate(name = paste0(clean_name, "_", row_number())) %>%
      ungroup()

    # 4. 转回 GRanges 对象
    final_gr <- makeGRangesFromDataFrame(final_df, keep.extra.columns = TRUE)
    export(final_gr,con=file.path(outdir, paste("brain.selected_region_overlap",current_class,"bed",sep=".")))

    # 查看结果
    print(final_gr)
  }
}


##################################
## Generate complete loss table
loops_list <- list()
for (samplename in regions) {
  load(file.path(indir, paste("multiple_result",samplename,"replicable_hubs_loops.RData",sep=".")))
  loops_list[[samplename]] <- observed_hubs_per_promoter_sub
}
print(paste0("hub number in each sample."))
print(sapply(loops_list,function(hublist){length(unique(hublist$promoter_id))}))

build.complete.loss.table(loops_list, outdir=indir, replicates="union", min_hub_size=3, verbose=TRUE) #Use it if there are NO replicates in the loop list
#build.complete.loss.table(loops_list, outdir=indir, replicates="intersect", min_hub_size=3, verbose=TRUE) #Use it if there are replicates in the loop list

##################################
# table used for plotting
load(file.path(indir, "brain.complete_loss_table.union.RData"))
complete_loss_table <- complete_loss_table_list$complete_loss_table
# build.circos.plot(complete_loss_table, regions <- c("Caudate","Hippocampus","MiddleFrontalGyrus","ParietalLobe","SubstantiaNigra","SuperiorTemporalGyri"), outdir = gsub("ProcessedData","Output",indir))
build.heatmap(complete_loss_table, regions <- c("SubstantiaNigra","MiddleFrontalGyrus","Caudate","Hippocampus","ParietalLobe","SuperiorTemporalGyri"), outdir = gsub("ProcessedData","Output",indir))
build.barplot(complete_loss_table, regions <- c("SubstantiaNigra","MiddleFrontalGyrus","Caudate","Hippocampus","ParietalLobe","SuperiorTemporalGyri"), outdir = gsub("ProcessedData","Output",indir))
build.jaccard(complete_loss_table, regions <- c("SubstantiaNigra","MiddleFrontalGyrus","Caudate","Hippocampus","ParietalLobe","SuperiorTemporalGyri"), outdir = gsub("ProcessedData","Output",indir))
build.go.txt(complete_loss_table, regions <- c("SubstantiaNigra","MiddleFrontalGyrus","Caudate","Hippocampus","ParietalLobe","SuperiorTemporalGyri"), outdir = gsub("ProcessedData","Output",indir))

##################################
# txt to bed
bed_list <- list()
for (samplename in regions) {
  bed_list[[samplename]] <- import(file.path(indir, paste("multiple_result",samplename,"replicable_hubs_bin_union.bed",sep=".")))
  bed_list[[samplename]]$region <- samplename
}
build.go.txt.to.bed(bed_list,outdir = gsub("ProcessedData","Output",indir))

##################################
## gene expression data
gene_expression <- read.delim("~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged/RNAseq/human_TPM.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene_sample <- read.delim("~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged/RNAseq/human_samples_name.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
regions <- c("SubstantiaNigra","MiddleFrontalGyrus","Caudate","Hippocampus","ParietalLobe","SuperiorTemporalGyri")

for (target_region in regions) {
  pairs_df <- read.delim(file.path(indir, paste("multiple_result",paste0(target_region),"pairwise_bin_union","txt",sep=".")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  hubs_df <- read.delim(file.path(indir, paste("multiple_result",paste0(target_region),"replicable_hubs_bin_intersect","txt",sep=".")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  build.gene.expression.in.plot(gene_expression, gene_sample, hubs=unique(gsub("\\.[0-9]+", "", unlist(strsplit(na.omit(hubs_df$promoter_gene_id), ",")))), pairs=unique(gsub("\\.[0-9]+", "", unlist(strsplit(na.omit(pairs_df$promoter_gene_id), ",")))), regions=regions, target_region=target_region, expr_threshold = 1,outdir = gsub("ProcessedData","Output",indir))
}






# for (target_region in c(unique(complete_loss_table$present_regions),"all")) {
#   build.gene.expression.out.plot(gene_expression, gene_sample, complete_loss_table, regions=regions, target_region=target_region, expr_threshold = 2,outdir = gsub("ProcessedData","Output",indir))
# }


# # 5) 产出 B：Hub signature gain/loss table
# hub_signature_one <- function(loops_df, sample_id, min_hub_size = 2) {
#   loops_df %>%
#     filter(pair_type == "PE") %>%
#     group_by(promoter_id) %>%
#     summarise(
#       sample = sample_id,
#       enhancers = list(sort(unique(enhancer_id))),
#       hub_size = length(enhancers[[1]]),
#       hub_strength = sum(counts, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     filter(hub_size >= min_hub_size) %>%
#     mutate(
#       signature = paste0(promoter_id, "|", vapply(enhancers, paste, collapse = ",", FUN.VALUE = character(1))),
#       region = regions
#     )
# }

# hub_sig_by_sample <- imap_dfr(loops_list, ~ hub_signature_one(.x, .y, min_hub_size = 2))

# hub_sig_by_region <- hub_sig_by_sample %>%
#   group_by(region, signature, promoter_id) %>%
#   summarise(
#     present = TRUE,
#     hub_size = max(hub_size),
#     hub_strength = sum(hub_strength),
#     .groups = "drop"
#   )

# hub_sig_presence <- hub_sig_by_region %>%
#   select(signature, promoter_id, region, present) %>%
#   pivot_wider(names_from = region, values_from = present, values_fill = FALSE)

# # 6) 产出 C：Rewiring table
# jaccard <- function(a, b) length(intersect(a,b)) / length(union(a,b))

# rewiring_pair <- function(regionA, regionB, hub_sig_by_region) {
#   A <- hub_sig_by_region %>% filter(region == regionA) %>% select(promoter_id, enhancers)
#   B <- hub_sig_by_region %>% filter(region == regionB) %>% select(promoter_id, enhancers)
#   AB <- inner_join(A, B, by="promoter_id", suffix=c("_A","_B"))
#   AB %>%
#     mutate(
#       jaccard = map2_dbl(enhancers_A, enhancers_B, jaccard),
#       shared_n = map2_int(enhancers_A, enhancers_B, ~ length(intersect(.x,.y))),
#       size_A = map_int(enhancers_A, length),
#       size_B = map_int(enhancers_B, length)
#     ) %>%
#     select(promoter_id, jaccard, shared_n, size_A, size_B)
# }


# # loops_list: named list, e.g. list(PFC=loops_pfc, HIP=loops_hip)
# make_promoter_table <- function(loops_df, sample_name, min_hub_size = 2) {
#   loops_df %>%
#     group_by(promoter_id, promoter_gene_id) %>%
#     summarise(
#       hub_size = n_distinct(enhancer_id),
#       hub_strength = sum(counts, na.rm = TRUE),
#       has_hub = hub_size >= min_hub_size,
#       .groups = "drop"
#     ) %>%
#     mutate(sample = sample_name)
# }

# promoter_by_sample <- imap_dfr(loops_list, ~ make_promoter_table(.x, .y, min_hub_size = 2))

# # 两个脑区示例：PFC vs HIP
# loss_table_pair <- promoter_by_sample %>%
#   select(sample, promoter_id, promoter_gene_id, has_hub, hub_size, hub_strength) %>%
#   pivot_wider(
#     names_from = sample,
#     values_from = c(has_hub, hub_size, hub_strength),
#     values_fill = 0
#   ) %>%
#   mutate(
#     classification = case_when(
#       has_hub_PFC == 1 & has_hub_HIP == 0 ~ "loss_in_HIP",
#       has_hub_PFC == 0 & has_hub_HIP == 1 ~ "loss_in_PFC",
#       has_hub_PFC == 1 & has_hub_HIP == 1 ~ "retained",
#       TRUE ~ "none"
#     )
#   )






















# ###############################################################################
# ### 1. 读取 hicdcplus differential 结果，生成 loop_DE
# ###############################################################################

# #' Read hicdcplus differential bedpe file
# #'
# #' @param diff_file  path to hicdcplus differential bedpe
# #' @param logFC_col  column name for log2FC (or other signed effect size)
# #' @param p_col      column name for p-value
# #' @param q_col      column name for FDR / q-value
# #'
# #' @return data.frame with: loop_id, chr1, start1, end1, chr2, start2, end2,
# #'         log2FC, pvalue, FDR, direction
# #'
# read_hicdc_diff <- function(diff_file,
#                             logFC_col = "log2FC",
#                             p_col    = "pvalue",
#                             q_col    = "qvalue") {
  
#   df <- read.table(diff_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
#   # 根据你的实际列名修改这里
#   chr1   <- df[["chr1"]]
#   start1 <- df[["start1"]]
#   end1   <- df[["end1"]]
#   chr2   <- df[["chr2"]]
#   start2 <- df[["start2"]]
#   end2   <- df[["end2"]]
  
#   log2FC <- df[[logFC_col]]
#   pval   <- df[[p_col]]
#   qval   <- df[[q_col]]
  
#   loop_id <- paste(chr1, start1, chr2, start2, sep = "_")
  
#   loop_DE <- tibble(
#     loop_id = loop_id,
#     chr1    = chr1,
#     start1  = start1,
#     end1    = end1,
#     chr2    = chr2,
#     start2  = start2,
#     end2    = end2,
#     log2FC  = log2FC,
#     pvalue  = pval,
#     FDR     = qval
#   ) %>%
#     mutate(
#       direction = case_when(
#         !is.na(FDR) & FDR < 0.05 & log2FC > 0  ~ "up",
#         !is.na(FDR) & FDR < 0.05 & log2FC < 0  ~ "down",
#         TRUE                                   ~ "ns"
#       )
#     )
  
#   return(loop_DE)
# }


# ###############################################################################
# ### 2. hub_A / hub_B 基于 promoter 做结构匹配（Jaccard、gain/loss）
# ###############################################################################

# #' Match hubs between two conditions by promoter_id
# #'
# #' @param hub_A GRanges for condition A hubs (each row = a locus in a hub)
# #' @param hub_B GRanges for condition B hubs
# #' @param hub_col      column name for hub id (same semantic in A/B, e.g. "hub_id")
# #' @param promoter_col column for promoter_id
# #' @param is_promoter_col logical column marking promoter loci
# #'
# #' @return data.frame: one row per promoter, with enhancer sets in A/B and Jaccard
# #'
# match_hubs_between_conditions <- function(hub_A,
#                                           hub_B,
#                                           hub_col         = "hub_id",
#                                           promoter_col    = "promoter_id",
#                                           is_promoter_col = "is_promoter") {
  
#   # 确保 locus_id 存在
#   if (!"locus_id" %in% colnames(mcols(hub_A))) {
#     mcols(hub_A)$locus_id <- paste0("A_", seq_along(hub_A))
#   }
#   if (!"locus_id" %in% colnames(mcols(hub_B))) {
#     mcols(hub_B)$locus_id <- paste0("B_", seq_along(hub_B))
#   }
  
#   mA <- mcols(hub_A)
#   mB <- mcols(hub_B)
  
#   promoters_common <- intersect(mA[[promoter_col]], mB[[promoter_col]])
  
#   out_list <- list()
  
#   for (p in promoters_common) {
    
#     subA <- hub_A[mA[[promoter_col]] == p]
#     subB <- hub_B[mB[[promoter_col]] == p]
    
#     # enhancer 集合（排除 promoter 自身）
#     enhA <- mcols(subA)$locus_id[!as.logical(mcols(subA)[[is_promoter_col]])]
#     enhB <- mcols(subB)$locus_id[!as.logical(mcols(subB)[[is_promoter_col]])]
    
#     enhA <- unique(enhA)
#     enhB <- unique(enhB)
    
#     nA <- length(enhA)
#     nB <- length(enhB)
#     n_shared <- length(intersect(enhA, enhB))
#     J <- if ((nA + nB - n_shared) == 0) NA_real_ else n_shared / (nA + nB - n_shared)
    
#     out_list[[length(out_list) + 1]] <- tibble(
#       promoter = p,
#       enhA     = list(enhA),
#       enhB     = list(enhB),
#       nA       = nA,
#       nB       = nB,
#       n_shared = n_shared,
#       Jaccard  = J
#     )
#   }
  
#   hub_pairs <- bind_rows(out_list)
#   return(hub_pairs)
# }


# ###############################################################################
# ### 3. 将 loop_DE 映射到 promoter/hub（anchor 落在 hub 的 loci 上）
# ###############################################################################

# # 思路：对每个 promoter 的 hub，把它在 A/B 中所有 loci 合并成一个 union GRanges，
# # 然后看哪些 loops 的两个锚点都落在这个 union 上，就认为该 loop 属于这个 promoter 的 hub 对。

# #' Build loop → promoter mapping using hub loci (union of A and B)
# #'
# #' @param loop_DE data.frame from read_hicdc_diff
# #' @param hub_A   GRanges (condition A hubs)
# #' @param hub_B   GRanges (condition B hubs)
# #' @param promoter_col column for promoter_id
# #'
# #' @return data.frame(loop_id, promoter)
# #'
# map_loops_to_promoters <- function(loop_DE,
#                                    hub_A,
#                                    hub_B,
#                                    promoter_col = "promoter_id") {
  
#   mA <- mcols(hub_A)
#   mB <- mcols(hub_B)
  
#   promoters_common <- intersect(mA[[promoter_col]], mB[[promoter_col]])
  
#   # 先把每个 promoter 相关的 loci（A+B）做 union
#   promoter_union_list <- lapply(promoters_common, function(p) {
#     loci_A <- hub_A[mA[[promoter_col]] == p]
#     loci_B <- hub_B[mB[[promoter_col]] == p]
#     c(loci_A, loci_B)
#   })
#   names(promoter_union_list) <- promoters_common
  
#   # 把 loops 的两个 anchor 变成 GRanges
#   gr1 <- GRanges(
#     seqnames = loop_DE$chr1,
#     ranges   = IRanges(loop_DE$start1, loop_DE$end1),
#     loop_id  = loop_DE$loop_id,
#     anchor   = "anchor1"
#   )
#   gr2 <- GRanges(
#     seqnames = loop_DE$chr2,
#     ranges   = IRanges(loop_DE$start2, loop_DE$end2),
#     loop_id  = loop_DE$loop_id,
#     anchor   = "anchor2"
#   )
  
#   # 分别计算每个 anchor 落在哪些 promoter union loci 上
#   loop_to_promoter_list <- list()
  
#   for (p in promoters_common) {
#     hub_union <- promoter_union_list[[p]]
#     ov1 <- findOverlaps(gr1, hub_union, ignore.strand = TRUE)
#     ov2 <- findOverlaps(gr2, hub_union, ignore.strand = TRUE)
    
#     loop_ids1 <- mcols(gr1)$loop_id[queryHits(ov1)]
#     loop_ids2 <- mcols(gr2)$loop_id[queryHits(ov2)]
    
#     # 两个 anchor 都在该 promoter hub union 内的 loops
#     loops_both <- intersect(loop_ids1, loop_ids2)
    
#     if (length(loops_both) > 0) {
#       loop_to_promoter_list[[length(loop_to_promoter_list) + 1]] <-
#         tibble(loop_id = loops_both, promoter = p)
#     }
#   }
  
#   if (length(loop_to_promoter_list) == 0) {
#     warning("No loops mapped to any promoter hub.")
#     return(tibble(loop_id = character(), promoter = character()))
#   }
  
#   loop_to_promoter <- bind_rows(loop_to_promoter_list) %>%
#     distinct()
  
#   return(loop_to_promoter)
# }


# ###############################################################################
# ### 4. Hub-level 汇总：结构 + loop 差异 → 分类
# ###############################################################################

# #' Summarize differential hubs using structure (Jaccard, gain/loss) and loop_DE
# #'
# #' @param hub_pairs        data.frame from match_hubs_between_conditions
# #' @param loop_DE          data.frame from read_hicdc_diff
# #' @param loop_to_promoter mapping from map_loops_to_promoters
# #'
# #' @return data.frame per promoter/hub pair with differential summary
# #'
# summarize_hub_differential <- function(hub_pairs,
#                                        loop_DE,
#                                        loop_to_promoter) {
  
#   loop_ann <- loop_DE %>%
#     inner_join(loop_to_promoter, by = "loop_id")
  
#   out_list <- list()
  
#   for (i in seq_len(nrow(hub_pairs))) {
#     p    <- hub_pairs$promoter[i]
#     enhA <- unlist(hub_pairs$enhA[i])
#     enhB <- unlist(hub_pairs$enhB[i])
    
#     nA   <- hub_pairs$nA[i]
#     nB   <- hub_pairs$nB[i]
#     n_sh <- hub_pairs$n_shared[i]
#     J    <- hub_pairs$Jaccard[i]
    
#     n_gain <- length(setdiff(enhB, enhA))
#     n_loss <- length(setdiff(enhA, enhB))
    
#     loops_p <- loop_ann %>% filter(promoter == p)
    
#     # 统计该 promoter hub 对里的 up/down loops 数
#     n_up   <- sum(loops_p$direction == "up", na.rm = TRUE)
#     n_down <- sum(loops_p$direction == "down", na.rm = TRUE)
    
#     out_list[[length(out_list) + 1]] <- tibble(
#       promoter      = p,
#       Jaccard       = J,
#       nA            = nA,
#       nB            = nB,
#       n_shared      = n_sh,
#       n_gain        = n_gain,
#       n_loss        = n_loss,
#       n_up_loops    = n_up,
#       n_down_loops  = n_down
#     )
#   }
  
#   hub_diff <- bind_rows(out_list)
  
#   # 简单的分类规则（你可以后期微调）
#   hub_diff <- hub_diff %>%
#     mutate(
#       category = case_when(
#         (n_gain + n_loss) >= 2 & !is.na(Jaccard) & Jaccard < 0.5 ~ "Rewired",
#         n_gain >= 1 & n_up_loops > n_down_loops                  ~ "Gained",
#         n_loss >= 1 & n_down_loops > n_up_loops                  ~ "Lost",
#         TRUE                                                     ~ "Stable"
#       )
#     )
  
#   return(hub_diff)
# }


# ###############################################################################
# ### 5. 主 Pipeline 函数：把所有步骤串起来
# ###############################################################################

# #' Differential hub pipeline based on hicdcplus differential output and hubs
# #'
# #' @param hicdc_diff_file path to hicdcplus differential bedpe
# #' @param hub_A           GRanges for condition A hubs
# #' @param hub_B           GRanges for condition B hubs
# #'
# #' @return list(loop_DE, hub_pairs, loop_to_promoter, hub_diff)
# #'
# run_differential_hub_with_hicdc <- function(hicdc_diff_file,
#                                             hub_A,
#                                             hub_B,
#                                             logFC_col = "log2FC",
#                                             p_col    = "pvalue",
#                                             q_col    = "qvalue",
#                                             hub_col         = "hub_id",
#                                             promoter_col    = "promoter_id",
#                                             is_promoter_col = "is_promoter") {
  
#   # 1) 读 hicdcplus differential
#   loop_DE <- read_hicdc_diff(
#     diff_file = hicdc_diff_file,
#     logFC_col = logFC_col,
#     p_col     = p_col,
#     q_col     = q_col
#   )
  
#   # 2) hub structure matching (A vs B)
#   hub_pairs <- match_hubs_between_conditions(
#     hub_A,
#     hub_B,
#     hub_col         = hub_col,
#     promoter_col    = promoter_col,
#     is_promoter_col = is_promoter_col
#   )
  
#   # 3) map loops to promoter hubs (union of A + B loci)
#   loop_to_promoter <- map_loops_to_promoters(
#     loop_DE   = loop_DE,
#     hub_A     = hub_A,
#     hub_B     = hub_B,
#     promoter_col = promoter_col
#   )
  
#   # 4) summarize hub differential
#   hub_diff <- summarize_hub_differential(
#     hub_pairs        = hub_pairs,
#     loop_DE          = loop_DE,
#     loop_to_promoter = loop_to_promoter
#   )
  
#   return(list(
#     loop_DE          = loop_DE,
#     hub_pairs        = hub_pairs,
#     loop_to_promoter = loop_to_promoter,
#     hub_diff         = hub_diff
#   ))
# }

# # hub_A / hub_B 是你已经鉴定好的两个 cell line 的 hub GRanges
# # hicdc_diff_file 是 hicdcplus 的 differential bedpe

# res <- run_differential_hub_with_hicdc(
#   hicdc_diff_file = "GM12878_vs_NSC.hicdcplus.diff.bedpe",
#   hub_A           = hub_A,
#   hub_B           = hub_B,
#   logFC_col       = "log2FC",   # 根据你自己的列名改
#   p_col           = "pvalue",
#   q_col           = "qvalue",
#   hub_col         = "hub_id",
#   promoter_col    = "promoter_id",
#   is_promoter_col = "is_promoter"
# )

# # 最重要的结果：
# head(res$hub_diff)


# # 假设与约定（请对照你自己的对象改名）
# # hicdcplus differential bedpe（hicdc_diff_file）至少有列：
# # chr1, start1, end1, chr2, start2, end2,
# # log2FC, pvalue, qvalue
# # 如果你的列名不同，在 read_hicdc_diff() 里改参数即可。
# # hub_A / hub_B 是 GRanges，每一行是 hub 中的一个 bin/locus，mcols 至少包含：
# # hub_id：这个 locus 属于哪个 hub（同一个 hub 多行）
# # promoter_id：这个 locus 对应的 promoter（一个 hub 对应一个 promoter）
# # is_promoter：逻辑型，是否为 promoter，FALSE 视为 enhancer/bin
# # locus_id：该 locus 的一个 ID（如果没有脚本会自动生成）


# # gene expression data compare between different gene set for each specific region
# build.gene.expression.in.plot = function(gene_expression, gene_sample, complete_loss_table, pairs, regions, target_region, expr_threshold = 1,outdir = ".", select_top_genes = FALSE) {  
#   print(target_region)
#   ## =========================
#   ## 参数
#   ## =========================
#   # target_region <- "MiddleFrontalGyrus"
#   # expr_threshold <- 1          # 你要的 TPM 阈值（按 max TPM 过滤）
#   set.seed(1000)

#   ## =========================
#   ## 1) gene_expression: 宽表 -> 长表，并与 gene_sample 合并
#   ##    - gene_expression: ensembl_gene_id, hgnc_symbol, + 很多 sample 列 (HA74, HC34, ...)
#   ##    - gene_sample$name1: 对应这些 sample 列名 (HA3, HB3, ...)
#   ## =========================
#   gene_sample_clean <- gene_sample %>%
#     filter(name3 %in% regions) %>%                    # 去掉 regions 里没有的
#     mutate(
#       name3 = factor(name3, levels = regions)          # 按指定顺序排序
#     ) %>%
#     arrange(name3) %>%  # 确保 sample1 按 regions 排序
#     select(2:4)

#   ge_long_ori <- gene_expression %>%
#     pivot_longer(
#       cols = -c(ensembl_gene_id, hgnc_symbol),
#       names_to = "sample",
#       values_to = "expr"
#     ) %>%
#     mutate(
#       ensembl_gene_id = as.character(ensembl_gene_id),
#       sample = as.character(sample),
#       expr = as.numeric(expr)
#     ) %>%
#     inner_join(
#       gene_sample_clean %>%
#         transmute(
#           sample = as.character(name1),
#           region = as.character(name3),
#           name2 = as.character(name2)
#         ),
#       by = "sample"
#     )

#   ge_long <- ge_long_ori %>%
#     mutate(
#       sample_grp = str_extract(sample, "^HA|^HB|^HC|^HD")
#     ) %>%
#     filter(!is.na(sample_grp)) %>%
#     group_by(ensembl_gene_id, hgnc_symbol, sample_grp, region, name2) %>% 
#     summarise(
#       expr = sum(expr, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     rename(sample = sample_grp)

#   ## =========================
#   ## 2) complete_loss_table: promoter_gene_id 拆分成单基因并去版本号
#   ##    例如 "ENSG...20" -> "ENSG..."; 多基因用逗号分隔
#   ## =========================
#   cl_genes <- complete_loss_table %>%
#     mutate(row_id = row_number()) %>%
#     transmute(
#       row_id,
#       promoter_id,
#       present_regions,
#       promoter_gene_id_raw = promoter_gene_id
#     ) %>%
#     mutate(promoter_gene_id_raw = str_replace_all(promoter_gene_id_raw, "\\s+", "")) %>%
#     separate_rows(promoter_gene_id_raw, sep = ",") %>%
#     mutate(
#       ensembl_gene_id = str_replace(promoter_gene_id_raw, "\\..*$", "")  # 去掉版本号 .20
#     ) %>%
#     select(-promoter_gene_id_raw)

  # ## =========================
  # ## 2.5) （新增）同一 promoter_id 下只保留在 target_region 中最高表达的那个基因
  # ##      present/absent 仍然只由 promoter 的 present_regions 决定
  # ##      如果放开这个，就要把下面的2.5节注释掉
  # ## =========================

  # # 先把 promoter 的 has_target 定义在 promoter 维度（不是 gene 维度）
  # promoter_has_target <- cl_genes %>%
  #   distinct(promoter_id, present_regions) %>%
  #   mutate(
  #     present_regions2 = paste0(",", present_regions, ","),
  #     has_target = str_detect(present_regions2, paste0(",", target_region, ","))
  #   ) %>%
  #   select(promoter_id, present_regions, has_target)

  # # 在 target_region 的表达里，为每个 promoter_id 的候选基因打分（汇总所有 HA/HB/HC/HD 样本）
  # # 这里用 sum(expr) 作为“总体表达量”，你也可以改成 mean(expr)
  # promoter_gene_scores <- cl_genes %>%
  #   distinct(promoter_id, ensembl_gene_id) %>%
  #   inner_join(
  #     ge_long %>% filter(region == target_region),
  #     by = "ensembl_gene_id"
  #   ) %>%
  #   group_by(promoter_id, ensembl_gene_id) %>%
  #   summarise(total_expr = sum(expr, na.rm = TRUE), .groups = "drop")

  # # 每个 promoter 选 total_expr 最大的那个基因
  # top_gene_per_promoter <- promoter_gene_scores %>%
  #   group_by(promoter_id) %>%
  #   slice_max(total_expr, n = 1, with_ties = FALSE) %>%
  #   ungroup()

  # # 用 top gene 生成 present/absent 两组基因（present/absent 只看 promoter_has_target）
  # genes_target_present <- top_gene_per_promoter %>%
  #   inner_join(promoter_has_target, by = "promoter_id") %>%
  #   filter(has_target) %>%
  #   pull(ensembl_gene_id) %>%
  #   unique()

  # genes_target_absent <- top_gene_per_promoter %>%
  #   inner_join(promoter_has_target, by = "promoter_id") %>%
  #   filter(!has_target) %>%
  #   pull(ensembl_gene_id) %>%
  #   unique()


#   ## =========================
#   ## 0) 准备：全基因集合（来自表达矩阵）与 complete_loss_table 基因集合
#   ## =========================
#   genes_in_cl <- cl_genes %>%
#     distinct(ensembl_gene_id) %>%
#     pull(ensembl_gene_id)

#   genes_in_expr <- ge_long %>%
#     distinct(ensembl_gene_id) %>%
#     pull(ensembl_gene_id)

#   ## =========================
#   ## 2.5) 定义：target-present vs target-absent（都限定在 complete_loss_table 的基因范围内）
#   ## =========================
#   # present_regions 可能包含多个区域；这里做“逗号边界”的包含匹配
#   cl_gene_regions <- cl_genes %>%
#     distinct(ensembl_gene_id, present_regions) %>%
#     mutate(
#       present_regions2 = paste0(",", present_regions, ","),
#       has_target = str_detect(present_regions2, paste0(",", target_region, ",")),
#       has_overlap = str_detect(present_regions2, paste0(",", "MiddleFrontalGyrus", ","))
#     ) %>%
#     group_by(ensembl_gene_id) %>%
#     summarise(has_target = any(has_target), has_overlap = any(has_overlap), .groups = "drop")

#   genes_target_present <- cl_gene_regions %>%
#     filter(has_target) %>%
#     pull(ensembl_gene_id)

#   genes_target_absent = pairs

#   # genes_target_absent <- cl_gene_regions %>%
#   #   filter(!has_target & !has_overlap) %>%
#   #   pull(ensembl_gene_id)

#   ## =========================
#   ## 2) gene_expression_rest：max(TPM) >= expr_threshold 且不与 complete_loss_table 重合
#   ## =========================
#   expr_gene_max <- ge_long %>%
#     group_by(ensembl_gene_id) %>%
#     summarise(
#       max_tpm = max(expr, na.rm = TRUE),
#       .groups = "drop"
#     )

#   genes_pass_threshold <- expr_gene_max %>%
#     filter(is.finite(max_tpm), max_tpm >= expr_threshold) %>%
#     pull(ensembl_gene_id)

#   genes_rest <- setdiff(genes_pass_threshold, genes_in_cl)

#   gene_expression_rest <- ge_long %>%
#     filter(ensembl_gene_id %in% genes_rest)

#   ## =========================
#   ## 3) gene_expression_random：从 rest 随机抽 1500 个基因
#   ## =========================
#   n_random <- 1500
#   if (length(genes_rest) < n_random) {
#     warning(sprintf("genes_rest 只有 %d 个，少于 %d；将抽取全部 genes_rest。", length(genes_rest), n_random))
#     genes_random <- genes_rest
#   } else {
#     genes_random <- sample(genes_rest, n_random, replace = FALSE)
#   }

#   gene_expression_random <- ge_long %>%
#     filter(ensembl_gene_id %in% genes_random)

#   ## =========================
#   ## 4) 只在 target_region 样本中取表达，并打上分组标签
#   ##    比较四组：target-present / target-absent / random-1500 / rest
#   ## =========================
#   expr_in_target_samples <- ge_long %>%
#     filter(region == target_region) %>%
#     mutate(
#       group = case_when(
#         ensembl_gene_id %in% genes_target_present ~ "Hubs",
#         ensembl_gene_id %in% genes_target_absent  ~ "Pairs",
#         ensembl_gene_id %in% genes_random         ~ "Random",
#         ensembl_gene_id %in% genes_rest           ~ "All",
#         TRUE ~ NA_character_
#       )
#     ) %>%
#     filter(!is.na(group))

#   # 为了展示更稳健，建议用 log1p(TPM)
#   expr_in_target_samples <- expr_in_target_samples %>%
#     mutate(expr_log = log1p(expr))

#   p5  <- quantile(expr_in_target_samples$expr_log, 0.05, na.rm = TRUE)
#   p95 <- quantile(expr_in_target_samples$expr_log, 0.95, na.rm = TRUE)

#   expr_in_target_samples <- expr_in_target_samples %>%
#     mutate(
#       expr_log_w = pmin(pmax(expr_log, p5), p95)
#     )
#   ## =========================
#   ## 5) 画图：四组在 target_region 样本里的表达分布
#   ##    这里是把 “gene × sample” 的点作为分布（最贴近你原始 boxplot 逻辑）
#   ## =========================
#   expr_in_target_samples$group <- factor(
#     expr_in_target_samples$group,
#     levels = c(
#       "Hubs",
#       "Pairs",
#       "Random",
#       "All"
#     )
#   )

#   ## =========================
#   ## 6) （可选但强烈推荐）每个 sample 单独比较：facet by sample
#   ## =========================
#   grp_colors <- c("#DE1A58","#E97F4A", "#B7B7B7", "#1B211A")
#   expr_in_target_samples$sample <- gsub("H", "Human_", expr_in_target_samples$sample)
#   ref <- levels(expr_in_target_samples$group)[1]

#   comparisons_list <- lapply(
#     levels(expr_in_target_samples$group)[-1],
#     function(g) c(ref, g)
#   )

#   p2 <- ggplot(expr_in_target_samples, aes(x = group, y = expr_log_w, fill = group, color = group)) +
#     geom_boxplot(alpha = 0.8) +
#     facet_wrap(~ sample, ncol = 6, scales = "free_y") +
#     scale_fill_manual(values = grp_colors) +
#     scale_color_manual(values = grp_colors) +
#     theme_classic(base_family = "Helvetica") +
#     theme(axis.text.x = element_text(angle = 25, hjust = 1),
#     text = element_text(color = "black"),
#     axis.text = element_text(color = "black"),
#     axis.title = element_text(color = "black"),
#     strip.text = element_text(color = "black", face = "bold"),
#     legend.position = "none",
#     axis.line = element_line(color = "black"),
#     plot.title = element_text(hjust = 0.5)
#     ) +
#     labs(x = "", y = "log1p(TPM)", title = paste0(target_region))+
#     # coord_cartesian(ylim = c(quantile(expr_in_target_samples$expr_log, 0.05, na.rm = TRUE), quantile(expr_in_target_samples$expr_log, 0.95, na.rm = TRUE))) +
#     stat_compare_means(
#       comparisons = comparisons_list,
#       method = "wilcox.test",
#       label = "p.format",   # 或 "p.signif"
#       hide.ns = FALSE
#     )

#   ggsave(plot=p2, file = file.path(outdir, paste("brain.selected_region_expression",target_region,"pdf",sep=".")), width = 7, height = 4)
#   write.table(gene_expression[gene_expression$ensembl_gene_id %in% genes_target_present & !is.na(gene_expression$hgnc_symbol) & gene_expression$hgnc_symbol != "","hgnc_symbol",drop=FALSE], 
#     file = file.path(outdir, paste("brain.selected_region_expression",target_region,"txt",sep=".")),row.names=FALSE, quote=FALSE, col.names=FALSE, sep="\t")
# }


