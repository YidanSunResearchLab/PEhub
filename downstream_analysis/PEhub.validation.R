library(GenomicRanges)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(IRanges)
library(ggplot2)
library(rlang)
library(ggsci)
library(colorspace)


library(rtracklayer)
library(skitools)
library(chromunity)
library(arrow)
library(MASS)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(parallel)
library(Gviz)
library(GenomicInteractions)
library(trackViewer)
library(data.table)
library(GenomicRanges)
library(gUtils)
chr.names = c(1:22) #,"X" 
samplenames=c("Human_GM12878_NlaIII_reps")
indir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer"
outdir = gsub("ProcessedData", "Output", indir)
prewd = "/mnt/citadelb/syidan/Data/Processed/HumanEncodePorechip/Porec/IntegrativeAnalysis/"
# indir = "/home/syidan/syidan/Projects/SnakeHichip/snakehichip/data/example_data/"

##Main idea: compare HiChIP hub with porec data
##Hub level statistics:
#1. the percentage of hubs identified overlap Porec data
#2. the percentage of hub bin coverage overlap Porec data: overlapped bins / total bins in a hub
#3. overlapped porec reads multi-way concatemer distribution
#4. overlapped porec reads number distribution per hub
##Bin level statistics:
#1. the percentage of hub bins identified overlap Porec data 
#2. overlapped porec reads number distribution in all hubs
#3. overlapped porec reads number distribution per hub
#4. bin reads number correlation with weight_raw
#5. EE interaction overlap with porec data, the percerage of EE interaction supported by Porec data, and reads number distribution


#Overlap the hichip hub with porec data
build_porec_overlap_DT <- function(gr1,
                             gr2,
                             hubid_col = "promoter_id", 
                             min_reads=1) {

  # Fast overlap
  hits <- findOverlaps(gr1, gr2, ignore.strand = TRUE)
  qh <- queryHits(hits)     # index into gr1 (hub bins)
  sh <- subjectHits(hits)   # index into gr2 (PoRe-C fragments)

  # Statistics for original hub which do not have clusters
  DT <- data.table(
    hubid    = mcols(gr1)[[hubid_col]][qh],
    read_idx = mcols(gr2)[["read_idx"]][sh],
    peakid   = mcols(gr1)[["enhancer_id"]][qh],
    total_enhancer = mcols(gr1)[["n_enh"]][qh]
  )

  # same hubid and read_indx, the peak bins is merged. 
  DT <- DT[!is.na(total_enhancer)]
  setorder(DT, hubid, read_idx, peakid)
  DTu <- unique(DT[, .(hubid, read_idx, peakid, total_enhancer)])

  DT_reads_support <- DTu[
    ,
    .(
      peaks      = list(peakid),        # 已经唯一且有序
      total_enhancer = first(na.omit(total_enhancer)),
      k          = .N
    ),
    by = .(hubid, read_idx)
  ]
  DT_reads_support <- DT_reads_support[k >= 3]

  #unique the bins from all reads, This step identified hubs using porec data, which we take as true hubs, 
  # calculate each bin in each hub supported reads number.
  DT_hub_porec <- DT_reads_support[
    ,
    .(
      peak           = unique(unlist(peaks)),
      n_unique_peaks = uniqueN(unlist(peaks)),
      total_enhancer = as.numeric(unique(total_enhancer)),
      n_reads        = .N,
      k_mean         = mean(k),
      k_max          = max(k)
    ),
    by = hubid
  ][n_reads >= min_reads]
  
  # 每个 hub 被多少条 “k>=3” read 支撑,plot hub summary, how many hubs has >=3 supporting reads, max k how many, hub bin coverage percentage 
  hub_summary <- DT_reads_support[, .(
    total_enhancer = as.numeric(unique(total_enhancer)),
    n_peaks = uniqueN(unlist(peaks)),
    n_reads_support3 = .N,
    max_k = max(k)
  ), by = hubid][order(-n_reads_support3, -max_k)]

  # 统计每个 peak 被多少条 “k>=3” read 支撑, plot peak bin overlapped reads number summary, later combine with the weight info from gr1
  DT_support_peak <- unique(DT[read_idx %in% DT_reads_support$read_idx, .(hubid, peakid, read_idx)])
  peak_summary <- DT_support_peak[, .(n_reads = .N), by = .(hubid, peakid)][order(-n_reads)]

  return(list(
    DT = DT,  
    DT_reads_support = DT_reads_support,  
    DT_hub_porec = DT_hub_porec,
    hub_summary = hub_summary,  
    peak_summary = peak_summary  
  ))
}

# Function to plot histogram with mean line
plot_hist <- function(
    data,
    col,
    outdir,
    filename = "",
    max_value = NULL,
    binwidth = NULL,
    bins = NULL,
    fill = "grey80",
    color = "white"
  ) {
    col_quo <- enquo(col)

    df <- data %>%
      mutate(.x = as.numeric(!!col_quo)) %>%
      filter(!is.na(.x)) 

    if (!is.null(max_value)) {
      df <- df %>% mutate(.x = pmin(.x, max_value[2]))
    }

    med <- mean(df$.x)

    p = ggplot(df, aes(x = .x)) +
      geom_histogram(
        binwidth = binwidth,
        bins = bins,
        color = color,
        fill = fill
      ) +
      geom_vline(
        xintercept = med,
        linetype = "dashed",
        linewidth = 0.6
      ) +
      labs(x = "", y = "Count") +
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

    if (!is.null(max_value)) {
      p = p + coord_cartesian(xlim = c(max_value[1], max_value[2]+0.1))
    }


    ggsave(plot=p, file = file.path(outdir, paste("GM12878.porec_validation.histogram",filename,"pdf",sep=".")), width = 2.5, height = 2.5)
  }


#hub level statistics
build.hub.summary <- function(observed_hubs_per_promoter_ori, porec_raw, i, index, outdir, min_reads=1){
  observed_hubs_per_promoter_ori$distancekb = observed_hubs_per_promoter_ori$D / 1000
  hub_summary_filt <- porec_raw$hub_summary[n_reads_support3 >= min_reads]
  hub_summary_filt$enhancer_bin_coverage <- hub_summary_filt$n_peaks/hub_summary_filt$total_enhancer * 100
  hub_summary_filt$n_reads_support3_hubsize <- hub_summary_filt$n_reads_support3 / hub_summary_filt$n_peaks
  print(paste("Percentage of hubs supported by PoreC data:", nrow(hub_summary_filt)/length(unique(observed_hubs_per_promoter_ori$promoter_id))))

  print("Percentage of hubs bin coverage supported by PoreC data using hub_summary: bins supported by porec / all bins in a hub") 
  print(summary(hub_summary_filt$n_peaks/hub_summary_filt$total_enhancer))

  print("Percentage of hubs bin coverage supported by PoreC data using porecoverlap summary reads>=min_reads: bins supported by porec / all bins in a hub") 
  print(summary(porec_raw$DT_hub_porec$n_unique_peaks/porec_raw$DT_hub_porec$total_enhancer))

  print("Summary of k multi-way for per hub supported by PoreC data") 
  print(summary(hub_summary_filt$max_k))

  print("Summary of reads number per hub supported by PoreC data using hub_summary: bins supported by porec / all bins in a hub") 
  print(summary(hub_summary_filt$n_reads_support3))
  
  select.color <- ifelse(index == "pvalue", "#FF0087", ifelse(index == "hub1", "#00F7FF", "grey70"))
  plot_hist(hub_summary_filt, n_reads_support3, outdir = outdir, filename = paste("hub_n_reads_support3", i, index, sep="."), max_value = c(0, 200), binwidth = 10, fill = select.color) #, xlim = c(0, 500)
  plot_hist(hub_summary_filt, n_reads_support3_hubsize, outdir = outdir, filename = paste("hub_n_reads_support3_hubsize", i, index, sep="."), max_value = c(0, 5), binwidth = 0.25, fill = select.color) #, xlim = c(0, 500)
  plot_hist(hub_summary_filt, max_k, outdir = outdir, filename = paste("hub_max_k", i, index, sep="."), max_value = c(3, 10), binwidth = 1, fill = select.color) #, xlim = c(0, 500)
  plot_hist(hub_summary_filt, total_enhancer, outdir = outdir, filename = paste("hub_total_enhancer", i, index, sep="."), max_value = c(3, 30), binwidth = 1, fill = select.color) #, xlim = c(0, 500)
  plot_hist(hub_summary_filt, n_peaks, outdir = outdir, filename = paste("hub_n_peaks", i, index, sep="."), max_value = c(4, 30), binwidth = 1, fill = select.color) #, xlim = c(0, 500)
  plot_hist(hub_summary_filt, enhancer_bin_coverage, outdir = outdir, filename = paste("hub_enhancer_bin_coverage", i, index, sep="."), binwidth = 10, fill = select.color) #, xlim = c(0, 500)
  plot_hist(observed_hubs_per_promoter_ori, distancekb, outdir = outdir, filename = paste("hub_distance", i, index, sep="."), max_value = c(0, 1000), binwidth = 20, fill = select.color) #, xlim = c(0, 500)
  plot_hist(porec_raw$peak_summary, n_reads, outdir = outdir, filename = paste("bin_n_reads_support3", i, index, sep="."), max_value = c(0, 100), binwidth = 5, fill = select.color) #, xlim = c(0, 500)
  
  summary_df <- data.frame(number=c(
    hub_percentage=nrow(hub_summary_filt)/length(unique(observed_hubs_per_promoter_ori$promoter_id)),
    summary(hub_summary_filt$n_peaks/hub_summary_filt$total_enhancer)
  ),
  weight_method = i,
  index = index,
  stat = c(
    "hub_percentage",
    "Minimum",
    "1st Quantile",
    "Median",
    "Mean",
    "3rd Quantile",
    "Maximum"
  )
  )

  return(summary_df)
}

plot_weight_method_bar <- function(df, outdir = outdir, filename = "barplot", index_value = "pvalue", bar_stat = "hub_percentage") {
  wm_order <- df %>%
    filter(index == "pvalue", stat == "hub_percentage") %>%
    arrange(number) %>%
    pull(weight_method)
  
  dat <- df %>%
    filter(index == index_value, stat == bar_stat) %>%
    mutate(weight_method = factor(weight_method, levels = wm_order),
            rank_val = as.numeric(weight_method),
            number = as.numeric(number) * 100)                 # 用于渐变映射（按排序）

  

  p <- ggplot(dat, aes(x = weight_method, y = number, fill = weight_method)) +
    geom_col(width = 0.8, color = NA, alpha = 0.9) +
    geom_text(
      aes(label = sprintf("%.0f", number)),
      vjust = -0.3,
      size = 2.8,
      color = "black"
    ) +
    # scale_fill_gradientn(colors = c("grey80", "orange", "red3")) +
    scale_fill_d3("category20") +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "weight_method", y = "Percentage (%)") +
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
    ggsave(plot=p, file = file.path(outdir, paste("GM12878.porec_validation.methods_comparison",filename,index_value,"pdf",sep=".")), width = 3.3, height = 3)
}

plot_weight_method_summary_line <- function(df, outdir = outdir, filename = "lineplot", index_value = "pvalue", bar_stat = "hub_percentage") {

  base <- df %>%
    filter(index == index_value)
  # 先拿“排序规则”：按 order_stat 的 number 对 weight_method 排序
  wm_order <- df %>%
    filter(index == "pvalue", stat == "hub_percentage") %>%
    arrange(number) %>%
    pull(weight_method)

  # 需要画线的统计量顺序
  stat_order <- c("Minimum", "1st Quantile", "Median", "3rd Quantile", "Maximum")

  dat <- base %>%
    filter(stat %in% stat_order) %>%
    mutate(
      stat = factor(stat, levels = stat_order),
      weight_method = factor(weight_method, levels = wm_order),
      rank_val = as.numeric(weight_method),  # 渐变映射（沿用 barplot 排序）
      number = as.numeric(number) * 100
    )

  cols <- setNames(
    colorRampPalette(c("grey80", "orange", "red3"))(length(wm_order)),
    wm_order
    )

  p=ggplot(dat, aes(x = stat, y = number, group = weight_method, color = weight_method)) +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    geom_point(size = 1.8, alpha = 0.7) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_color_d3("category20") +
    # scale_color_manual(values = cols) +
    labs(x = NULL, y = "Percentage (%)") +
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
  ggsave(plot=p, file = file.path(outdir, paste("GM12878.porec_validation.methods_comparison",filename,index_value,"pdf",sep=".")), width = 2.7, height = 2.5)
}

##Porec plot
build.porec.plot <- function(gr2, outdir = outdir) {
  ################################
  # Use GenomicInteractions (Gviz) to visualize the genomic regions for all the binned concatemers (chrommunity)
  ################################
  library(GenomicInteractions)
  dir.create(file.path(outdir, "examples"))

  load(file.path(indir,paste("multiple_result.exampleGM12878.porec_raw", "log_minmax", "bin_log_ratio_sig", "pvalue", "RData", sep=".")))
  observed_hubs_per_promoter_sub_pvalue <- observed_hubs_per_promoter_sub %>% filter(!is.na(hub_p_adj_global), hub_p_adj_global <= 0.05, reproducibility_rate >= 0.5)
  gr1 = makeGRangesFromDataFrame(observed_hubs_per_promoter_sub_pvalue, keep.extra.columns = TRUE) 
  hits <- findOverlaps(gr1, gr2, ignore.strand = TRUE)
  gr2.overlap = gr2[unique(subjectHits(hits))]
  gr2.overlap = gr2.overlap[gr2.overlap$read_idx %in% porec_raw$DT_reads_support$read_idx]


  # 1. Split the GRanges object by 'cid' (handles multiple groups, here only one for simplicity)
  for (peak_name in unique(gr1$promoter_id) ){
      # peak_name = "peak_chr1_100034_100040"
      gr = gr1[gr1$promoter_id == peak_name]
      hits <- findOverlaps(gr, gr2.overlap, ignore.strand = TRUE)
      gr.reads = gr2.overlap[unique(subjectHits(hits))]
      chr=unique(as.character(seqnames(gr)))
      
      # 3. Create the GenomicInteractions object from the anchors
      interactions <- GenomicInteractions(rep(granges(gr[gr$pair_type=="PP"]), length(gr)-1), granges(gr[gr$pair_type=="PE"]), counts = log1p(gr[gr$pair_type=="PE"]$counts))
      interaction_track <- InteractionTrack(interactions, name = "HiChIP", chromosome = chr )
      displayPars(interaction_track) = list(col.interactions="#D91656", 
                                            col.anchors.fill ="black",
                                            col.anchors.line = "black",
                                            anchor.height = 0.2)
      

      
      ################################
      # Use Gviz to visualize the genomic regions for all the reads
      ################################
      # Create the AnnotationTrack
      gr.reads <- gr.reads[seq_len(min(1000, length(gr.reads)))]
      track <- AnnotationTrack(range = gr.reads, 
                              group = gr.reads$read_idx, # Group the ranges by cid
                              type = "lollipop",
                              name = paste0(chr,":", min(start(gr.reads)), "-", max(end(gr.reads))), 
                              genome = "hg38",
                              col.line = "grey",
                              col = "#640D5F",           # Line color
                              fill = "#640D5F",          # Box color
                              chromosome = chr)
        
      # Create an IdeogramTrack for chromosome visualization (optional)
      ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = chr)

      
      # Plot both the GeneRegionTrack and the AnnotationTrack
      pdf(file.path(outdir, "examples", paste("GM12878",peak_name,"track.porec.plot.pdf", sep="-")),width = 10, height = 20)
      track_heights <- c(0.2, 1, 5)  # Adjust heights as needed for each track
      plotTracks(list(ideoTrack, interaction_track, track),from = min(start(gr)), to = max(end(gr)), sizes = track_heights) #
      dev.off()
  }

}


##the porec data
#load(file.path(indir,"multiple_result_exampleGM12878.prephichip.RData"))
raw_reads_granges = readRDS(file.path(prewd, "TestData", paste0(samplenames[1],".contacts.nlaiii.granges.rds") ))
raw_reads_granges <- keepSeqlevels(raw_reads_granges, paste0("chr",chr.names), pruning.mode="coarse")
gr2 = raw_reads_granges

# ##Start processing each method
# for (method in c("log_minmax")){ #, "log_maxnorm", "log_zscore"
#   print(paste0("Starting method:", method))
#   # for (i in c("bin_percentile_plus_sig", "bin_log_ratio_sig", "bin_diff_global", "bin_diff_binmax")) {
#   for (i in c("distance_only", "count_only", "sig_only", "count_sig", "zscore_residual", "count_sig_plus_dist_linear", "bin_percentile_plus_sig", "bin_log_ratio_sig", "bin_diff_global", "bin_diff_binmax", "bin_percentile")) {
#     print(paste0("Starting method:", i))
#     load(file.path(indir,paste("multiple_result.exampleGM12878.hub", method, i, "preprocess.RData",sep=".")))
#     observed_hubs_per_promoter_ori <- observed_hubs_per_promoter$hubs
#     for (index in c("pvalue")){ #"all","hub1","others",
#       print(paste0("Starting method:", index))
#       if(index=="all"){
#         observed_hubs_per_promoter_sub <- observed_hubs_per_promoter_ori
#       } else if(index=="hub1"){
#         observed_hubs_per_promoter_sub <- observed_hubs_per_promoter_ori[!is.na(observed_hubs_per_promoter_ori$hub_index) & observed_hubs_per_promoter_ori$hub_index=="hub1",]
#       } else if(index=="others"){
#         observed_hubs_per_promoter_sub <- observed_hubs_per_promoter_ori[!is.na(observed_hubs_per_promoter_ori$hub_index) & observed_hubs_per_promoter_ori$hub_index!="hub1",]
#       } else if(index=="pvalue"){
#         load(file.path(indir,paste("multiple_result.exampleGM12878.hub", method, i, "RData",sep=".")))
#         observed_hubs_per_promoter_sub <- observed_hubs_per_promoter_sub_hub %>% filter(!is.na(hub_p_adj_global), hub_p_adj_global <= 0.05, reproducibility_rate >= 0.5) #hub_p_value_global
#       }
#       gr1 = makeGRangesFromDataFrame(observed_hubs_per_promoter_sub, keep.extra.columns = TRUE) 
#       gr1$is_promoter = gr1$pair_type == "PP"

#       porec_raw <- build_porec_overlap_DT(gr1 = gr1, gr2 = gr2)
#       save(observed_hubs_per_promoter_sub, porec_raw, file = file.path(indir,paste("multiple_result.exampleGM12878.porec_raw", method, i, index, "RData", sep=".")))
#     }}}

# ##Start calculating each method 
# for (method in c("log_minmax")){ #, "log_maxnorm", "log_zscore"
#   hub_level_summary = data.frame()
#   # for (i in c("bin_log_ratio_sig")) {
#   for (i in c("distance_only", "count_only", "sig_only", "count_sig", "zscore_residual", "count_sig_plus_dist_linear", "bin_percentile_plus_sig", "bin_log_ratio_sig", "bin_diff_global", "bin_diff_binmax", "bin_percentile")) {
#     for (index in c("hub1","others","pvalue")){ #"all"
#       load(file.path(indir,paste("multiple_result.exampleGM12878.porec_raw", method, i, index, "RData", sep=".")))
#       print(paste("Starting method:", method, i, index, sep=";"))
#       hub_level_summary <- rbind(hub_level_summary, build.hub.summary(observed_hubs_per_promoter_sub, porec_raw, i, index, outdir, min_reads=1))
#     }}
#   write.table(hub_level_summary, file = file.path(outdir, paste("GM12878.porec_validation.hub_level_summary", method, "tsv", sep=".")), sep="\t", quote=FALSE, row.names=FALSE)
#   for (index in c("hub1","others","pvalue")){ #"all"
#     plot_weight_method_bar(df=hub_level_summary, outdir =  outdir, filename = "barplot", index_value = index, bar_stat = "hub_percentage")
#     plot_weight_method_summary_line(df=hub_level_summary, outdir =  outdir, filename = "lineplot", index_value = index, bar_stat = "hub_percentage")
#   }
# }


##start porec plot
build.porec.plot(gr2, outdir=outdir)








##Not used anymore
# for (method in c("log_minmax", "log_maxnorm", "log_zscore")){
#   # for (i in c("bin_percentile_plus_sig", "bin_log_ratio_sig", "bin_diff_global", "bin_diff_binmax")) {
#   for (i in c("distance_only", "count_only", "sig_only", "count_sig", "zscore_residual", "count_sig_plus_dist_linear", "bin_percentile_plus_sig", "bin_log_ratio_sig", "bin_diff_global", "bin_diff_binmax", "bin_percentile")) {
#     for (index in c("hub1","others","pvalue")){ #"all"
#       load(file.path(indir,paste("multiple_result.exampleGM12878.porec_raw", method, i, index, "RData", sep=".")))
#       print(paste0("Starting method:", method, i))
#       bin_level_summary <- build.bin.level(observed_hubs_per_promoter_sub, porec_raw)
#     }}
# }

# for (method in c("log_minmax", "log_maxnorm", "log_zscore")){
#   # for (i in c("bin_percentile_plus_sig", "bin_log_ratio_sig", "bin_diff_global", "bin_diff_binmax")) {
#   for (i in c("distance_only", "count_only", "sig_only", "count_sig", "zscore_residual", "count_sig_plus_dist_linear", "bin_percentile_plus_sig", "bin_log_ratio_sig", "bin_diff_global", "bin_diff_binmax", "bin_percentile")) {
#     for (index in c("all","hub1","others")){
#       load(file.path(indir,paste("multiple_result.exampleGM12878.hub", method, i,  "RData", sep=".")))
#       load(file.path(indir,paste("multiple_result.exampleGM12878.porec_raw", method, i, index, "RData", sep=".")))
#       ee_level_summary <- build.ee.level(prep_hichip_all, observed_hubs_per_promoter_sub, porec_raw)
#     }}
#   }





#bin level statistics
build.bin.level <- function(observed_hubs_per_promoter_ori, porec_raw){
  ########################################################################
  #We separate HiChIP hub bins into porec overlap and nonoverlap to check how many bins were overlapped by porec data
  DT_matched <- porec_raw$DT_hub_porec %>%
    left_join(
      observed_hubs_per_promoter_ori,
      by = c(
        "hubid" = "promoter_id",
        "peak"  = "enhancer_id"
      )
    )

  DT_unmatched <- observed_hubs_per_promoter_ori %>%
    anti_join(
      porec_raw$DT_hub_porec,
      by = c(
        "promoter_id" = "hubid",
        "enhancer_id" = "peak"
      )
    ) %>%
    mutate(
      peak = enhancer_id,
      hubid = promoter_id
    )

  print("How many Porec identified hubs defined in the HiChIP identified hubs, should be really high")
  print(table(DT_matched$hub_index))
  print("How many Porec unsupported hubs defined in the HiChIP dataset, should be really low")
  print(table(DT_unmatched$hub_index))

  ########################################################################
  #each bin overlapped reads number summary
  observed_hubs_per_promoter_merged <-
    observed_hubs_per_promoter_ori %>%
    left_join(
      porec_raw$peak_summary,
      by = c(
        "promoter_id" = "hubid",
        "enhancer_id" = "peakid"
      )
    )

  hub_reads_summary <- observed_hubs_per_promoter_merged %>%
  group_by(hub_index) %>%
  summarise(
    n_hubs = n(),
    total_n_reads = sum(n_reads, na.rm = TRUE),
    min_n_reads   = min(n_reads, na.rm = TRUE),
    q25_n_reads   = quantile(n_reads, 0.25, na.rm = TRUE),
    median_n_reads = median(n_reads, na.rm = TRUE),
    q75_n_reads   = quantile(n_reads, 0.75, na.rm = TRUE),
    max_n_reads   = max(n_reads, na.rm = TRUE),
    mean_n_reads = mean(n_reads, na.rm = TRUE),
    iqr_n_reads   = IQR(n_reads, na.rm = TRUE),
    # min_n_reads   = min(n_reads, na.rm = TRUE),
    # q25_n_reads   = quantile(n_reads, 0.25, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(as.integer(str_remove(hub_index, "hub")))
  print("In sub hub level, bin reads statistics")
  print(hub_reads_summary, width=Inf)

  #each bin overlapped reads number summary (in each hub)
  hub_reads_summary_eachhub <- observed_hubs_per_promoter_merged %>%
  group_by(hub_id) %>%
  summarise(
    n_hubs = n(),
    total_n_reads = sum(n_reads, na.rm = TRUE),
    min_n_reads   = min(n_reads, na.rm = TRUE),
    q25_n_reads   = quantile(n_reads, 0.25, na.rm = TRUE),
    mean_n_reads = mean(n_reads, na.rm = TRUE),
    median_n_reads = median(n_reads, na.rm = TRUE),
    q75_n_reads   = quantile(n_reads, 0.75, na.rm = TRUE),
    max_n_reads   = max(n_reads, na.rm = TRUE),
    iqr_n_reads   = IQR(n_reads, na.rm = TRUE),
    # min_n_reads   = min(n_reads, na.rm = TRUE),
    # q25_n_reads   = quantile(n_reads, 0.25, na.rm = TRUE),
    .groups = "drop"
  )
  print("In each sub hub level, bin reads statistics")
  print(hub_reads_summary_eachhub, width=Inf)  

  ########################################################################
  #supporting reads per bin correlation with weight_raw
  print("Single peak bin level reads statistics, not including bin cooccurrences")
  print(paste("method", i, "weight_raw and porec overlapped reads per bin correlation pearson and spearman: \n", 
        cor(observed_hubs_per_promoter_merged$weight_raw.x, observed_hubs_per_promoter_merged$n_reads, method = "pearson", use = "complete.obs"), "\n",
        cor(observed_hubs_per_promoter_merged$weight_raw.x, observed_hubs_per_promoter_merged$n_reads, method = "spearman", use = "complete.obs"), "\n") )

  list(
    DT_matched = DT_matched, 
    DT_unmatched = DT_unmatched, 
    observed_hubs_per_promoter_merged = observed_hubs_per_promoter_merged,
    hub_reads_summary_eachhub = hub_reads_summary_eachhub  
  )
}

build.ee.level <- function(prep_hichip_all, observed_hubs_per_promoter_ori, porec_raw){
  ########################################################################
  #Check if Enhancer enhancer interactions in Hichip is supported by porec data, so we know how informative they are. 
  EE = prep_hichip_all$loops[prep_hichip_all$loops$pair_type=="EE",]

  DT_matched <- porec_raw$DT_hub_porec %>%
    left_join(
      observed_hubs_per_promoter_ori,
      by = c(
        "hubid" = "promoter_id",
        "peak"  = "enhancer_id"
      )
    )

  DT_unmatched <- observed_hubs_per_promoter_ori %>%
    anti_join(
      porec_raw$DT_hub_porec,
      by = c(
        "promoter_id" = "hubid",
        "enhancer_id" = "peak"
      )
    ) %>%
    mutate(
      peak = enhancer_id,
      hubid = promoter_id
    )

  print("How many Porec identified hubs are supported by HiChIP enhancer enhancer interactions, should be really high true")
  print(table(DT_matched$peak %in% c(EE$promoter_id, EE$enhancer_id)))
  print("How many Porec unsupported hubs are supported by HiChIP enhancer enhancer interactions, should be really high false")
  print(table(DT_unmatched$enhancer_id %in% c(EE$promoter_id, EE$enhancer_id)))
  
  setDT(DT_matched); setDT(EE); setDT(DT_unmatched)
  EE[, `:=`(
    a = pmin(promoter_id, enhancer_id),
    b = pmax(promoter_id, enhancer_id)
  )]
  #matched
  ans <- DT_matched[, .(peak=unique(peak)), by=hubid][ , {p=peak; if(length(p)<2) return(NULL); m=t(combn(p,2)); .(a=pmin(m[,1],m[,2]), b=pmax(m[,1],m[,2]))}, by=hubid][ EE, on=.(a,b), nomatch=0 ][ , .(n_EE=.N, sum_counts=sum(counts), mean_counts=mean(counts), max_counts=max(counts)), by=hubid ]
  tot <- DT_matched[, .(peak=unique(peak)), by=hubid][, .(n_total=choose(.N,2)), by=hubid]
  EE_DT_matched <- tot[ans, on="hubid"][, frac := n_EE/n_total]
  #unmatched
  ans_un <- DT_unmatched[,.(peak=unique(peak)),by=hubid][,{p<-peak;if(length(p)<2)NULL else {m<-t(combn(p,2));.(a=pmin(m[,1],m[,2]),b=pmax(m[,1],m[,2]))}},by=hubid][EE,on=.(a,b),nomatch=0][,.(n_EE=.N,sum_counts=sum(counts),mean_counts=mean(counts),max_counts=max(counts)),by=hubid]
  tot_un <- DT_unmatched[, .(peak=unique(peak)), by=hubid][, .(n_total=choose(.N,2)), by=hubid]
  EE_DT_unmatched <- tot_un[ans_un, on="hubid"][, frac := n_EE/n_total]

  print("How many Porec identified hubs are supported by HiChIP enhancer enhancer interactions, supported number / total number, should be really high")
  print(summary(EE_DT_matched$frac))
  print("How many Porec identified hubs are supported by HiChIP enhancer enhancer interactions, supported EE counts number, should be really high")
  print(summary(EE_DT_matched$sum_counts))

  print("How many Porec unsupported hubs are supported by HiChIP enhancer enhancer interactions, supported number / total number, should be really low")
  print(summary(EE_DT_unmatched$frac))
  print("How many Porec unsupported hubs are supported by HiChIP enhancer enhancer interactions, supported EE counts number, should be really low")
  print(summary(EE_DT_unmatched$sum_counts))

  list(
    DT_matched = DT_matched, 
    DT_unmatched = DT_unmatched, 
    EE = EE,
    EE_DT_matched = EE_DT_matched, 
    EE_DT_unmatched = EE_DT_unmatched  
  )
}







##calling chromunity with enhancer and promoter regions (ajusted: resolution=5k, training/testing, k=5)
# if(TRUE){
# print("Calling chromunity with ATAC peaks")
# for (samplename in samplenames[1]){
#   for (method in c("atac")){ #"h3kall","his"
#   print(paste("Processing sample:", samplename))
#   set.seed(198)
#   resolution = 2.5e4 #2.5e4 ##
#   k=7
#   genome_length= read.delim("~/syidan/Genomes/GRCh38/release-47-index/genome_fasta/genome.fa.fai",  header = FALSE, stringsAsFactors = FALSE)
#   load(file.path(prewd, "ProcessedData",paste0("chromunity.covariats.",samplename,".RData") ))
#   raw_reads_granges = readRDS(file.path(prewd, "TestData", paste0(samplename,".contacts.nlaiii.granges.rds") ))
#   raw_reads_granges <- keepSeqlevels(raw_reads_granges, paste0("chr",chr.names), pruning.mode="coarse")

#   ## getting E-P annotations
#   # enhancer_region = import(file.path(upperdir, "externaldata", paste(samplename,method,"bed",sep="."))) #result include atacresult
#   enhancer_region = import("~/syidan/Data/Processed/HiCHIP_GSE/ATAC/merged.snakePipes.out/MACS2/Human_GM12878_unknown_WT_unknown_ATAC_standard_mergedSRR.filtered.short.BAM_summits.bed") #result include atacresult
#   enhancer_region <- keepSeqlevels(enhancer_region, paste0("chr",chr.names), pruning.mode="coarse")

#   for (chr in paste0("chr",c(0))){ #,chr.names
#     print(paste("Processing chromosome:", chr))
#     if (chr == "chr0"){
#       enhancer_region_chr = enhancer_region
#     } else {
#       enhancer_region_chr=enhancer_region[seqnames(enhancer_region) %in% chr]
#     }

#     ## Targets
#     enhancer_region_chr_extend = gr.reduce(enhancer_region_chr+resolution)

#     ## Subsample
#     chr_parq_gr = raw_reads_granges
#     #chr_parq_gr = raw_reads_granges %&% enhancer_region_chr_extend
#     raw_reads_chromunity_all = re_chromunity(concatemers = chr_parq_gr, windows = enhancer_region_chr_extend, piecewise = FALSE, shave = TRUE, resolution = resolution, mc.cores = 50)
#     print(paste("All binsets identified:", length(unique(raw_reads_chromunity$binsets$chid))))
#     save(raw_reads_chromunity_all, file=file.path(indir ,paste(samplename,chr,"chromunity",paste0(method,"results"),"original.RData",sep=".")))
#     }
#   }
#   rm(list=setdiff(ls(), basic.objects))
#   gc()
# }
# }


#' Match HiChIP hubs to Chromunity synergistic interactions
#'
#' @param gr1     GRanges: each row is a locus in a HiChIP hub
#'                      mcols must contain: hub_id (char), is_promoter (logical, or 0/1)
#' @param gr2 GRanges: each row is a locus in a Chromunity synergy set
#'                      mcols must contain: synergy_id (char)
#' @param hub_id_col    Column name in gr1 mcols for hub_id
#' @param syn_id_col    Column name in gr2 mcols for synergy_id
#' @param promoter_col  Column name in gr1 mcols, logical, marking promoter loci
#' @param min_shared    Minimum number of shared loci between hub and synergy
#' @param min_cov_hub   Minimum fraction of hub loci covered by synergy (n_shared / size_hub)
#' @param min_jaccard   Minimum Jaccard index between hub loci and synergy loci
#'
#' @return A list with:
#'   - matches: data.frame summarizing matched hub–synergy pairs and overlap metrics
#'   - overlaps_raw: data.frame of per-locus overlaps (hub_id, synergy_id, locus indices)
#'
#' 
# match_hubs_chromunity <- function(
#   gr1,
#   gr2,
#   hub_id_col   = "hubid",
#   syn_id_col   = "chid",
#   promoter_col = "is_promoter",
#   min_shared   = 3,
#   min_cov_hub  = 0.5,
#   min_jaccard  = 0.3
# ) {
#   # ---- 0. 基本检查 ----
#   stopifnot(
#     inherits(gr1, "GRanges"),
#     inherits(gr2, "GRanges")
#   )
  
#   hits <- findOverlaps(gr1, gr2, ignore.strand = TRUE)
#   overlaps_gr1 <- gr1[unique(queryHits(hits)), c("D","counts", "qvalue","hubid","is_promoter", "peakoverlap1","peakoverlap2")]
#   overlaps_gr2 <- gr1[-unique(queryHits(hits)), c("D","counts", "qvalue","hubid","is_promoter", "peakoverlap1","peakoverlap2")]

#   # For overlaps_gr1
#   gr1_df <- as.data.frame(overlaps_gr1) %>%
#     tibble::as_tibble() %>%
#     mutate(
#       seqnames = as.character(seqnames),
#       strand = as.character(strand)
#     )

#   write_tsv(gr1_df, "/home/syidan/syidan/Projects/SnakeHichip/snakehichip/scripts/overlaps_gr1.tsv")

#   # For overlaps_gr2
#   gr2_df <- as.data.frame(overlaps_gr2) %>%
#     tibble::as_tibble() %>%
#     mutate(
#       seqnames = as.character(seqnames),
#       strand = as.character(strand)
#     )

#   write_tsv(gr2_df, "/home/syidan/syidan/Projects/SnakeHichip/snakehichip/scripts/overlaps_gr2.tsv")

#   h_mcols <- mcols(gr1)
#   c_mcols <- mcols(gr2)
  
#   if (!all(c(hub_id_col, promoter_col) %in% colnames(h_mcols))) {
#     stop("gr1 mcols must contain columns: ",
#          hub_id_col, " and ", promoter_col)
#   }
#   if (!syn_id_col %in% colnames(c_mcols)) {
#     stop("gr2 mcols must contain column: ", syn_id_col)
#   }
  
#   # ---- 1. 确保有 locus_id（如果没有就自动生成） ----
#   if (!"locus_id" %in% colnames(h_mcols)) {
#     mcols(gr1)$locus_id <- paste0("H_", seq_along(gr1))
#   }
#   if (!"locus_id" %in% colnames(c_mcols)) {
#     mcols(gr2)$locus_id <- paste0("C_", seq_along(gr2))
#   }
  
#   # 为了简化后面代码，拎出列名
#   h_hub_id   <- h_mcols[[hub_id_col]]
#   h_promoter <- as.logical(h_mcols[[promoter_col]])
#   h_locus_id <- mcols(gr1)[["locus_id"]]
  
#   c_syn_id   <- c_mcols[[syn_id_col]]
#   c_locus_id <- mcols(gr2)[["locus_id"]]
  
#   # ---- 2. 用 findOverlaps 在 locus 层面建立 Hub–Synergy 的连接 ----
#   seqinfo(gr1) <- seqinfo(gr2)
  
#   if (length(hits) == 0L) {
#     warning("No overlaps between gr1 and gr2.")
#     return(
#       list(
#         matches      = tibble(),
#         overlaps_raw = tibble()
#       )
#     )
#   }
  
#   overlaps_raw <- tibble(
#     h_query_idx = queryHits(hits),
#     c_subj_idx  = subjectHits(hits),
    
#     hub_id      = h_hub_id[queryHits(hits)],
#     synergy_id  = c_syn_id[subjectHits(hits)],
    
#     h_locus_id  = h_locus_id[queryHits(hits)],
#     c_locus_id  = c_locus_id[subjectHits(hits)],
    
#     is_promoter = h_promoter[queryHits(hits)],
#     counts = h_mcols[["counts"]][queryHits(hits)],
#     qvalue = h_mcols[["qvalue"]][queryHits(hits)]
    
#   )
  
  
#   # ---- 3. 计算每个 hub 与每个 synergy 的交集大小等 ----
#   # hub 总大小
#   hub_size <- tibble(
#     hub_id    = h_hub_id,
#     h_locus_id = h_locus_id
#   ) %>%
#     distinct() %>%
#     count(hub_id, name = "size_hub")
  
#   # synergy 总大小
#   syn_size <- tibble(
#     synergy_id = c_syn_id,
#     c_locus_id = c_locus_id
#   ) %>%
#     distinct() %>%
#     count(synergy_id, name = "size_syn")
  
#   # 汇总交集
#   overlap_summary <- overlaps_raw %>%
#     group_by(hub_id, synergy_id) %>%
#     summarise(
#       n_shared       = n_distinct(h_locus_id),
#       n_shared_prom  = n_distinct(h_locus_id[is_promoter]),
#       .groups = "drop"
#     ) %>%
#     left_join(hub_size, by = "hub_id") %>%
#     left_join(syn_size, by = "synergy_id") %>%
#     mutate(
#       jaccard = n_shared / (size_hub + size_syn - n_shared),
#       cov_hub = n_shared / size_hub,
#       cov_syn = n_shared / size_syn，
#       relative_enrichment = (n_shared / size_hub) /   # 比 random 高多少倍
#           (length(gr2) * mean(width(gr2)) / sum(as.numeric(seqlengths(gr1))))
#     )
  
#   # ---- 4. 按 promoter + 多 enhancer 的逻辑过滤匹配 ----
#   matches <- overlap_summary %>%
#     # 至少 promoter 一致（同一 promoter-centered hub）
#     filter(n_shared_prom >= 1) %>%
#     # 至少共享 min_shared 个 loci（一般是 promoter + ≥2 enhancers）
#     filter(n_shared >= min_shared) %>%
#     # 至少覆盖一定比例的 hub
#     filter(cov_hub >= min_cov_hub) %>%
#     # Jaccard 相似度也要过一个阈值
#     filter(jaccard >= min_jaccard) %>%
#     arrange(desc(jaccard), desc(cov_hub), desc(n_shared))
  
#   return(
#     list(
#       matches      = matches,
#       overlaps_raw = overlaps_raw
#     )
#   )
# }



match_hubs_readsgranges <- function(
  gr1,
  gr2,
  hub_id_col   = "hubid",
  syn_id_col   = "chid",
  promoter_col = "is_promoter",
  min_shared   = 3,
  min_cov_hub  = 0.5,
  min_jaccard  = 0.3
) {
  # hub_table: data.frame with columns hub_id and locus_id
  # multiway_contacts: list of character vectors, each vector is the locus_ids hit by one concatemer

  
  hits <- findOverlaps(gr1, gr2, ignore.strand = TRUE)
  overlaps_gr1 <- gr1[unique(queryHits(hits)), c("D","counts", "qvalue","hubid","is_promoter", "peakoverlap1","peakoverlap2")]
  overlaps_gr2 <- gr1[-unique(queryHits(hits)), c("D","counts", "qvalue","hubid","is_promoter", "peakoverlap1","peakoverlap2")]

  qh <- queryHits(hits)     # index into porec_chr
  sh <- subjectHits(hits)   # index into hub_chr
  hubid <- mcols(gr1)$hubid[qh]
  peakid <- mcols(gr1)$enhancer_id[qh]
  peakid[is.na(peakid)] <- mcols(gr1)$promoter_id[qh][is.na(peakid)]
  rid   <- mcols(gr2)$read_idx[sh]

  DT <- data.table(read_idx = rid, hubid = hubid, peakid = peakid, sh = sh)
  # k = 该 read 在该 hub 中 overlap 的不同 hub member 数
  supp <- DT[, .(peaks = list(sort(unique(peakid))),k = uniqueN(peakid)), by = .(hubid, read_idx)]
  support3 <- supp[k >= 3]
  
  hub_summary <- support3[, .(
    n_reads_support3 = .N,
    max_k = max(k)
  ), by = hubid][order(-n_reads_support3, -max_k)]


  # 统计每个 peak 被多少条 “k>=3” read 支撑
  reads_supporting <- support3$read_idx
  DT_support <- unique(DT[read_idx %in% reads_supporting, .(read_idx, peakid)])
  peak_support_counts <- DT_support[, .(n_reads = .N), by = peakid][order(-n_reads)]


  # 从 gr1 抽取坐标（去重）
  top_peaks <- peak_support_counts[n_reads >= 1, peakid]   # 举例：至少被 5 条 multi-way read 支撑
  mcols(gr1)$enhancer_id[is.na(mcols(gr1)$enhancer_id)] <- mcols(gr1)$promoter_id[is.na(mcols(gr1)$enhancer_id)]
  hub_bins_supported <- gr1[gr1$enhancer_id %in% top_peaks]
  hub_bins_supported_not <- gr1[!gr1$enhancer_id %in% top_peaks]
  
  return(hub_support_df)
}





#' Compute Jaccard distance between two sets (character vectors)
jaccard_vec <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  inter <- length(intersect(a, b))
  uni   <- length(union(a, b))
  if (uni == 0) return(NA_real_)
  inter / uni
}

compute_differential_hubs <- function(
  gr1,
  gr2,
  hub_id_col   = "chid",
  promoter_col = "is_promoter",
  locus_id_col = "locus_id",
  n_perm       = 1000,
  min_size     = 3
) {
  stopifnot(inherits(gr1, "GRanges"), inherits(gr2, "GRanges"))
  
  m1 <- mcols(gr1)
  m2 <- mcols(gr2)
  
  if (!hub_id_col %in% colnames(m1) || !hub_id_col %in% colnames(m2)) {
    stop("Both gr1 and gr2 must have hub_id_col = ", hub_id_col)
  }
  
  # 确保 locus_id 存在
  if (is.null(locus_id_col) || !locus_id_col %in% colnames(m1)) {
    mcols(gr1)$locus_id <- paste0("L1_", seq_along(gr1))
    locus_id_col <- "locus_id"
  }
  if (!locus_id_col %in% colnames(m2)) {
    mcols(gr2)$locus_id <- paste0("L2_", seq_along(gr2))
  }
  
  # 构建 per-hub locus 集合
  df1 <- tibble(
    hub_id   = m1[[hub_id_col]],
    locus_id = mcols(gr1)[[locus_id_col]]
  ) %>%
    distinct()
  
  df2 <- tibble(
    hub_id   = m2[[hub_id_col]],
    locus_id = mcols(gr2)[[locus_id_col]]
  ) %>%
    distinct()
  
  # 仅对两个条件都存在的 hub 做 differential（交集）
  common_hubs <- intersect(df1$hub_id, df2$hub_id)
  if (length(common_hubs) == 0L) {
    warning("No common hub IDs between gr1 and gr2 for differential test.")
    return(tibble())
  }
  
  # locus universe 用 gr1+gr2 的所有 locus_id
  universe_loci <- unique(c(df1$locus_id, df2$locus_id))
  n_universe    <- length(universe_loci)
  
  # 预先按 hub 分组
  loci_by_hub1 <- df1 %>%
    filter(hub_id %in% common_hubs) %>%
    group_by(hub_id) %>%
    summarise(loci1 = list(locus_id), .groups = "drop")
  
  loci_by_hub2 <- df2 %>%
    filter(hub_id %in% common_hubs) %>%
    group_by(hub_id) %>%
    summarise(loci2 = list(locus_id), .groups = "drop")
  
  # 合并两个条件的列表
  hub_list <- loci_by_hub1 %>%
    inner_join(loci_by_hub2, by = "hub_id")
  
  # 对每个 hub 计算 observed D + permutation p-value
  res_list <- lapply(seq_len(nrow(hub_list)), function(i) {
    hub_id_i <- hub_list$hub_id[i]
    A        <- unique(hub_list$loci1[[i]])
    B        <- unique(hub_list$loci2[[i]])
    
    size1 <- length(A)
    size2 <- length(B)
    
    # 过滤太小的 hub
    if (size1 < min_size || size2 < min_size) {
      return(
        tibble(
          hub_id       = hub_id_i,
          size1        = size1,
          size2        = size2,
          n_shared     = length(intersect(A, B)),
          jaccard_obs  = NA_real_,
          D_obs        = NA_real_,
          pvalue       = NA_real_
        )
      )
    }
    
    J_obs <- jaccard_vec(A, B)
    D_obs <- 1 - J_obs
    n_shared_obs <- length(intersect(A, B))
    
    # permutation null：同样的 size，在 universe 上随机抽
    perm_D <- numeric(n_perm)
    for (k in seq_len(n_perm)) {
      A_rand <- sample(universe_loci, size1, replace = FALSE)
      B_rand <- sample(universe_loci, size2, replace = FALSE)
      J_rand <- jaccard_vec(A_rand, B_rand)
      perm_D[k] <- 1 - J_rand
    }
    
    pval <- (1 + sum(perm_D >= D_obs)) / (1 + n_perm)
    
    tibble(
      hub_id       = hub_id_i,
      size1        = size1,
      size2        = size2,
      n_shared     = n_shared_obs,
      jaccard_obs  = J_obs,
      D_obs        = D_obs,
      pvalue       = pval
    )
  })
  
  res_df <- bind_rows(res_list)
  # FDR
  res_df <- res_df %>%
    mutate(
      FDR = p.adjust(pvalue, method = "BH")
    )
  
  res_df
}
