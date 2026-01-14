#######GWAS#########
sed 's/^#//; s/ ([^)]*)//g' ~/syidan/Data/Processed/GWAS_file/gwas_disease_data.bed |sort -u > cleaned_gwas_disease.bed

######overlap with hub file########
############################################
# Activate environment
############################################
conda activate /mnt/citadel2/research/syidan/miniconda3/envs/bedtools && cd ~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer

############################################
# Paths
############################################
BASE_DIR=~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged
GWAS_file=~/syidan/Data/Processed/GWAS_file/cleaned_gwas_disease.clean.bed
OUT_BASE=${BASE_DIR}/Output/regionSpecific_GWAS
HUB_BASE=~/syidan/Projects/SnakeHichipResult/Output/multiple_enhancer

############################################
# Brain regions
############################################
REGIONS=(
  class1_MiddleFrontalGyrus
  class1_SubstantiaNigra
  class2_MidSub
  class3_common3
)

############################################
# Loop over regions
############################################
for REGION in "${REGIONS[@]}"; do
echo "=============================="
echo "Processing ${REGION}"
echo "=============================="

########################################
# Setup
########################################
cd "${OUT_BASE}" || exit 1

HUBBED="${HUB_BASE}/brain.selected_region_overlap.${REGION}.bed"

########################################
# Intersect hubs and GWAS
########################################
bedtools intersect \
-a "${HUBBED}" \
-b "${GWAS_file}" \
-wa -wb \
| awk '{print $NF "\t" $0}' \
| sort -k1,1n \
| cut -f2- \
| awk 'BEGIN{OFS="\t"} {for (i=1; i<=NF; i++) if (i<4 || i>6) printf "%s%s", $i, (i==NF?ORS:OFS)}' \
| awk 'BEGIN{OFS="\t"} {
    # keep first 6 columns
    out = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6
    
    # merge columns from $7 to end with "_"
    merged=""
    for(i=7;i<=NF;i++){
        merged = (merged=="" ? $i : merged"_"$i)
    }
    
    # print first 6 columns + merged column
    print out, merged
}' | sort -k7,7 |uniq > "${REGION}_hubs_dif_Region_GWAS_overlap.bed"

echo "Finished ${REGION}"
done







#########FDR calculate##########
cd ~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged/Output/ && conda activate generalR_figures 


suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(stringr)
})

# ===============================
# Define all regions to process
# ===============================
REGIONS <- c(
  "class1_MiddleFrontalGyrus",
  "class1_SubstantiaNigra",
  "class2_MidSub",
  "class3_common3"
)

# Fixed GWAS file (same for all regions)
GWAS_FILE <- "~/syidan/Data/Processed/GWAS_file/cleaned_gwas_disease.clean.bed"

# Pattern for hub BED files
HUB_PATTERN <- "~/syidan/Projects/SnakeHichipResult/Output/multiple_enhancer/brain.selected_region_overlap.%s.bed"

# Single shared output directory
OUT_DIR <- "~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged/Output/regionSpecific_GWAS"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load and preprocess GWAS once (shared across regions)
gwas <- import(GWAS_FILE, extraCols = c(trait = "character"))

# Clean trait names once
mcols(gwas)$trait <- mcols(gwas)$trait %>%
  str_replace_all("\\s*\\(.*?\\)", "") %>%
  str_trim()

# Harmonize chromosome style once
seqlevelsStyle(gwas) <- "UCSC"

# ===============================
# Loop over each region
# ===============================
all_results <- list()  # Optional: collect all for a combined file at the end

for (REGION in REGIONS) {
  
  cat("\n=== Processing region:", REGION, "===\n")
  
  hub_file <- sprintf(HUB_PATTERN, REGION)
  
  if (!file.exists(hub_file)) {
    warning(paste("Hub file not found:", hub_file, "- Skipping", REGION))
    next
  }
  
  # Load hubs
  hubs <- import(hub_file)
  seqlevelsStyle(hubs) <- "UCSC"
  
  # Common chromosomes
  common_chr <- intersect(seqlevels(hubs), seqlevels(gwas))
  if (length(common_chr) == 0) {
    warning(paste("No common chromosomes for", REGION, "- Skipping"))
    next
  }
  
  hubs <- keepSeqlevels(hubs, common_chr, pruning.mode = "coarse")
  gwas_region <- keepSeqlevels(gwas, common_chr, pruning.mode = "coarse")
  
  # Overlap counts
  gwas_hits <- countOverlaps(gwas_region, hubs) > 0
  n_total_snps <- length(gwas_region)
  n_snps_in_hubs <- sum(gwas_hits)
  n_snps_outside_hubs <- n_total_snps - n_snps_in_hubs
  
  # Split by trait
  gwas_by_trait <- split(gwas_region, mcols(gwas_region)$trait)
  gwas_by_trait <- gwas_by_trait[sapply(gwas_by_trait, length) >= 5]
  
  if (length(gwas_by_trait) == 0) {
    message(paste("No traits with >=2 SNPs for", REGION))
    next
  }
  
  # Fisher test
  results <- lapply(names(gwas_by_trait), function(trait_name) {
    idx <- mcols(gwas_region)$trait == trait_name
    a <- sum(gwas_hits & idx)
    if (a == 0) return(NULL)
    
    b <- sum(idx) - a
    c <- n_snps_in_hubs - a
    d <- n_snps_outside_hubs - b
    
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    res <- fisher.test(mat, alternative = "greater")
    
    data.frame(
      trait = trait_name,
      region = REGION,  # Add region column
      snps_in_trait = a + b,
      snps_in_hubs = a,
      odds_ratio = unname(res$estimate),
      pvalue = res$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  region_df <- bind_rows(results)
  
  if (nrow(region_df) == 0) {
    message(paste("No enriched traits for", REGION))
    next
  }
  
  region_df <- region_df %>%
    mutate(FDR = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue)
  
  # Save region-specific files (with region name in filename)
  prefix <- paste0("gwas_fisher_", REGION)
  
  write.csv(
    region_df,
    file.path(OUT_DIR, paste0(prefix, "_all.csv")),
    row.names = FALSE
  )
  write.csv(
    region_df %>% filter(FDR < 0.05),
    file.path(OUT_DIR, paste0(prefix, "_FDR0.05.csv")),
    row.names = FALSE
  )
  
  # collect for combined file
  all_results[[REGION]] <- region_df
  
  cat("✅ Completed:", REGION,
      "- Tested", nrow(region_df), "traits,",
      sum(region_df$FDR < 0.05), "at FDR<0.05\n")
}




# After the for loop ends
all_sig_df <- bind_rows(all_results)

if (nrow(all_sig_df) > 0) {
  write.csv(
    all_sig_df,
    file.path(OUT_DIR, "gwas_fisher_all_regions_combined.csv"),
    row.names = FALSE
  )
  write.csv(
    all_sig_df %>% filter(FDR < 0.05),
    file.path(OUT_DIR, "gwas_fisher_all_regions_FDR0.05_combined.csv"),
    row.names = FALSE
  )
  cat("\nCombined significant results saved:", 
      nrow(all_sig_df), "enrichments across all regions.\n")
} else {
  cat("\nNo significant enrichments (FDR < 0.05) in any region.\n")
}










# 
# ###########make heatmap and radar plot based on FDR 0.05#################
# cd ~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged/Output/ && conda activate generalR_figures 
# 
# ########brain disease############
# ###############
# # Load required libraries
# library(dplyr)
# library(stringr)
# library(ggplot2)
# 
# all_region_df <- read.csv("~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged/Output/regionSpecific_GWAS/gwas_fisher_all_regions_FDR0.05_combined.csv",header = TRUE, stringsAsFactors = FALSE)
# 
# Neuro_Brain = c("Major depressive disorder or stress-related disorder", "Memory decline in normal cognition", "Insomnia", "Phoneme awareness", "General factor of neuroticism", "Neuroticism", "Short sleep duration", "Depressed affect", "Non-word reading", "Intelligence", "Fornix white matter microstructure", "Schizophrenia", "Anterior amygdaloid area volume", "Memory decline", "Subjective well-being", "Baseline memory in impaired cognition x sex interaction", "Decaffeinated coffee consumption and/or neuroticism")
# Lipids_Metabolites = c("Phosphatidylglycerol_[M+OAc]1- levels", "High density lipoprotein cholesterol levels", "Apolipoprotein A levels", "HDL cholesterol", "Omega-6 fatty acid levels", "Lipoprotein levels", "Lipoprotein A levels", "Polyunsaturated fatty acid levels", "Linoleic acid levels", "Valylleucine levels", "N4-acetylcytidine levels", "Pseudouridine levels", "Beta-endorphin levels", "Phosphatidylethanolamine_[M+H]1+ levels")
# Proteins_Level = c("IFNGR2 protein levels", "Tyrosine-protein phosphatase non-receptor type substrate 1 levels", "SIRPA protein levels", "Contactin-2 levels", "Beta-Ala-His dipeptidase levels", "CNTN2 protein levels", "CPB2 protein levels", "CNDP1 protein levels", "IL1R2 protein levels", "LEFTY2 protein levels", "ERAP1 protein levels", "Carboxypeptidase B2 levels", "CLUL1 protein levels", "CD5 levels", "Epidermal growth factor-like protein 6 levels", "CSTB protein levels", "NCAM2 protein levels", "ITGA2 protein levels", "Semaphorin-5A levels", "PON2 protein levels", "HSBP1 protein levels", "MOCS2 protein levels", "LRIG1 protein levels", "CST7 protein levels", "FGFBP2 protein levels", "BST1 protein levels", "ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase 2 levels", "GPR37 protein levels", "BTD protein levels", "Proprotein convertase subtilisin/kexin type 7 levels", "TNFRSF11A protein levels", "Biotinidase levels", "PCDH9 protein levels", "GRN protein levels", "NCAM1 protein levels", "HS1BP3 protein levels", "RNF41/WWP2 protein level ratio", "GLRX protein levels", "DPT protein levels", "CDKN2D/MANF protein level ratio", "Adhesion G protein-coupled receptor B3 levels", "Neogenin levels", "MELTF protein levels", "Dermatopontin levels", "Vitronectin levels", "A disintegrin and metalloproteinase with thrombospondin motifs 5 levels", "Semaphorin-3A levels", "Interferon alpha/beta receptor 1 levels", "IL10RB protein levels", "TMPRSS5 protein levels","ADAMTS8 protein levels", "GSTM4 protein levels", "Granulins levels", "LHPP protein levels", "KLB protein levels", "MAM domain-containing glycosylphosphatidylinositol anchor protein 1 levels", "CXCL5 levels", "GMPR protein levels", "CELSR2 protein levels", "vascular endothelial growth factor D levels", "Dihydropteridine reductase levels","NT3 levels")
# 
# Body_Physiology = c("Red blood cell erythrocyte count", "Corneal endothelial cell shape", "Monocyte percentage", "Visceral adipose tissue volumes to abdominal adipose tissue volumes ratio", "Visceral adipose tissue volumes", "Height", "Vertical cup-disc ratio", "Diastolic blood pressure", "Heel bone mineral density", "Waist circumference adjusted for BMI", "monocyte")
# Others = c("Normal pressure hydrocephalus", "Raw vegetable consumption", "Low myopia", "Alanine aminotransferase level after methotrexate initiation in rheumatoid arthritis", "Asthma-chronic obstructive  pulmonary disease overlap syndrome", "Polyneuropathy in diabetes", "Alcohol-related disorders", "Abdominal adipose tissue volumes to gluteofemoral adipose tissue volumes ratio", "Abdominal aortic aneurysm", "Heart valve disorders", "Metabolic syndrome", "Ventral hernia", "Coronary artery / coronary heart disease", "Heart rate response to beta blockers", "ER positive breast cancer x age at menarche interaction", "Substance use disorder", "Hand Osteoarthritis", "Breast cancer", "Cataracts", "Diverticulitis", "Pain intensity", "Pain intensity in opioid-treated advanced cancer", "Atrial fibrillation","BRCA1/2-negative high-risk breast cancer", "Thyrotoxic hypokalemic periodic paralysis and Graves disease", "Blood cell traits latent factor 5", "Blood cell traits latent factor 22", "Relative abundance of the human milk microbiota Pseudomonas hunanensis", "Relative abundance of the human milk microbiota Pseudomonas viridiflava", "Brugada syndrome","Occlusion and stenosis of precerebral arteries")
# 
# trait_cat_df <- data.frame(
#   trait = c(Lipids_Metabolites,Proteins_Level,Neuro_Brain,Body_Physiology,Others),
#   category = c(
#     rep("Lipids_Metabolites", length(Lipids_Metabolites)),
#     rep("Proteins_Level",length(Proteins_Level)),  # adjust counts as needed
#     rep("Neuro_Brain", length(Neuro_Brain)),
#     rep("Body_Physiology", length(Body_Physiology)),
#     rep("Others", length(Others))
#   ),
#   stringsAsFactors = FALSE
# )
# 
# all_region_df <- merge(all_region_df, trait_cat_df, by = "trait", all.x = TRUE)
# 
# all_region_df$category <- factor(all_region_df$category,
#                                  levels = c("Body_Physiology",
#                                             "Others",
#                                             "Lipids_Metabolites",
#                                             "Proteins_Level","Neuro_Brain"))
# 
# 
# # --- Summarize: number of unique traits per region and category ---
# donut_data <- all_region_df %>%
#   group_by(region, category) %>%
#   summarise(n_traits = n_distinct(trait), .groups = "drop") %>%
#   filter(n_traits > 0) %>%
#   group_by(region) %>%
#   mutate(percentage = round(100 * n_traits / sum(n_traits), 1),
#          total_traits = sum(n_traits)) %>%
#   ungroup()
# 
# # --- Colors ---
# cat_colors <- c(
#   "Lipids_Metabolites" = "#F4CE14",
#   "Proteins_Level"     = "#56B4E9",
#   "Neuro_Brain"        = "#FC5185",
#   "Body_Physiology"    = "#3FC1C9",
#   "Others"             = "#999999"
# )
# 
# # --- Create folder ---
# dir.create("region_donut_plots_pdf", showWarnings = FALSE)
# 
# # --- Generate one donut plot per region ---
# regions <- unique(donut_data$region)
# 
# for (reg in regions) {
#   data_reg <- donut_data %>% filter(region == reg)
#   
#   # Inner label: percentage + count
#   data_reg <- data_reg %>%
#     mutate(inner_label = paste0(percentage, "%\n(", n_traits, ")"))
#   
#   # Compute position for outer category labels (midpoint of each arc)
#   data_reg <- data_reg %>%
#     arrange(desc(category)) %>%
#     mutate(y_pos = cumsum(n_traits) - n_traits / 2,
#            # Angle for outer label (in radians)
#            angle = 2 * pi * y_pos / sum(n_traits) - pi / 2,
#            # Horizontal alignment based on position
#            hjust = ifelse(cos(angle) > 0, 0, 1),
#            # X position for outer label
#            x_pos = 2.8)  # slightly outside the donut
#   
#   p <- ggplot(data_reg, aes(x = 2, y = n_traits, fill = category)) +
#     geom_col(width = 1, color = "white", size = 1.2) +
#     coord_polar(theta = "y") +
#     xlim(0.5, 2.5) +  # reduced a bit since no outer labels needed
#     
#     # Inner labels: % and count
#     geom_text(aes(label = inner_label),
#               position = position_stack(vjust = 0.5),
#               color = "black",
#               size = 4.2,
#               fontface = "bold") +
#     
#     scale_fill_manual(values = cat_colors) +
#     
#     # Legend with category (trait) names
#     guides(fill = guide_legend(title = "Trait Category",   # you can change the title or remove it
#                                override.aes = list(color = "white", size = 1.2))) +
#     
#     theme_void() +
#     theme(
#       legend.position = "right",           # or "bottom" if you prefer
#       legend.title = element_text(face = "bold", size = 12),
#       legend.text = element_text(size = 11),
#       plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#       plot.subtitle = element_text(size = 12, hjust = 0.5),
#       plot.margin = margin(30, 30, 30, 30)  # less right margin needed
#     )
#   
#   # Optional titles (uncomment if desired)
#   # + labs(title = paste("Number of Unique Traits:", reg),
#   #        subtitle = paste("Total traits =", data_reg$total_traits[1]))
#   
#   # Save as PDF
#   safe_name <- str_replace_all(reg, "[^A-Za-z0-9]", "_")
#   filename <- paste0("regionSpecific_GWAS/donut_", safe_name, ".pdf")
#   
#   ggsave(filename, p, width = 6, height = 5, device = "pdf")  # widened a bit for legend space
#   cat("Saved:", filename, "\n")
# }
# 
# 
# 
# # =============================================
# # Heatmap: SNPs in hubs for Protein_Level traits across regions
# # =============================================
# library(dplyr)
# library(tidyr)
# library(stringr)
# library(pheatmap)  # Make sure pheatmap is installed: install.packages("pheatmap")
# 
# # --- Filter only protein traits ---
# protein_df <- all_region_df %>%
#   filter(category == "Proteins_Level") %>%
#   select(trait, region, snps_in_hubs) %>%
#   filter(!is.na(snps_in_hubs) & snps_in_hubs > 0)  # optional: remove zeros
# 
# # Check data
# cat("Found", length(unique(protein_df$trait)), "protein traits across",
#     length(unique(protein_df$region)), "regions.\n")
# 
# if (nrow(protein_df) == 0) {
#   stop("No data for Proteins_Level category.")
# }
# 
# # --- Prepare wide matrix ---
# heat_wide <- protein_df %>%
#   pivot_wider(
#     names_from = region,
#     values_from = snps_in_hubs,
#     values_fill = 0
#   )
# 
# # --- Convert to matrix using base R (no tibble needed) ---
# mat <- as.matrix(heat_wide[ , -1])              # all columns except trait
# rownames(mat) <- heat_wide$trait                # set trait names as row names
# 
# # --- Order rows and columns by total SNPs ---
# row_order <- order(rowSums(mat), decreasing = TRUE)
# mat <- mat[row_order, ]
# 
# #col_order <- order(colSums(mat), decreasing = TRUE)
# col_order <- c("class1_SubstantiaNigra", "class1_MiddleFrontalGyrus", "class2_MidSub", "class3_common3" )
# 
# mat <- mat[, col_order]
# 
# # --- Generate and save the heatmap ---
# pdf("regionSpecific_GWAS/Proteins_snps_in_hubs_heatmap.pdf", width = 18, height = 32)
# 
# pheatmap(mat,
#          color = colorRampPalette(c("white","#FF7DB0", "#FF0087"))(100),
#          cluster_rows = T,
#          cluster_cols = FALSE,
#          show_rownames = TRUE,
#          show_colnames = TRUE,
#          fontsize_row = 24,
#          fontsize_col = 28,
#          border_color = "grey60",
#          legend = TRUE,
#          fontsize = 20,
#          cellwidth = 30,
#          cellheight = 30)
# 
# dev.off()
# 
# 
# 
# 
# 
# # =============================================
# # Heatmap: SNPs in hubs for Lipids_Metabolites traits across regions
# 
# Lipids_Metabolites_df <- all_region_df %>%
#   filter(category == "Lipids_Metabolites") %>%
#   select(trait, region, snps_in_hubs) %>%
#   filter(!is.na(snps_in_hubs) & snps_in_hubs > 0)  # optional: remove zeros
# 
# if (nrow(Lipids_Metabolites_df) == 0) {
#   stop("No data for Lipids_Metabolites category.")
# }
# 
# # --- Prepare wide matrix ---
# heat_wide <- Lipids_Metabolites_df %>%
#   pivot_wider(
#     names_from = region,
#     values_from = snps_in_hubs,
#     values_fill = 0
#   )
# 
# # --- Convert to matrix using base R (no tibble needed) ---
# mat <- as.matrix(heat_wide[ , -1])              # all columns except trait
# rownames(mat) <- heat_wide$trait%>%
#   gsub(" levels$", "", .)                
# 
# # --- Order rows and columns by total SNPs ---
# row_order <- order(rowSums(mat), decreasing = TRUE)
# mat <- mat[row_order, ]
# 
# col_order <- c("class1_SubstantiaNigra", "class1_MiddleFrontalGyrus", "class2_MidSub", "class3_common3" )
# mat <- mat[, col_order]
# 
# 
# 
# # mat[mat > 50] <- 50
# # --- Generate and save the heatmap ---
# pdf("regionSpecific_GWAS/lipid_snps_in_hubs_heatmap.pdf", width = 6, height = 6)
# 
# pheatmap(mat,
#          color = colorRampPalette(c("white","#FF7DB0","#FF0087","#FF0087"))(100),
#          cluster_rows = T,
#          cluster_cols = FALSE,
#          show_rownames = TRUE,
#          show_colnames = TRUE,
#          fontsize_row = 12,
#          fontsize_col = 12,
#          border_color = "grey60",
#          legend = TRUE,
#          fontsize = 16,
#          cellwidth = 20,
#          cellheight = 20)
# 
# dev.off()
# 
# 
# 
# 
# # =============================================
# # Heatmap: SNPs in hubs for brain traits across regions
# 
# Neuro_Brain_df <- all_region_df %>%
#   filter(category == "Neuro_Brain") %>%
#   select(trait, region, snps_in_hubs) %>%
#   filter(!is.na(snps_in_hubs) & snps_in_hubs > 0)  # optional: remove zeros
# 
# if (nrow(Neuro_Brain_df) == 0) {
#   stop("No data for Neuro_Brain category.")
# }
# 
# # --- Prepare wide matrix ---
# heat_wide <- Neuro_Brain_df %>%
#   pivot_wider(
#     names_from = region,
#     values_from = snps_in_hubs,
#     values_fill = 0
#   )
# 
# # --- Convert to matrix using base R (no tibble needed) ---
# mat <- as.matrix(heat_wide[ , -1])              # all columns except trait
# rownames(mat) <- heat_wide$trait%>%
#   gsub(" levels$", "", .)                
# 
# # --- Order rows and columns by total SNPs ---
# row_order <- order(rowSums(mat), decreasing = TRUE)
# mat <- mat[row_order, ]
# 
# col_order <- c("class1_MiddleFrontalGyrus", "class1_SubstantiaNigra", "class2_MidSub", "class3_common3" )
# 
# mat <- mat[, col_order]
# 
# 
# 
# # mat[mat > 50] <- 50
# # --- Generate and save the heatmap ---
# pdf("regionSpecific_GWAS/Neuro_Brain_snps_in_hubs_heatmap.pdf", width = 10, height = 12)
# 
# pheatmap(mat,
#          color = colorRampPalette(c("white","#8C00FF50","#8C00FF"))(100),
#          cluster_rows = T,
#          cluster_cols = FALSE,
#          show_rownames = TRUE,
#          show_colnames = TRUE,
#          fontsize_row = 12,
#          fontsize_col = 12,
#          border_color = "grey60",
#          legend = TRUE,
#          fontsize = 16,
#          cellwidth = 20,
#          cellheight = 20)
# 
# dev.off()
# 
# 
















###########make heatmap and radar plot based on Pvalue 0.01 odds_ratio #################
cd ~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged/Output/ && conda activate generalR_figures 

########brain disease############
###############
# Load required libraries
library(dplyr)
library(stringr)
library(ggplot2)

all_region_df <- read.csv("~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged/Output/regionSpecific_GWAS/gwas_fisher_all_regions_combined.csv",header = TRUE, stringsAsFactors = FALSE)
all_region_df = subset(all_region_df, all_region_df$pvalue <= 0.01&odds_ratio>=5)
all_region_df <- all_region_df[order(all_region_df$odds_ratio, decreasing = TRUE), ]


# unique(subset(all_region_df,all_region_df$region == "class1_SubstantiaNigra")$trait)
# unique(all_region_df$trait)


Neuro_Brain = c("Major depressive disorder or stress-related disorder", "Memory decline in normal cognition", "Insomnia", "General factor of neuroticism", "Neuroticism", "Short sleep duration", "Depressed affect", "Numerical cognitive ability", "Cerebral cortical growth", "Brain age gap from white matter microstructure", "Habitual snoring", "Phoneme awareness", "Non-word reading", "Fornix white matter microstructure", "Anterior amygdaloid area volume", "Memory decline", "Subjective well-being", "Baseline memory in impaired cognition x sex interaction", "Decaffeinated coffee consumption and/or neuroticism", "Age at dementia onset in PSEN1 E280A mutation carriers", "Aggressiveness in attention deficit hyperactivity disorder", "Feeling tense", "Obesity class II and Attention deficit hyperactivity disorder or Anorexia nervosa or Major depressive disorder or Obsesive-compulsive disorder or Schizophrenia", "Autism spectrum disorder, attention deficit-hyperactivity disorder, bipolar disorder, major depressive disorder, and schizophrenia", "Irritable bowel syndrome or major depressive disorder", "Affective disorder and suicide attempts", "Anxiety disorders", "Normal pressure hydrocephalus", "Polyneuropathy in diabetes", "Frontotemporal dementia", "Pediatric central nervous system tumors", "Schizophrenia", "Schizophrenia vs ADHD)","Coffee difference liking", "Cake icing liking", "Globe artichoke liking", "Alcohol-related disorders", "Opiates and related narcotics causing adverse effects in therapeutic use", "Spicy food liking", "Cerebrospinal fluid AB1-42 levels", "wg rh intensity-contrast parsopercularis", "Vertical cup-disc ratio", "Temporal pole thickness", "Cerebellum cortex volume change rate x age interaction")
Lipids_Metabolites = c("Phosphatidylglycerol_[M+OAc]1- levels", "Valylleucine levels", "Leucylglycine levels", "N4-acetylcytidine levels", "Pseudouridine levels", "3-methylcytidine levels", "Nervonoylcarnitine levels", "Dihomo-linolenoylcarnitine levels", "Succinylcarnitine levels", "Acylcarnitine levels", "Beta-endorphin levels", "Valylglycine levels", "Gamma-glutamyl-2-aminobutyrate levels", "Lysoalkylphosphatidylcholine levels", "Sulfatide levels", "Cholesteryl ester levels in large VLDL", "Cholesteryl ester levels in very large VLDL", "Total lipid levels in small VLDL", "Total lipids in VLDL", "Phospholipid levels in VLDL", "Large HDL particle concentration", "Cholesteryl esters in medium HDL", "Cholesterol in medium HDL", "Phosphatidylcholine_32:2_[M+H]1+/Phosphatidylethanolamine_35:2_[M+H]1+/Phosphatidate_37:3_[M+NH4]1+ levels", "Arabonate/xylonate levels", "1,2-dilinoleoyl-GPC levels", "Pregnenediol disulfate levels", "Symmetrical dimethylarginine levels", "Phosphatidylethanolamine_38:6_[M-H]1- levels", "Phosphatidylethanolamine_36:4_[M-H]1- levels", "Lipoprotein levels", "Lipoprotein A levels")

Proteins_Level = c("Dihydropteridine reductase levels", "Contactin-2 levels", "HS1BP3 protein levels", "Beta-Ala-His dipeptidase levels", "Epidermal growth factor-like protein 6 levels", "Dermatopontin levels", "RNF41/WWP2 protein level ratio", "Tyrosine-protein phosphatase non-receptor type substrate 1 levels", "Biotinidase levels", "Platelet-derived growth factor receptor-like protein levels", "A disintegrin and metalloproteinase with thrombospondin motifs 5 levels", "Adhesion G protein-coupled receptor B3 levels", "Semaphorin-5A levels", "Neogenin levels", "CLEC1B/PPP1R2 protein level ratio", "NADH-cytochrome b5 reductase 2 levels", "Osteopontin levels", "Protocadherin-9 levels", "GMPR protein levels", "IFNGR2 protein levels", "ITGA2 protein levels", "Interferon alpha/beta receptor 1 levels", "ERAP1 protein levels", "CDKN2D/MANF protein level ratio", "CD69/CLEC1B protein level ratio", "CHMP1A/EREG protein level ratio", "Cerebral dopamine neurotrophic factor levels", "Dickkopf-related protein 3 levels", "EREG/MPI protein level ratio", "FLI1/SH2B3 protein level ratio", "LMNB2 protein levels", "MANF/PPIB protein level ratio", "NFATC1/SPRY2 protein level ratio", "SRPK2/VASH1 protein level ratio", "Sarcalumenin levels", "Carboxypeptidase B2 levels", "vascular endothelial growth factor D levels", "QDPR protein levels", "DPT protein levels", "Proprotein convertase subtilisin/kexin type 7 levels", "CYB5R2 protein levels", "CD5 levels", "SIRPA protein levels", "ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase 2 levels", "DOK2/MESD protein level ratio", "ANGPT1/APP protein level ratio", "AXIN1/IRAK4 protein level ratio", "Beta-klotho levels", "CALCOCO1/INPPL1 protein level ratio", "CD69/MPIG6B protein level ratio", "DFFA protein levels", "GP1BA/HBEGF protein level ratio", "MED18/SH2B3 protein level ratio", "PFKFB2 protein levels", "Nidogen-1 levels", "Carboxypeptidase Q levels", "CPB2 protein levels", "Granulins levels", "Neural cell adhesion molecule 1, 120 kDa isoform levels", "Glutathione S-transferase Mu 3 levels", "Vascular cell adhesion protein 1 levels", "Testican-3 levels", "Vitronectin levels", "CNDP1 protein levels", "MMP9 protein levels", "CCN3 protein levels", "GRN protein levels", "CXCL5 levels", "Scavenger receptor cysteine-rich type 1 protein M130 levels", "CLGN protein levels", "NFASC protein levels", "BTD protein levels", "PCDH9 protein levels", "CELSR2 protein levels", "LHPP protein levels", "IL10RB protein levels", "CLUL1 protein levels", "MAM domain-containing glycosylphosphatidylinositol anchor protein 1 levels", "AXIN1/HEXIM1 protein level ratio", "EREG protein levels", "NTRK2 protein levels", "IL1R2 protein levels", "PON2 protein levels", "GPR37 protein levels", "LEFTY2 protein levels", "CSTB protein levels", "NCAM2 protein levels", "SERPINB8 protein levels", "NCAM1 protein levels", "GSTM4 protein levels", "BST1 protein levels", "HSBP1 protein levels", "VSIG10 protein levels", "HPR protein levels", "PCDH17 protein levels", "ICOSLG protein levels", "CD83 protein levels", "KLB protein levels", "BPIFB1 protein levels", "CD38 protein levels", "DKK3 protein levels", "ROBO1 protein levels", "PON1 protein levels", "MEGF9 protein levels", "FGFBP2 protein levels", "CST7 protein levels", "ITGBL1 protein levels", "PROC protein levels", "CD36 protein levels", "CLEC4G protein levels", "TMPRSS5 protein levels", "SEMA3G protein levels", "FUT8 protein levels", "ENPP2 protein levels", "ALDH3A1 protein levels", "CCDC80 protein levels", "CLEC1B/SNAP29 protein level ratio","CTSO protein levels","Semaphorin-3A levels" ,"CNTN2 protein levels","GLRX protein levels", "LRIG1 protein levels", "MELTF protein levels", "MOCS2 protein levels","NT3 levels","Occlusion and stenosis of precerebral arteries", "TNFRSF11A protein levels","ADAMTS8 protein levels")

Body_Physiology = c("Subscapular skin fold thickness", "Corneal endothelial cell shape", "Bone marrow fat fraction of thoracic vertebra 10", "Platelet side scatter", "Immature platelet count", "Platelet side fluorescence", "Erythrocyte sedimentation rate", "Physical activity", "Red blood cell erythrocyte count", "Monocyte percentage", "Visceral adipose tissue volumes to abdominal adipose tissue volumes ratio", "Visceral adipose tissue volumes", "Height", "Heel bone mineral density", "Waist circumference adjusted for BMI", "Monocyte", "Abdominal adipose tissue volumes to gluteofemoral adipose tissue volumes ratio", "Posterior urethral valves", "Abdominal adipose tissue volumes", "Tonometry", "Ideal cardiovascular health score", "Relative abundance of the human milk microbiota Sphingobacterium multivorum", "Relative abundance of the human milk microbiota Pseudomonas viridiflava", "Relative abundance of the human milk microbiota Pseudomonas hunanensis", "C-reactive protein and white blood cell count", "Gut microbial network clusters x Older Siblings interaction","circulating leptin levels", "Energy expenditure")
Others = c("ER positive breast cancer x age at menarche interaction", "Asthma-chronic obstructive  pulmonary disease overlap syndrome", "Anemia of chronic disease", "Ingrowing nail", "Chronic ulcer of leg or foot", "Molar-incisor hypomineralization", "Neoplasm of unspecified nature of digestive system", "OTU99_30 abundance", "Postoperative survival time in hepatocellular carcinoma", "Cardiovascular risk factors", "Central serous chorioretinopathy", "Chronic ulcer of skin", "Edema", "Chronic laryngitis", "Midgestational circulating levels of PCBs", "Beef consumption", "Dermatophytosis of nail", "Heart failure with reduced EF [Systolic or combined heart failure]", "Congestive heart failure; nonhypertensive", "Hypertensive heart and/or renal disease", "Childhood asthma exacerbations in long-acting beta2-agonist treatment", "Lesions of stomach and duodenum", "Response to tofacitinib treatment in rheumatoid arthritis", "SDF-1 levels in metastatic colorectal cancer", "Ventral hernia", "Arm fat ratio or colorectal cancer", "Hair colour: Black", "Coronary artery / coronary heart disease", "Hand Osteoarthritis", "Breast cancer in BRCA2 mutation carriers", "Thyrotoxic hypokalemic periodic paralysis and Graves disease", "Cutaneous leishmaniasis", "Melanomas of skin", "Heart valve disorders", "Heart rate response to beta blockers", "Alanine aminotransferase level after methotrexate initiation in rheumatoid arthritis", "Raw vegetable consumption", "Diverticulitis", "Brugada syndrome", "Mucinous adenocarcinoma in colorectal cancer", "Severe insulin-deficient type 2 diabetes", "Keloid", "X-23593 levels", "Response to methotrexate in juvenile idiopathic arthritis", "Pelvic organ prolapse x age interaction", "B-cell acute lymphoblastic leukaemia", "Thiazide-induced adverse metabolic effects in hypertensive patients", "Airflow obstruction", "Heart attack", "Ankle injury", "Thyroid peroxidase antibody levels in pregnancy", "Other chronic ischemic heart disease, unspecified", "Clear cell renal cell carcinoma", "Polycystic ovary syndrome", "Obstructive sleep apnea", "Ovarian cancer", "Calcific aortic valve stenosis", "Endometriosis or asthma", "Angina", "Cataracts", "Clonal hematopoiesis", "IgG glycosylation patterns", "X-23593 levels", "X-21383 levels", "X-13728 levels","Cough","Insulin resistance/response","Left-hemisphere salience/ventral attention network to right-hemisphere salience/ventral attention network white-matter structural connectivity","Coronary artery disease or tissue plasminogen activator levels","Low myopia", "Pain intensity in opioid-treated advanced cancer","BRCA1/2-negative high-risk breast cancer")

trait_cat_df <- data.frame(
  trait = c(Lipids_Metabolites,Proteins_Level,Neuro_Brain,Body_Physiology,Others),
  category = c(
    rep("Lipids_Metabolites", length(Lipids_Metabolites)),
    rep("Proteins_Level",length(Proteins_Level)),  # adjust counts as needed
    rep("Neuro_Brain", length(Neuro_Brain)),
    rep("Body_Physiology", length(Body_Physiology)),
    rep("Others", length(Others))
  ),
  stringsAsFactors = FALSE
)

all_region_df <- merge(all_region_df, trait_cat_df, by = "trait", all.x = TRUE)

all_region_df$category <- factor(all_region_df$category,
                                 levels = c("Body_Physiology",
                                            "Others",
                                            "Lipids_Metabolites",
                                            "Proteins_Level","Neuro_Brain"))


# --- Summarize: number of unique traits per region and category ---
donut_data <- all_region_df %>%
  group_by(region, category) %>%
  summarise(n_traits = n_distinct(trait), .groups = "drop") %>%
  filter(n_traits > 0) %>%
  group_by(region) %>%
  mutate(percentage = round(100 * n_traits / sum(n_traits), 1),
         total_traits = sum(n_traits)) %>%
  ungroup()

# --- Colors ---
cat_colors <- c(
  "Lipids_Metabolites" = "#F4CE14",
  "Proteins_Level"     = "#56B4E9",
  "Neuro_Brain"        = "#FC5185",
  "Body_Physiology"    = "#3FC1C9",
  "Others"             = "#999999"
)

# --- Create folder ---
dir.create("region_donut_plots_pdf", showWarnings = FALSE)

# --- Generate one donut plot per region ---
regions <- unique(donut_data$region)

for (reg in regions) {
  data_reg <- donut_data %>% filter(region == reg)
  
  # Inner label: percentage + count
  data_reg <- data_reg %>%
    mutate(inner_label = paste0(percentage, "%\n(", n_traits, ")"))
  
  # Compute position for outer category labels (midpoint of each arc)
  data_reg <- data_reg %>%
    arrange(desc(category)) %>%
    mutate(y_pos = cumsum(n_traits) - n_traits / 2,
           # Angle for outer label (in radians)
           angle = 2 * pi * y_pos / sum(n_traits) - pi / 2,
           # Horizontal alignment based on position
           hjust = ifelse(cos(angle) > 0, 0, 1),
           # X position for outer label
           x_pos = 2.8)  # slightly outside the donut
  
  p <- ggplot(data_reg, aes(x = 2, y = n_traits, fill = category)) +
    geom_col(width = 1, color = "white", size = 1.2) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +  # reduced a bit since no outer labels needed
    
    # Inner labels: % and count
    geom_text(aes(label = inner_label),
              position = position_stack(vjust = 0.5),
              color = "black",
              size = 4.2,
              fontface = "bold") +
    
    scale_fill_manual(values = cat_colors) +
    
    # Legend with category (trait) names
    guides(fill = guide_legend(title = "Trait Category",   # you can change the title or remove it
                               override.aes = list(color = "white", size = 1.2))) +
    
    theme_void() +
    theme(
      legend.position = "right",           # or "bottom" if you prefer
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.margin = margin(30, 30, 30, 30)  # less right margin needed
    )
  
  # Optional titles (uncomment if desired)
  # + labs(title = paste("Number of Unique Traits:", reg),
  #        subtitle = paste("Total traits =", data_reg$total_traits[1]))
  
  # Save as PDF
  safe_name <- str_replace_all(reg, "[^A-Za-z0-9]", "_")
  filename <- paste0("regionSpecific_GWAS/donut_", safe_name, ".pdf")
  
  ggsave(filename, p, width = 6, height = 5, device = "pdf")  # widened a bit for legend space
  cat("Saved:", filename, "\n")
}



# =============================================
# Heatmap: SNPs in hubs for Protein_Level traits across regions
# =============================================
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)  # Make sure pheatmap is installed: install.packages("pheatmap")

# --- Filter only protein traits ---
protein_df <- all_region_df %>%
  filter(category == "Proteins_Level") %>%
  select(trait, region, odds_ratio) %>%
  filter(!is.na(odds_ratio) & odds_ratio > 0)  # optional: remove zeros

# Check data
cat("Found", length(unique(protein_df$trait)), "protein traits across",
    length(unique(protein_df$region)), "regions.\n")

if (nrow(protein_df) == 0) {
  stop("No data for Proteins_Level category.")
}

# --- Prepare wide matrix ---
heat_wide <- protein_df %>%
  pivot_wider(
    names_from = region,
    values_from = odds_ratio,
    values_fill = 0
  )

# --- Convert to matrix using base R (no tibble needed) ---
mat <- as.matrix(heat_wide[ , -1])              # all columns except trait
rownames(mat) <- heat_wide$trait                # set trait names as row names

# --- Order rows and columns by total SNPs ---
row_order <- order(rowSums(mat), decreasing = TRUE)
mat <- mat[row_order, ]

#col_order <- order(colSums(mat), decreasing = TRUE)
col_order <- c("class1_SubstantiaNigra", "class1_MiddleFrontalGyrus", "class2_MidSub", "class3_common3" )

mat <- mat[, col_order]

# --- Generate and save the heatmap ---
pdf("regionSpecific_GWAS/Proteins_snps_in_hubs_heatmap.pdf", width = 15, height = 58)

pheatmap(mat,
         color = colorRampPalette(c("white","#FF7DB0", "#FF0087", "#FF0087"))(100),
         cluster_rows = T,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 26,
         fontsize_col = 28,
         border_color = "grey60",
         legend = TRUE,
         fontsize = 20,
         cellwidth = 30,
         cellheight = 30)

dev.off()





# =============================================
# Heatmap: SNPs in hubs for Lipids_Metabolites traits across regions

Lipids_Metabolites_df <- all_region_df %>%
  filter(category == "Lipids_Metabolites") %>%
  select(trait, region, odds_ratio) %>%
  filter(!is.na(odds_ratio) & odds_ratio > 0)  # optional: remove zeros

if (nrow(Lipids_Metabolites_df) == 0) {
  stop("No data for Lipids_Metabolites category.")
}

# --- Prepare wide matrix ---
heat_wide <- Lipids_Metabolites_df %>%
  pivot_wider(
    names_from = region,
    values_from = odds_ratio,
    values_fill = 0
  )

# --- Convert to matrix using base R (no tibble needed) ---
mat <- as.matrix(heat_wide[ , -1])              # all columns except trait
rownames(mat) <- heat_wide$trait%>%
  gsub(" levels$", "", .)                

# --- Order rows and columns by total SNPs ---
row_order <- order(rowSums(mat), decreasing = TRUE)
mat <- mat[row_order, ]

col_order <- c("class1_SubstantiaNigra", "class1_MiddleFrontalGyrus", "class2_MidSub", "class3_common3" )
mat <- mat[, col_order]



# mat[mat > 50] <- 50
# --- Generate and save the heatmap ---
pdf("regionSpecific_GWAS/lipid_snps_in_hubs_heatmap.pdf", width = 6, height = 10)

pheatmap(mat,
         color = colorRampPalette(c("white","#FF7DB0","#FF0087","#FF0087","#FF0087"))(100),
         cluster_rows = T,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 12,
         fontsize_col = 12,
         border_color = "grey60",
         legend = TRUE,
         fontsize = 16,
         cellwidth = 20,
         cellheight = 20)

dev.off()





# =============================================
# Heatmap: SNPs in hubs for brain traits across regions

Neuro_Brain_df <- all_region_df %>%
  filter(category == "Neuro_Brain") %>%
  select(trait, region, odds_ratio) %>%
  filter(!is.na(odds_ratio) & odds_ratio > 0)  # optional: remove zeros

if (nrow(Neuro_Brain_df) == 0) {
  stop("No data for Neuro_Brain category.")
}

# --- Prepare wide matrix ---
heat_wide <- Neuro_Brain_df %>%
  pivot_wider(
    names_from = region,
    values_from = odds_ratio,
    values_fill = 0
  )

# --- Convert to matrix using base R (no tibble needed) ---
mat <- as.matrix(heat_wide[ , -1])              # all columns except trait
rownames(mat) <- heat_wide$trait%>%
  gsub(" levels$", "", .)                

# --- Order rows and columns by total SNPs ---
row_order <- order(rowSums(mat), decreasing = TRUE)
mat <- mat[row_order, ]

col_order <- c("class1_MiddleFrontalGyrus", "class1_SubstantiaNigra", "class2_MidSub", "class3_common3" )

mat <- mat[, col_order]


# mat[mat > 50] <- 50
# --- Generate and save the heatmap ---
pdf("regionSpecific_GWAS/Neuro_Brain_snps_in_hubs_heatmap.pdf", width = 14, height = 14)

pheatmap(mat,
         color = colorRampPalette(c("white","#8C00FF50","#8C00FF"))(100),
         cluster_rows = T,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 12,
         fontsize_col = 12,
         border_color = "grey60",
         legend = TRUE,
         fontsize = 16,
         cellwidth = 20,
         cellheight = 20)

dev.off()




########### Radar plots per region using Neuro_Brain_df (original style, top 10 traits) #################

library(dplyr)
library(ggplot2)

# ===============================
# Define colors for regions
# ===============================
region_colors <- c(
  class1_MiddleFrontalGyrus = "#8C00FF",
  class1_SubstantiaNigra    = "#0BA6DF",
  class2_MidSub             = "#DD0303",
  class3_common3            = "#FA891A"
)

# ===============================
# Loop over each region
# ===============================
regions <- unique(Neuro_Brain_df$region)

for (REGION in regions) {
  
  # Subset data for this region
  region_df <- all_region_df %>%
  filter(category == "Neuro_Brain") %>%
  select(trait, region, snps_in_hubs) %>%
  filter(region == REGION) %>%
  arrange(desc(snps_in_hubs)) %>%        # sort by snps_in_hubs
  slice_head(n = 10)                  # top 10 traits
  
  # Only keep traits present
  disease_trait_order <- region_df$trait
  top_n <- length(disease_trait_order)
  
  # ===============================
  # Per-region scaling
  # ===============================
  max_val <- max(region_df$snps_in_hubs, na.rm = TRUE)+0.5
  ylim_max <- max_val * 1.3
  label_radius <- max_val * 1.15
  
  # ===============================
  # Prepare radar data (including center for original style)
  # ===============================
  radar_data <- region_df %>%
    mutate(
      trait = factor(trait, levels = disease_trait_order),
      id = row_number(),
      angle = (id - 1) * 2 * pi / top_n,
      value_norm = snps_in_hubs
    )
  
  # Line data to connect points back to center
  line_data <- bind_rows(radar_data, mutate(radar_data, value_norm = 0))
  
  # Radial lines (spokes)
  radial_lines <- radar_data %>%
    transmute(angle = angle, y_start = 0, y_end = max_val)
  
  col <- region_colors[REGION]
  
  # ===============================
  # Build radar plot (original style)
  # ===============================
  p <- ggplot() +
    
    # Spokes
    geom_segment(data = radial_lines,
                 aes(x = angle, xend = angle, y = y_start, yend = y_end),
                 color = "grey30", linetype = "dashed", linewidth = 0.5) +
    
    # Main polygon and points (connect to center)
    geom_line(data = line_data,
              aes(x = angle, y = value_norm, group = id),
              color = col, linewidth = 2) +
    geom_point(data = radar_data,
               aes(x = angle, y = value_norm),
               size = 4, color = col) +
    
    # Trait labels
    geom_text(data = radar_data,
              aes(x = angle, y = label_radius, label = trait),
              size = 4.5, fontface = "bold", hjust = 0, vjust = 0.5) +
    
    # Value labels
    geom_text(data = radar_data,
              aes(x = angle, y = value_norm + max_val * 0.05, label = round(snps_in_hubs,1)),
              size = 6, fontface = "bold", color = col) +
    
    # Polar coordinates
    coord_polar(theta = "x", start = -pi / top_n) +
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      plot.margin = margin(60, 60, 60, 60)
    ) +
    ggtitle(REGION) +
    ylim(0, ylim_max) +
    
    # Concentric circles
    annotate("path", x = seq(0, 2*pi, length.out = 240), y = rep(max_val*0.25, 240),
             linetype = "dashed", color = "grey30") +
    annotate("path", x = seq(0, 2*pi, length.out = 240), y = rep(max_val*0.50, 240),
             linetype = "dashed", color = "grey30") +
    annotate("path", x = seq(0, 2*pi, length.out = 240), y = rep(max_val*0.75, 240),
             linetype = "dashed", color = "grey30") +
    annotate("path", x = seq(0, 2*pi, length.out = 240), y = rep(max_val, 240),
             color = "grey10", linewidth = 1)
  
  # ===============================
  # Save plot
  # ===============================
  ggsave(
    filename = paste0("regionSpecific_GWAS/radar_plot_", REGION, "_top10_originalStyle.pdf"),
    plot = p,
    width = 6,
    height = 6,
    dpi = 300
  )
  
  cat("✅ Radar plot generated for region:", REGION, "(top 10 traits, original style)\n")
}












###########make radar plot #################
cd ~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged/Output/ && conda activate generalR_figures 



#odds_ratio
library(dplyr)
library(stringr)
library(ggplot2)

# ===============================
# Define regions and colors
# ===============================
REGIONS <- c(
  "Caudate",
  "Hippocampus",
  "MiddleFrontalGyrus",
  "ParietalLobe",
  "SubstantiaNigra",
  "SuperiorTemporalGyri"
)

# Assign a unique color to each region (you can customize these)
region_colors <- c(
  Caudate = "#E37434",              
  Hippocampus = "#8C00FF",          
  MiddleFrontalGyrus = "#1FAB89",   
  ParietalLobe = "#0BA6DF",         
  SubstantiaNigra = "#DD0303",      
  SuperiorTemporalGyri = "#00B7B5"  
)

# Path template
BASE_PATH <- "~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged/Output/GWAS_enrichment"

# Read combined significant results to define trait order and max value
all_region_df <- read.csv(file.path(BASE_PATH, "gwas_fisher_all_regions_combined.csv"),header = TRUE, stringsAsFactors = FALSE)

# Define trait order (unique, in desired display order)
disease_trait_order <- unique(c("General factor of neuroticism","Educational attainment","Intelligence","Depressed affect","Maternal history of Alzheimer's disease","Attention deficit hyperactivity disorder or autism spectrum disorder or intelligence","Life satisfaction","Subjective well-being","Memory decline in normal cognition","Schizophrenia","Neuroticism","Insomnia"))

top_n <- length(disease_trait_order)

# Compute global max for consistent scaling across all plots
overlap_all_traits <- all_region_df %>% filter(trait %in% disease_trait_order)
max_val <- max(overlap_all_traits$odds_ratio, 0)
label_radius <- max_val * 1.15
ylim_max <- max_val * 1.30



#subset(overlap_all_traits,overlap_all_traits$region=="SuperiorTemporalGyri")


# ===============================
# Loop over each region
# ===============================
for (REGION in REGIONS) {
  cat("Generating radar plot for:", REGION, "\n")
  
  # Read region-specific full results
  region_file <- file.path(BASE_PATH, paste0("gwas_fisher_", REGION, "_all.csv"))
  if (!file.exists(region_file)) {
    warning(paste("File not found:", region_file))
    next
  }
  region_df <- read.csv(region_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Subset to traits of interest and complete missing traits
  radar_data <- data.frame(trait = disease_trait_order, stringsAsFactors = FALSE) %>%
    left_join(region_df %>% select(trait, odds_ratio), by = "trait") %>%
    mutate(
      odds_ratio = ifelse(is.na(odds_ratio) | odds_ratio < 1, 1, odds_ratio),  # missing or no enrichment → OR = 1
      value_norm = odds_ratio,  # use raw odds ratio as radial value
      trait = factor(trait, levels = disease_trait_order)
    ) %>%
    arrange(trait) %>%
    mutate(
      id = row_number(),
      angle = (id - 1) * (2 * pi / top_n)
    )
  
  # Line data (connect to center)
  line_data <- bind_rows(radar_data, mutate(radar_data, value_norm = 1))  # close at OR=1, not 0
  
  # Radial lines (spokes) - extend to global max
  radial_lines <- radar_data %>%
    transmute(angle = angle, y_start = 1, y_end = max_val)  # start from OR=1 line
  
  # Get color for this region
  col <- region_colors[REGION]
  
  # Build plot
  p <- ggplot() +
    # Concentric circles
    annotate("path", x = seq(0, 2*pi, length.out = 240), y = rep(1, 240), 
             color = "grey30", linewidth = 1) +  # highlight OR=1
    annotate("path", x = seq(0, 2*pi, length.out = 240), y = rep(max_val*0.33, 240), 
             linetype = "dashed", color = "grey30") +
    annotate("path", x = seq(0, 2*pi, length.out = 240), y = rep(max_val*0.66, 240), 
             linetype = "dashed", color = "grey30") +
    annotate("path", x = seq(0, 2*pi, length.out = 240), y = rep(max_val, 240), 
             color = "grey10", linewidth = 1) +
    
    # Spokes
    geom_segment(data = radial_lines,
                 aes(x = angle, xend = angle, y = y_start, yend = y_end),
                 color = "grey30", linetype = "dashed", linewidth = 0.6) +
    
    # Main polygon and points
    geom_line(data = line_data,
              aes(x = angle, y = value_norm, group = id),
              color = col, linewidth = 2) +
    geom_point(data = radar_data,
               aes(x = angle, y = value_norm),
               size = 4, color = col) +
    
    # Trait labels (better alignment)
    geom_text(data = radar_data,
              aes(x = angle, y = label_radius, label = trait,
                  hjust = ifelse(angle > pi, 1, 0)),  # left/right align
              size = 6, fontface = "bold", vjust = 0.5) +
    
    # Value labels (only show if OR > 2 to reduce clutter)
    geom_text(data = radar_data %>% filter(odds_ratio > 2),
              aes(x = angle, y = value_norm + max_val * 0.05, 
                  label = round(odds_ratio, 1)),
              size = 5, fontface = "bold", color = col) +
    
    # Polar coordinates
    coord_polar(theta = "x", start = -pi / top_n) +
    theme_void() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 30)),
      plot.margin = margin(60, 60, 60, 60)
    ) +
    ggtitle(REGION) +
    ylim(1, ylim_max)  # start from OR=1
  
  # Save
  plot_file <- file.path(BASE_PATH, paste0("radar_plot_", REGION, "_odds_ratio.pdf"))
  ggsave(plot_file, p, width = 10, height = 10, device = "pdf", dpi = 300)
  
  cat("Saved:", plot_file, "\n")
}
cat("All radar plots generated successfully!\n")