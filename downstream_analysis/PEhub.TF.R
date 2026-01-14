#!/bin/bash
set -euo pipefail

############################################
# Activate environment
############################################
conda activate /mnt/citadel2/research/syidan/miniconda3/envs/bedtools

############################################
# Paths
############################################
BASE_DIR=~/syidan/Data/Processed/HiCHIP_brain_dif_part_GSE147672_softlink_merged
TF_DIR=${BASE_DIR}/scPrinter
OUT_BASE=${BASE_DIR}/Output
HUB_BASE=~/syidan/Projects/SnakeHichipResult/ProcessedData/multiple_enhancer
GENE_LIST=${BASE_DIR}/Output/expressed_gene_list

############################################
# Brain regions
############################################
REGIONS=(
  Caudate
  Hippocampus
  MiddleFrontalGyrus
  ParietalLobe
  SubstantiaNigra
  SuperiorTemporalGyri
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
OUTDIR=${OUT_BASE}/${REGION}
mkdir -p ${OUTDIR}
cd ${OUTDIR}

TFCSV=${TF_DIR}/${REGION}_TFBS_scores.csv
HUBBED=${HUB_BASE}/multiple_result.${REGION}.replicable_hubs_bin_union.bed

########################################
# 1. TF BED files
########################################
awk -F',' 'NR>1 {print $1,$2,$3,$4}' OFS="\t" \
${TFCSV} > ${REGION}_TFBS_scores.bed

awk -F',' 'NR>1 {
        mean=($5+$6)/2;
        printf "%s\t%d\t%d\t%s\t%.8g\n",$1,$2,$3,$4,mean
    }' ${TFCSV} > ${REGION}_TFBS_with_mean_scores.bed

awk '$5>0.1' \
${REGION}_TFBS_with_mean_scores.bed \
| grep -w -F -f "${GENE_LIST}/${REGION}_gene_list.txt" \
> ${REGION}_TFBS_with_mean_scores0.1.bed

########################################
# 2. Extend hubs (±1kb)
########################################
awk 'BEGIN{OFS="\t"}{$2-=1000;$3+=1000;print}' \
${HUBBED} > ${REGION}_hubs_extended1000bp.bed

########################################
# 3. Intersect hubs and TFs
########################################
bedtools intersect \
-a ${HUBBED} \
-b ${REGION}_TFBS_with_mean_scores0.1.bed \
-wa -wb > ${REGION}_hubs_tfs_overlap.bed

########################################
# 4. Hub–TF matrix
########################################
awk '{
    split($4, a, "_"); 
    hub = a[1]"_"a[2]"_"a[3]"_"a[4]; 
    print hub "\t" $10
}' ${REGION}_hubs_tfs_overlap.bed > hub_tf.txt

sort hub_tf.txt | uniq -c | awk '{print $2"\t"$3"\t"$1}' \
> hub_tf_count.txt

awk '
    {
        hub=$1; tf=$2; count=$3
        hubs[hub]; tfs[tf]; data[hub,tf]=count
    }
    END{
        n=asorti(tfs,TF)
        m=asorti(hubs,HUB,"@ind_num_asc")
        printf "hub"
        for(i=1;i<=n;i++) printf "\t%s",TF[i]
        print ""
        for(i=1;i<=m;i++){
            printf "%s",HUB[i]
            for(j=1;j<=n;j++)
                printf "\t%d",((HUB[i],TF[j]) in data ? data[HUB[i],TF[j]] : 0)
            print ""
        }
    }' hub_tf_count.txt > hub_tf_matrix_sorted.txt

########################################
# 5. TF hub percentage + mean score
########################################
awk '
    BEGIN{OFS="\t"}
    {
        split($4,a,"_"); hub=a[1]"_"a[2]"_"a[3]"_"a[4]
        tf=$10; score=$11
        hubs[hub]++
        seen[hub":"tf]=1
        sum[tf]+=score; count[tf]++
    }
    END{
        for(h in hubs) total++
        for(k in seen){
            split(k,a,":")
            hub_tf[a[2]]++
        }
        for(tf in hub_tf){
            printf "%s\t%.4f\t%.6f\n",
                   tf,
                   hub_tf[tf]/total*100,
                   sum[tf]/count[tf]
        }
    }' ${REGION}_hubs_tfs_overlap.bed \
| sort -k2,2nr > tf_summary.tsv

########################################
# 6. Synergistic hubs – fragment number
########################################
awk -F'\t' '
    {
        split($4,a,"_"); hub = a[1]"_"a[2]"_"a[3]"_"a[4]
        region=$1":"$2"-"$3":"$4
        hubs[hub][region]++
        lines[NR]=$0
        hub_line[NR]=hub
    }
    END{
        for(h in hubs){
            c=0
            for(r in hubs[h]) c++
            hub_count[h]=c
        }
        for(i=1;i<=NR;i++)
            print lines[i]"\thub_"hub_count[hub_line[i]]"fragments"
    }' ${HUBBED} > Synergistic_hubs_add_fragments_num.bed

########################################
# 7. Intersect synergistic hubs with TFs
########################################
bedtools intersect \
-a Synergistic_hubs_add_fragments_num.bed \
-b ${REGION}_TFBS_with_mean_scores0.1.bed \
-wa -wb | \
awk '{print $1,$2,$3,$4,$7,$11}' OFS="\t" | \
sort -u > Synergistic_hubs_${REGION}_TFBS_unique.bed

########################################
# 8. Boxplot input
########################################
awk -F'\t' '
    {
        split($4,a,"_"); hub = a[1]"_"a[2]"_"a[3]"_"a[4]
        if(!(hub in first)){
            first[hub]=$1"\t"$2"\t"$3
            frag[hub]=$5
        }
        tf[hub][$6]++
    }
    END{
        for(h in tf){
            c=0
            for(t in tf[h]) c++
            print first[h]"\t"h"\t"frag[h]"\t"c
        }
    }' Synergistic_hubs_${REGION}_TFBS_unique.bed \
| sort -t $'\t' -k1,1 -k2,2n \
> Synergistic_hubs_${REGION}_TFBS_boxplot.bed

echo "Finished ${REGION}"
done



