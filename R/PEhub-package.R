#' PEhub: promoter-centric detection of multi-way enhancer hubs
#'
#' PEhub reconstructs multi-way enhancer hubs from pairwise chromatin-interaction
#' data by modelling enhancer-promoter interactions as promoter-anchored weighted
#' networks, quantifying enhancer co-membership, partitioning enhancers into
#' candidate hubs via Leiden community detection, and evaluating hubs against a
#' distance-matched null model.
#'
#' @keywords internal
"_PACKAGE"

## Imports for the core functions, which were kept verbatim from the original
## script and therefore use unqualified calls. Users never need library().
##
## Note: data.table, dplyr and igraph export several functions with the same
## names (first, last, between, union, groups, as_data_frame). The resulting
## "replacing previous import" messages at load time are harmless: the later
## import wins, which is the same resolution the original script had when it
## attached these packages in the same order.
#'
#' @import data.table
#' @import dplyr
#' @import igraph
#' @importFrom Matrix Matrix sparseMatrix
#' @importFrom tidyr unnest fill
#' @importFrom purrr map map2 map_dbl map_lgl map_int map_dfr pmap reduce keep compact
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyselect all_of any_of
#' @importFrom rlang .data
#' @importFrom stats quantile p.adjust ppois sd median setNames var rpois ave
#' @importFrom utils head tail write.table read.delim
#' @importFrom methods as is
#' @importFrom GenomicRanges GRanges promoters findOverlaps mcols seqnames start end makeGRangesFromDataFrame countOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
NULL
