## =============================================================================
## Pre-flight power/depth diagnostics for pooled low-depth / single-cell input.
##
## PEhub's background model, weighting schemes and null-model permutation
## testing were developed and validated on deeply sequenced bulk HiChIP and
## Pore-C data (e.g. 9 paired human-heart H3K27ac HiChIP samples, GM12878
## Pore-C), where "pooled" effectively means millions of input cells. When
## the input instead comes from pooling a modest number of single cells or
## nuclei (e.g. single-cell/single-nucleus Hi-C-type assays such as
## snm3C-seq), the resulting interaction counts can be far sparser than what
## the background model and permutation testing were calibrated for,
## especially at longer genomic distances. This file adds an explicit,
## opt-out pre-flight check that surfaces this as a diagnostic rather than
## leaving it undocumented; it does not change any downstream statistics.
## =============================================================================

#' Diagnose statistical power for pooled low-depth / single-cell input
#'
#' PEhub's background model and significance testing were developed and
#' validated on deeply sequenced bulk HiChIP / Pore-C data (effectively
#' millions of input cells). When the input instead comes from pooling a
#' modest number of single cells or nuclei (e.g. single-cell/single-nucleus
#' Hi-C-type assays), interaction counts can be far sparser than what the
#' background model and null-model permutation testing were calibrated for,
#' especially at longer genomic distances. This function performs a
#' pre-flight check on the full interaction set (the same file normally
#' passed as \code{loop_file_all}) and flags depth ranges where hub calling
#' may not be reliable.
#'
#' The default thresholds are not arbitrary: \code{min_total_reads = 2e6} is
#' the depth below which a real single-donor snm3C-seq calibration run
#' (two Microglia donors, 658K vs. 5.35M valid contact pairs, otherwise
#' identical pipeline) showed hub sizes collapsing to the \code{k_min} floor
#' rather than forming a real distribution, and \code{min_bin_reads = 1e5}
#' reflects the point at which a background model fit to pooled single-donor
#' data was found to become numerically degenerate (expected count and
#' variance collapsing toward machine epsilon) specifically beyond ~50kb.
#' Both are deliberately conservative defaults for a pooled-single-cell
#' regime, not estimates of a universal minimum for all assays.
#'
#' @param loop_file_all Path to the full (background) interaction set, in
#'   the same BEDPE-like format used by \code{\link{preprocess_hichip}}
#'   (first six columns \code{chr1,start1,end1,chr2,start2,end2}; a
#'   \code{counts} column is used as the per-row read weight if present,
#'   otherwise each row counts as one).
#' @param n_cells Optional. The number of single cells/nuclei pooled to
#'   produce \code{loop_file_all}, if known. Included in the report only;
#'   total read/contact depth, not raw cell count, is what determines
#'   statistical power, since per-cell contact recovery varies by assay.
#' @param breaks Distance bin edges, in bp, used only for the diagnostic
#'   report. Defaults match \code{\link{compute_weights}}'s default
#'   distance bins.
#' @param min_total_reads Minimum total contact count (summed \code{counts},
#'   or row count if no \code{counts} column) below which the whole sample
#'   is flagged low-power. Default \code{2e6}.
#' @param min_bin_reads Minimum total contact count within a single distance
#'   bin below which that bin specifically is flagged unreliable, even if
#'   the sample's total depth clears \code{min_total_reads}. Default
#'   \code{1e5}.
#' @param long_range_from_bin Index (1-based) of the first distance bin,
#'   among \code{breaks}, considered "long-range" for the purposes of the
#'   \code{"long_range_unreliable"} flag. Default \code{4}, i.e. bins
#'   starting at the 4th break (50kb with the default \code{breaks}),
#'   matching the distance at which the background-model degeneracy above
#'   was observed to begin.
#'
#' @return Invisibly, a list with \code{total_reads}, \code{n_cells},
#'   \code{by_bin} (a data frame of per-distance-bin read counts and a
#'   \code{reliable} flag), \code{overall_flag} (one of \code{"ok"},
#'   \code{"low_power"}, \code{"long_range_unreliable"}) and \code{message}
#'   (the human-readable summary also printed as a side effect).
#'
#' @seealso \code{\link{pehub_prepare_interactions}}, which calls this check
#'   automatically (\code{check_power = TRUE} by default) and attaches its
#'   result to the returned object as \code{$power_check}.
#' @export
#' @examples
#' \dontrun{
#' pehub_check_power("all_interactions.bedpe.gz", n_cells = 66)
#' }
pehub_check_power <- function(loop_file_all,
                              n_cells = NA_integer_,
                              breaks = c(0, 1e4, 2.5e4, 5e4, 1e5, 2.5e5, 5e5, 1e6, 2e6),
                              min_total_reads = 2e6,
                              min_bin_reads = 1e5,
                              long_range_from_bin = 4) {
  d <- data.table::fread(loop_file_all, header = TRUE)
  if (ncol(d) < 6) {
    stop("'", loop_file_all, "' has fewer than 6 columns; expected a BEDPE-like ",
         "file with chr1,start1,end1,chr2,start2,end2 as the first six.", call. = FALSE)
  }
  data.table::setnames(d, 1:6, c("chr1", "start1", "end1", "chr2", "start2", "end2"))

  has_counts <- "counts" %in% names(d)
  if (!has_counts) {
    message("pehub_check_power(): no 'counts' column found; treating each row as one read/contact.")
  }
  d[, .pehub_weight := if (has_counts) counts else 1L]
  d[, .pehub_dist := abs((start1 + end1) / 2 - (start2 + end2) / 2)]
  d[, .pehub_bin := cut(.pehub_dist, breaks = breaks, include.lowest = TRUE, right = FALSE)]

  by_bin <- d[, .(n_pairs = .N, total_reads = sum(.pehub_weight, na.rm = TRUE)), by = .pehub_bin]
  data.table::setnames(by_bin, ".pehub_bin", "bin")
  by_bin <- by_bin[order(bin)]
  by_bin[, reliable := total_reads >= min_bin_reads]

  total_reads <- sum(d$.pehub_weight, na.rm = TRUE)

  long_range_bins <- by_bin$bin[seq_len(nrow(by_bin)) >= long_range_from_bin]
  long_range_unreliable <- any(!by_bin$reliable[by_bin$bin %in% long_range_bins])

  overall_flag <- if (total_reads < min_total_reads) {
    "low_power"
  } else if (long_range_unreliable) {
    "long_range_unreliable"
  } else {
    "ok"
  }

  msg <- paste0(
    "PEhub power check: ", format(total_reads, big.mark = ",", scientific = FALSE),
    " total reads/contacts",
    if (!is.na(n_cells)) paste0(" pooled from ", n_cells, " cells/nuclei") else "",
    ".\n",
    switch(overall_flag,
      low_power = paste0(
        "  ** LOW POWER: total depth is below the ", format(min_total_reads, big.mark = ",", scientific = FALSE),
        "-read\n     threshold at which hub calling was observed to degrade to the k_min floor in\n",
        "     real single-donor calibration data (see the vignette's 'Single-cell / pooled\n",
        "     low-depth input' section). Treat hub sizes and significance calls from this\n",
        "     sample with caution, and avoid comparing hub size/count across samples that\n",
        "     differ substantially in this total. **\n"),
      long_range_unreliable = paste0(
        "  ** LONG-RANGE CAUTION: one or more distance bins from ",
        format(breaks[long_range_from_bin], big.mark = ",", scientific = FALSE),
        "bp onward have fewer than\n     ", format(min_bin_reads, big.mark = ",", scientific = FALSE),
        " total reads. The background/null model is least reliable at long range\n",
        "     under low total depth; enhancer-hub-relevant (long-range) calls in the\n",
        "     flagged bins below should be treated as provisional. **\n"),
      ok = "  Depth looks adequate across all tested distance bins.\n"
    )
  )
  message(msg)

  invisible(list(
    total_reads = total_reads,
    n_cells = n_cells,
    by_bin = as.data.frame(by_bin),
    overall_flag = overall_flag,
    message = msg
  ))
}
