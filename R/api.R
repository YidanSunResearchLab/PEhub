## =============================================================================
## User-facing API.
##
## These are descriptive names for the four pipeline stages. Each one both
## writes its RData checkpoint (as before) AND returns the key objects, so the
## user does not have to load() anything to inspect results.
##
## The original build.* functions are retained as deprecated aliases.
## =============================================================================

#' Stage 1: prepare enhancer-promoter interactions
#'
#' Reads significant and background chromatin interactions, annotates each anchor
#' as promoter (P) or enhancer (E) using the supplied TSS annotation, and builds
#' the promoter-anchored enhancer-promoter (EP) interaction tables used by all
#' downstream stages.
#'
#' @param loop_file Path to the significant interactions in BEDPE format.
#' @param loop_file_all Path to the full (background) interaction set, used to
#'   build the distance-matched null pool.
#' @param outdir Output directory. Created if it does not exist.
#' @param tss_input A \code{GRanges} of transcription start sites, with a
#'   \code{gene_id} metadata column. See \code{\link{build_tss_from_txdb}}.
#' @param filename Sample name; used as a prefix for all output files.
#' @param promoter_window Half-width, in bp, of the promoter window around each
#'   TSS. \code{0} treats the TSS itself as the promoter (default).
#'
#' @return Invisibly, a list with \code{prep_hichip} (significant EP tables) and
#'   \code{all_ep_data} (background EP tables). Also writes
#'   \code{multiple_result.<filename>.hub.all.preprocess.RData} to \code{outdir}.
#'
#' @seealso \code{\link{pehub_detect_hubs}} for the next stage.
#' @export
pehub_prepare_interactions <- function(loop_file, loop_file_all, outdir,
                                       tss_input, filename,
                                       promoter_window = 0) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  build.run.preprocess(loop_file       = loop_file,
                       loop_file_all   = loop_file_all,
                       outdir          = outdir,
                       tss_input       = tss_input,
                       filename        = filename,
                       promoter_window = promoter_window)

  ## reload the checkpoint we just wrote, so we can hand the objects back
  e <- new.env(parent = emptyenv())
  load(file.path(outdir, paste("multiple_result", filename, "hub",
                               "all.preprocess.RData", sep = ".")), envir = e)
  res <- list(prep_hichip = e$prep_hichip, all_ep_data = e$all_ep_data)
  message("Stage 1 complete: ", nrow(res$prep_hichip$loops),
          " significant EP interactions across ",
          length(unique(res$prep_hichip$loops$promoter_id)), " promoters.")
  invisible(res)
}


#' Stage 2: detect candidate enhancer hubs
#'
#' Computes edge weights for each enhancer-promoter interaction, builds the
#' per-promoter enhancer co-membership network, partitions it into communities
#' with the Leiden algorithm, retains the highest-ranking community (by total
#' node strength) as the candidate hub for each promoter, and prepares the
#' distance-matched null pool used in stage 3.
#'
#' @param outdir Output directory used in stage 1.
#' @param filename Sample name used in stage 1.
#' @param weight_method Edge-weighting scheme. One of the eleven schemes; see
#'   \code{\link{pehub_weight_schemes}} for the full list. Default
#'   \code{"bin_log_ratio_sig"} (Distance-adjusted log-enrichment).
#' @param method Co-membership normalisation method. Default \code{"log_minmax"}.
#' @param k_min Minimum number of enhancers for a promoter to be considered.
#'   \code{3} is the definitional lower bound for multi-way assembly.
#' @param resolution Leiden resolution parameter. Default \code{1.0}; values
#'   above this progressively fragment hubs (see the resolution sensitivity
#'   analysis in the manuscript).
#' @param quantile_cutoff Sparsification threshold: co-membership edges below
#'   this quantile of non-zero edge weights are removed. Default \code{0.2}.
#' @param B_pvalue Number of permutations for the distance-matched null model.
#'   Default \code{1000}.
#' @param promoter_window Passed through for bookkeeping; default \code{0}.
#'
#' @return Invisibly, a list with \code{hubs} (all candidate hubs),
#'   \code{hub1} (the retained top-ranked hub per promoter) and
#'   \code{null_prepared} (the null-model inputs). Also writes
#'   \code{multiple_result.<filename>.hub.<method>.<weight_method>.preprocess.RData}.
#'
#' @seealso \code{\link{pehub_evaluate_hubs}} for the next stage.
#' @export
pehub_detect_hubs <- function(outdir, filename,
                              weight_method   = "bin_log_ratio_sig",
                              method          = "log_minmax",
                              k_min           = 3,
                              resolution      = 1.0,
                              quantile_cutoff = 0.2,
                              B_pvalue        = 1000,
                              promoter_window = 0) {
  build.run.all(outdir          = outdir,
                method          = method,
                weight_method   = weight_method,
                filename        = filename,
                promoter_window = promoter_window,
                k_min           = k_min,
                resolution      = resolution,
                quantile_cutoff = quantile_cutoff,
                B_pvalue        = B_pvalue)

  e <- new.env(parent = emptyenv())
  load(file.path(outdir, paste("multiple_result", filename, "hub",
                               method, weight_method, "preprocess.RData",
                               sep = ".")), envir = e)
  res <- list(
    hubs          = e$observed_hubs_per_promoter,
    hub1          = e$observed_hubs_per_promoter_sub_hub,
    null_prepared = e$observed_hubs_per_promoter_pvalue_prepare
  )
  message("Stage 2 complete: ", length(unique(res$hub1$promoter_id)), " candidate hubs (one per promoter).")
  invisible(res)
}


#' Stage 3: evaluate candidate hubs for significance and reproducibility
#'
#' Computes, for each candidate hub, (i) an empirical p-value against a
#' distance-matched null model and (ii) a reproducibility rate from bootstrap
#' re-discovery, then merges these statistics onto the hub table.
#'
#' @param outdir,filename,method,weight_method As in \code{\link{pehub_detect_hubs}}.
#' @param k_min,resolution,quantile_cutoff As in \code{\link{pehub_detect_hubs}};
#'   these must match stage 2, since hubs are re-discovered during bootstrapping.
#' @param B_pvalue Permutations for the distance-matched null. Default \code{1000}.
#' @param B_stability Bootstrap iterations for reproducibility. Default \code{10}.
#' @param null_mode Null construction. \code{"hist_matched"} (default) preserves
#'   each hub's distance-bin composition.
#' @param pvalue_cutoff,stability_cutoff Retained for interface symmetry; the
#'   thresholds are applied in \code{\link{pehub_export_hubs}}.
#'
#' @return Invisibly, a list with \code{hub1} (hubs annotated with
#'   \code{hub_p_value_global}, \code{hub_p_adj_global},
#'   \code{reproducibility_rate}, \code{OE_ratio_global}), plus the separate
#'   \code{pvalues} and \code{stability} tables. Also writes
#'   \code{multiple_result.<filename>.hub.<method>.<weight_method>.RData}.
#'
#' @seealso \code{\link{pehub_export_hubs}} for the next stage.
#' @export
pehub_evaluate_hubs <- function(outdir, filename,
                                weight_method    = "bin_log_ratio_sig",
                                method           = "log_minmax",
                                k_min            = 3,
                                resolution       = 1.0,
                                quantile_cutoff  = 0.2,
                                B_pvalue         = 1000,
                                B_stability      = 10,
                                null_mode        = "hist_matched",
                                pvalue_cutoff    = 0.05,
                                stability_cutoff = 0.5) {
  build.cutoff(outdir           = outdir,
               method           = method,
               weight_method    = weight_method,
               filename         = filename,
               k_min            = k_min,
               resolution       = resolution,
               quantile_cutoff  = quantile_cutoff,
               B_pvalue         = B_pvalue,
               B_stability      = B_stability,
               null_mode        = null_mode,
               pvalue_cutoff    = pvalue_cutoff,
               stability_cutoff = stability_cutoff)

  e <- new.env(parent = emptyenv())
  load(file.path(outdir, paste("multiple_result", filename, "hub",
                               method, weight_method, "RData",
                               sep = ".")), envir = e)
  res <- list(
    hub1      = e$observed_hubs_per_promoter_sub_hub,
    pvalues   = e$observed_hubs_per_promoter_pvalue,
    stability = e$observed_hubs_per_promoter_stability
  )
  n_sig <- sum(res$hub1$hub_p_value_global <= pvalue_cutoff &
               res$hub1$reproducibility_rate >= stability_cutoff, na.rm = TRUE)
  message("Stage 3 complete: ", length(unique(res$hub1$promoter_id)), " candidate hubs evaluated")
  invisible(res)
}


#' Stage 4: classify and export high-confidence hubs
#'
#' Applies the significance and reproducibility thresholds to classify candidate
#' hubs as high-confidence, and writes BED, BEDPE and tab-delimited exports for
#' both the high-confidence hubs and the remaining pairwise interactions.
#'
#' @param outdir,filename,method,weight_method As in \code{\link{pehub_detect_hubs}}.
#' @param pvalue_cutoff Significance threshold. Default \code{0.05}.
#' @param stability_cutoff Minimum reproducibility rate. Default \code{0.5}.
#' @param use_adjusted_p If \code{FALSE} (default) filter on the unadjusted
#'   p-value (\code{hub_p_value_global}), as in the manuscript. If \code{TRUE},
#'   filter on the Benjamini-Hochberg adjusted q-value instead.
#'
#' @return Invisibly, a list with \code{high_confidence} (the retained hubs),
#'   \code{pairwise} (interactions not assigned to a high-confidence hub) and
#'   \code{files} (paths of everything written).
#'
#' @export
pehub_export_hubs <- function(outdir, filename,
                              weight_method    = "bin_log_ratio_sig",
                              method           = "log_minmax",
                              pvalue_cutoff    = 0.05,
                              stability_cutoff = 0.5,
                              use_adjusted_p   = FALSE) {
  res <- build.postprocess(outdir           = outdir,
                    method           = method,
                    weight_method    = weight_method,
                    filename         = filename,
                    pvalue_cutoff    = pvalue_cutoff,
                    stability_cutoff = stability_cutoff)

  stem <- paste("multiple_result", filename, "hub", method, weight_method, sep = ".")
  pstem <- paste("multiple_result", filename, "pairwise", method, weight_method, sep = ".")
  files <- list(
    hub_bed    = file.path(outdir, paste0(stem, ".bed")),
    hub_txt    = file.path(outdir, paste0(stem, ".txt")),
    hub_bedpe  = file.path(outdir, paste0(stem, ".bedpe")),
    pair_bed   = file.path(outdir, paste0(pstem, ".bed")),
    pair_txt   = file.path(outdir, paste0(pstem, ".txt")),
    pair_bedpe = file.path(outdir, paste0(pstem, ".bedpe"))
  )

message("Stage 4 complete: ",
          length(unique(res$high_confidence_pval$promoter_id)), " high-confidence hubs (unadjusted p < ", pvalue_cutoff, "); ",
          length(unique(res$high_confidence_fdr$promoter_id)),  " high-confidence hubs (FDR < ", pvalue_cutoff, "). ",
          "Exported to ", outdir)
  invisible(list(high_confidence_pval = res$high_confidence_pval,
                 high_confidence_fdr  = res$high_confidence_fdr,
                 files                = files))
  }


#' Run the complete PEhub pipeline
#'
#' Convenience wrapper that runs all four stages in order.
#'
#' @inheritParams pehub_prepare_interactions
#' @inheritParams pehub_evaluate_hubs
#' @param use_adjusted_p See \code{\link{pehub_export_hubs}}.
#' @param workers Number of parallel workers for the null model. Requires the
#'   \pkg{future} and \pkg{furrr} packages; if absent, runs serially with
#'   identical results.
#'
#' @return Invisibly, the list returned by \code{\link{pehub_export_hubs}}.
#' @export
pehub_run <- function(loop_file, loop_file_all, outdir, tss_input, filename,
                      weight_method    = "bin_log_ratio_sig",
                      method           = "log_minmax",
                      promoter_window  = 0,
                      k_min            = 3,
                      resolution       = 1.0,
                      quantile_cutoff  = 0.2,
                      B_pvalue         = 1000,
                      B_stability      = 10,
                      null_mode        = "hist_matched",
                      pvalue_cutoff    = 0.05,
                      stability_cutoff = 0.5,
                      workers          = 1) {

  if (workers > 1) {
    if (requireNamespace("future", quietly = TRUE) &&
        requireNamespace("furrr", quietly = TRUE)) {
      future::plan(future::multisession, workers = workers)
      on.exit(future::plan(future::sequential), add = TRUE)
    } else {
      warning("future/furrr not installed; running serially.", call. = FALSE)
    }
  }

  pehub_prepare_interactions(loop_file, loop_file_all, outdir, tss_input,
                             filename, promoter_window)
  pehub_detect_hubs(outdir, filename, weight_method, method, k_min, resolution,
                    quantile_cutoff, B_pvalue, promoter_window)
  pehub_evaluate_hubs(outdir, filename, weight_method, method, k_min, resolution,
                      quantile_cutoff, B_pvalue, B_stability, null_mode,
                      pvalue_cutoff, stability_cutoff)
  out <- pehub_export_hubs(outdir, filename, weight_method, method,
                           pvalue_cutoff, stability_cutoff)
  pehub_shutdown()
  invisible(out)
}


#' Available edge-weighting schemes
#'
#' Returns a data frame documenting the eleven interaction weighting schemes,
#' their human-readable names, and what each one uses.
#'
#' @return A data frame with columns \code{code}, \code{label},
#'   \code{category} and \code{description}.
#' @export
#' @examples
#' pehub_weight_schemes()
pehub_weight_schemes <- function() {
  data.frame(
    code = c("count_only", "sig_only", "distance_only",
             "count_sig", "count_sig_plus_dist_linear",
             "bin_percentile", "bin_percentile_plus_sig",
             "bin_diff_global", "bin_diff_binmax", "bin_log_ratio_sig",
             "zscore_residual"),
    label = c("Count", "Significance", "Distance (control)",
              "Count-Significance", "Count-significance-distance",
              "Distance-adjusted percentile",
              "Distance-adjusted percentile-significance",
              "Distance-adjusted excess signal (global)",
              "Distance-adjusted excess signal (bin-normalized)",
              "Distance-adjusted log-enrichment",
              "Z-score residual"),
    category = c(rep("Single-feature", 3),
                 rep("Count-significance", 2),
                 rep("Distance-adjusted", 5),
                 "Model residual"),
    description = c(
      "Raw contact frequency",
      "Statistical confidence, -log10(q), capped",
      "Inverse log-distance; negative control",
      "Counts weighted by statistical significance",
      "Linear integration of count, significance and distance",
      "Percentile rank of count within distance bins",
      "Distance percentile weighted by significance",
      "Contact excess above distance-specific expectation",
      "Distance-adjusted excess, normalized within distance bins",
      "Log enrichment over distance-specific expectation, weighted by significance (default)",
      "Standardized residual from the HiC-DC+ background model"
    ),
    stringsAsFactors = FALSE
  )
}
