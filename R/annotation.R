## TSS annotation helpers.
##
## preprocess_hichip() takes a GRanges of TSS positions. Two convenience helpers
## are provided for building it. Neither is required: users may construct the
## GRanges however they like.

#' Read transcription start sites from a GTF file
#'
#' Parses a GTF/GFF file and returns a \code{GRanges} of transcription start
#' sites, in the form expected by \code{\link{preprocess_hichip}}. Uses
#' \pkg{rtracklayer} when available, and otherwise falls back to a lightweight
#' \pkg{data.table} parser with no additional dependencies.
#'
#' @param gtf_file Path to a GTF or GFF file. May be gzipped.
#' @param feature Feature type to extract TSS from. Default \code{"transcript"}.
#' @param attributes Character vector of GTF attributes to retain as metadata
#'   columns. Default keeps \code{gene_id}, \code{gene_name} and
#'   \code{transcript_id}.
#'
#' @return A \code{GRanges} of 1-bp TSS positions.
#'
#' @examples
#' \dontrun{
#' tss <- read_tss_from_gtf("genes.gtf.gz")
#' }
#' @export
read_tss_from_gtf <- function(gtf_file,
                              feature    = "transcript",
                              attributes = c("gene_id", "gene_name", "transcript_id")) {

  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    gtf <- rtracklayer::import(gtf_file)
    gtf <- gtf[as.character(gtf$type) == feature]
    tss <- GenomicRanges::resize(gtf, width = 1, fix = "start")
    return(tss)
  }

  ## ---- fallback: parse the GTF directly, no rtracklayer needed --------------
  message("rtracklayer not installed; using the built-in GTF parser.")
  dt <- data.table::fread(
    gtf_file, sep = "\t", header = FALSE, skip = "#",
    col.names = c("seqnames", "source", "type", "start", "end",
                  "score", "strand", "frame", "attribute"),
    showProgress = FALSE
  )
  dt <- dt[dt$type == feature]
  if (nrow(dt) == 0L) {
    stop("No records of type '", feature, "' found in ", gtf_file, call. = FALSE)
  }

  ## pull the requested attributes out of the 9th column
  extract_attr <- function(x, key) {
    m <- regmatches(x, regexpr(paste0(key, ' "[^"]*"'), x))
    out <- rep(NA_character_, length(x))
    hit <- lengths(m) > 0 | nzchar(m)
    out[hit] <- sub(paste0('^', key, ' "'), "", sub('"$', "", m[hit]))
    out
  }
  meta <- lapply(attributes, function(k) extract_attr(dt$attribute, k))
  names(meta) <- attributes

  ## TSS = start for + strand, end for - strand
  tss_pos <- ifelse(dt$strand == "-", dt$end, dt$start)

  gr <- GenomicRanges::GRanges(
    seqnames = dt$seqnames,
    ranges   = IRanges::IRanges(start = tss_pos, width = 1L),
    strand   = dt$strand
  )
  for (k in attributes) S4Vectors::mcols(gr)[[k]] <- meta[[k]]
  gr
}


#' Build a TSS GRanges from a TxDb annotation
#'
#' Convenience wrapper that extracts transcription start sites from a TxDb, in
#' the form expected by \code{\link{preprocess_hichip}}. Requires the relevant
#' annotation packages, which are in \code{Suggests}.
#'
#' @param txdb A TxDb object. Defaults to
#'   \code{TxDb.Hsapiens.UCSC.hg38.knownGene} if installed.
#'
#' @return A \code{GRanges} of TSS positions with a \code{gene_id} column.
#' @export
build_tss_from_txdb <- function(txdb = NULL) {
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop("Package 'GenomicFeatures' is required.\n",
         "Install with: BiocManager::install('GenomicFeatures')", call. = FALSE)
  }
  if (is.null(txdb)) {
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
      stop("Provide a 'txdb', or install 'TxDb.Hsapiens.UCSC.hg38.knownGene'.",
           call. = FALSE)
    }
    txdb <- get("TxDb.Hsapiens.UCSC.hg38.knownGene",
                envir = asNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene"))
  }
  tx  <- GenomicFeatures::transcripts(txdb, columns = c("gene_id", "tx_name"))
  GenomicRanges::resize(tx, width = 1, fix = "start")
}
