## Lightweight BED export, so that rtracklayer is not needed for the core
## workflow. rtracklayer pulls in XML, which frequently fails to compile on
## systems lacking libxml2 development headers; making it optional removes that
## barrier for users who only need the standard pipeline.

## Write a data.frame with chrom/start/end (plus an optional name column) as BED.
## Coordinates are converted from 1-based inclusive (R / GRanges convention) to
## 0-based half-open (BED convention), matching rtracklayer::export().
.write_bed <- function(df, con, name_col = "name") {
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  ## locate the coordinate columns, tolerating the usual spellings
  chrom_col <- intersect(c("seqnames", "chrom", "chr", "chr1"), names(df))[1]
  start_col <- intersect(c("start", "chromStart", "start1", "x1"), names(df))[1]
  end_col   <- intersect(c("end", "chromEnd", "end1", "x2"), names(df))[1]

  if (any(is.na(c(chrom_col, start_col, end_col)))) {
    stop("Cannot find chrom/start/end columns for BED export.", call. = FALSE)
  }

  bed <- data.frame(
    chrom = as.character(df[[chrom_col]]),
    start = as.integer(df[[start_col]]) - 1L,   # BED is 0-based, half-open
    end   = as.integer(df[[end_col]]),
    stringsAsFactors = FALSE
  )

  if (!is.null(name_col) && name_col %in% names(df)) {
    bed$name <- as.character(df[[name_col]])
  }

  ## strand, if present and informative
  if ("strand" %in% names(df)) {
    s <- as.character(df[["strand"]])
    if (any(s %in% c("+", "-"))) {
      if (!"name" %in% names(bed)) bed$name <- "."
      bed$score  <- 0L
      bed$strand <- ifelse(s %in% c("+", "-"), s, ".")
    }
  }

  utils::write.table(bed, file = con, sep = "\t",
                     quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(con)
}
