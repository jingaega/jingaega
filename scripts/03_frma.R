# fRMA normalisation (PLAN §1). Per-array: frozen probe effects from hgu133plus2frmavecs,
# so no parameter is fitted on this dataset. Usage: Rscript 03_frma.R <cel_dir> <out_tsv>
suppressMessages({library(affy); library(frma); library(hgu133plus2frmavecs)})
args <- commandArgs(trailingOnly = TRUE)
cel_dir <- args[1]; out <- args[2]
files <- sort(list.files(cel_dir, pattern = "CEL.gz$", full.names = TRUE))
cat("arrays:", length(files), "\n")
res <- list()
chunk <- split(files, ceiling(seq_along(files) / 20))   # memory; fRMA is per-array so chunking is exact
for (i in seq_along(chunk)) {
  ab <- ReadAffy(filenames = chunk[[i]])
  es <- frma(ab, summarize = "robust_weighted_average")
  res[[i]] <- exprs(es)
  cat("chunk", i, "of", length(chunk), "\n")
}
m <- do.call(cbind, res)
colnames(m) <- sub("_.*", "", basename(colnames(m)))      # GSM id
write.table(round(m, 5), out, sep = "\t", quote = FALSE, col.names = NA)
cat("R", R.version.string, "| affy", as.character(packageVersion("affy")), "| frma",
    as.character(packageVersion("frma")), "| vecs", as.character(packageVersion("hgu133plus2frmavecs")), "\n")
