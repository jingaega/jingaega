# In-fold limma (PLAN §4 step 1+3). Reads a job list; each job = training GSM ids + mode.
# Low-expression filter and moderated F are fitted on the job's training samples only.
# Usage: Rscript limma_batch.R <expr_tsv (genes x samples)> <samples_tsv> <jobs_tsv> <out_dir>
suppressMessages(library(limma))
a <- commandArgs(trailingOnly = TRUE)
expr <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
smp <- read.delim(a[2], row.names = 1)
jobs <- read.delim(a[3], stringsAsFactors = FALSE, colClasses = "character")   # key, mode, gsm, labels (comma-separated; empty = true labels)
for (j in seq_len(nrow(jobs))) {
  out <- file.path(a[4], paste0(jobs$key[j], ".tsv"))
  if (file.exists(out)) next
  g <- strsplit(jobs$gsm[j], ",")[[1]]
  x <- expr[, g, drop = FALSE]; s <- smp[g, ]
  if (!is.na(jobs$labels[j]) && nchar(jobs$labels[j]) > 0) s$diagnosis <- strsplit(jobs$labels[j], ",")[[1]]
  mu <- rowMeans(x); x <- x[mu > quantile(mu, 0.25), ]
  dx <- factor(s$diagnosis); reg <- factor(s$region)
  if (jobs$mode[j] == "cov") {
    design <- model.matrix(~ 0 + dx + reg + ph + rin + pmi + age + sex_expr, data = s)
  } else {
    design <- model.matrix(~ 0 + dx + reg)
  }
  colnames(design) <- make.names(colnames(design))
  dc <- duplicateCorrelation(x, design, block = s$donor)
  fit <- lmFit(x, design, block = s$donor, correlation = dc$consensus.correlation)
  lv <- paste0("dx", levels(dx))
  cm <- makeContrasts(contrasts = paste(lv[-1], "-", lv[1]), levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cm))
  tt <- topTableF(fit2, number = Inf, sort.by = "none")
  res <- data.frame(gene = rownames(tt), F = tt$F, P = tt$P.Value, train_mean = mu[rownames(tt)])
  if (ncol(cm) == 1) res$logFC <- fit2$coefficients[, 1]
  attr_line <- sprintf("# consensus_cor=%.4f limma=%s n=%d", dc$consensus.correlation,
                       as.character(packageVersion("limma")), length(g))
  writeLines(attr_line, out)
  suppressWarnings(write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
}
