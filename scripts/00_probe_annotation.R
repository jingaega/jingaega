# Probe -> gene symbol / chromosome annotation for HG-U133 Plus 2 (fixed annotation, no data fitting).
suppressMessages({library(hgu133plus2.db); library(AnnotationDbi)})
ids <- keys(hgu133plus2.db, keytype = "PROBEID")
sym <- AnnotationDbi::select(hgu133plus2.db, keys = ids, columns = c("SYMBOL", "CHR"), keytype = "PROBEID")
# one row per probe: collapse multiple mappings and flag them
agg <- aggregate(cbind(SYMBOL, CHR) ~ PROBEID, data = sym, na.action = na.pass,
                 FUN = function(x) paste(sort(unique(na.omit(x))), collapse = "|"))
write.table(agg, "data/processed/probe_annot.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("probes", nrow(agg), "hgu133plus2.db", as.character(packageVersion("hgu133plus2.db")), "\n")
