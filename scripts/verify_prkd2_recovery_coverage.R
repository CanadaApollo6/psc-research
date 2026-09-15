# Independent R calculation of the new marginal evidence-coverage diagnostics.
# Does not load the Python implementation and does not calculate an H4 posterior.
args <- commandArgs(trailingOnly=TRUE)
out <- if (length(args)) args[[1]] else "work/prkd2-recovery-independent-coverage.csv"
g <- read.csv(gzfile("data/derived/prkd2-expanded-inputs/GWAS-PRKD2.csv.gz"), stringsAsFactors=FALSE)
gf <- read.csv(gzfile("data/derived/prkd2-expanded-inputs/GWAS-full-evidence.csv.gz"), stringsAsFactors=FALSE)
v <- 1 / (2 * gf$N * gf$MAF * (1 - gf$MAF) * gf$case_fraction * (1 - gf$case_fraction))
r <- 0.04 / (0.04 + v)
logg <- (log1p(-r) + r * qnorm(gf$pvalue / 2)^2) / 2
norm <- function(x) { e <- exp(x - max(x)); e / sum(e) }
wg <- norm(logg)
files <- c("DICE-original", "BLUEPRINT-original", paste0("IBDverse-original-Myeloid-", c(0,1,3,7)))
results <- lapply(files, function(name) {
  path <- paste0("data/derived/prkd2-recovery-", name, "-alleles.csv.gz")
  if (!file.exists(path)) stop(paste("Missing complete original source", name))
  x <- read.csv(gzfile(path), stringsAsFactors=FALSE, na.strings=c("", "NA"))
  x$in_fixed_region <- toupper(as.character(x$in_fixed_region)) == "TRUE"
  eligible <- x[!is.na(x$snp) & x$in_fixed_region, ]
  common <- intersect(g$snp, eligible$snp)
  gm <- sum(wg[gf$source_row %in% g$source_row[g$snp %in% common]])
  distinct <- x[!duplicated(x[,setdiff(names(x), "source_row"),drop=FALSE]), ]
  qm <- NA_real_
  if (all(is.finite(distinct$source_beta)) && all(is.finite(distinct$source_se)) && all(distinct$source_se > 0)) {
    variance <- distinct$source_se^2
    shrink <- 0.15^2 / (0.15^2 + variance)
    logq <- (log1p(-shrink) + shrink * distinct$source_beta^2 / variance) / 2
    qm <- sum(norm(logq)[!is.na(distinct$snp) & distinct$snp %in% common & distinct$in_fixed_region])
  }
  data.frame(dataset_id=name, shared_edits=length(common), full_PSC_weight=gm, full_RNA_weight=qm,
             source_rows=nrow(x), distinct_source_rows=nrow(distinct))
})
dir.create(dirname(out), recursive=TRUE, showWarnings=FALSE)
write.csv(do.call(rbind, results), out, row.names=FALSE, na="")
cat("Independent R coverage calculations complete\n")

format_results <- lapply(paste0("IBDverse-original-Myeloid-", c(0,1,3,7)), function(name) {
  x <- read.delim(gzfile(paste0("data/derived/prkd2-recovery-queries/", name, "-PRKD2.tsv.gz")))
  statistic <- abs(x$slope / x$slope_se)
  errors <- lapply(1:500, function(df) abs((log(2) + pt(-statistic, df=df, log.p=TRUE) - log(x$pval_nominal))/log(10)))
  chosen <- which.min(vapply(errors, median, numeric(1)))
  normal <- abs((log(2) + pnorm(-statistic, log.p=TRUE) - log(x$pval_nominal))/log(10))
  data.frame(dataset_id=name, best_matching_integer_residual_df=chosen,
             median_abs_log10_P_error=median(errors[[chosen]]),
             p95_abs_log10_P_error=unname(quantile(errors[[chosen]], .95)),
             max_abs_log10_P_error=max(errors[[chosen]]),
             normal_median_abs_log10_P_error=median(normal), normal_max_abs_log10_P_error=max(normal))
})
write.csv(do.call(rbind, format_results), sub("coverage.csv$", "format.csv", out), row.names=FALSE)
cat("Independent R distribution-convention checks complete\n")
