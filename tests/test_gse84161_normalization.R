#!/usr/bin/env Rscript

# Synthetic regression fixture for the normalization helpers. It evaluates
# only function definitions from the runner; no CEL file is opened.

options(stringsAsFactors = FALSE, warn = 1)
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- all_args[grepl("^--file=", all_args)][1L]
if (is.na(file_arg)) stop("test script path is unavailable", call. = FALSE)
this_file <- sub("^--file=", "", file_arg)
root <- normalizePath(file.path(dirname(normalizePath(this_file)), ".."), mustWork = TRUE)
script_path <- file.path(root, "scripts", "normalize_gse84161.R")
expressions <- parse(file = script_path)
helper_env <- new.env(parent = globalenv())
for (expression in expressions) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        length(expression) >= 3L && is.call(expression[[3L]]) &&
        identical(expression[[3L]][[1L]], as.name("function"))) {
        eval(expression, envir = helper_env)
    }
}

stop_if <- function(condition, message) if (isTRUE(condition)) stop(message, call. = FALSE)
probs <- c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1)
finite <- matrix(c(0, 1, 2, 4, 5, 6), nrow = 3L,
                 dimnames = list(c("p1", "p2", "p3"), c("GSM1", "GSM2")))
raw <- helper_env$column_summary(finite, c("GSM1", "GSM2"), probs, raw = TRUE)
stop_if(!is.character(raw$sample_accession), "sample_accession must be an atomic character column")
stop_if(!identical(raw$sample_accession, c("GSM1", "GSM2")), "sample IDs changed during summary")
stop_if(raw$n_zero[1L] != 1 || raw$n_zero[2L] != 0, "zero counts are incorrect")
stop_if(raw$n_missing[1L] != 0 || raw$n_nonfinite[1L] != 0, "finite counts are incorrect")
stop_if(abs(raw$q50[1L] - 1) > 1e-12 || abs(raw$q50[2L] - 5) > 1e-12,
        "raw type-7 medians are incorrect")

temporary <- tempfile("gse84161-summary-")
dir.create(temporary)
raw_csv <- file.path(temporary, "raw.csv")
utils::write.csv(raw, raw_csv, row.names = FALSE)
round_trip <- utils::read.csv(raw_csv, check.names = FALSE)
stop_if(!identical(as.character(round_trip$sample_accession), c("GSM1", "GSM2")),
        "summary CSV did not preserve sample IDs")

with_invalid <- c(NA_real_, Inf, 6)
stop_if(sum(is.na(with_invalid)) != 1L || sum(!is.finite(with_invalid) & !is.na(with_invalid)) != 1L,
        "missing/nonfinite count definitions are incorrect")

probe_medians <- apply(finite, 1L, stats::median)
rle <- sweep(finite, 1L, probe_medians, FUN = "-")
rle_summary <- helper_env$column_summary(rle, c("GSM1", "GSM2"), probs, raw = FALSE)
stop_if(abs(rle_summary$q0[1L] + 2) > 1e-12 || abs(rle_summary$q50[1L] + 2) > 1e-12 ||
        abs(rle_summary$q100[1L] + 2) > 1e-12, "RLE type-7 quantiles are incorrect")

map_path <- file.path(temporary, "probe-map.csv")
utils::write.csv(data.frame(
    probe_id = c("p1", "p2", "p3"),
    entrez_gene_id = c("10", "10", "11"),
    eligible_probe = c("true", "true", "true")
), map_path, row.names = FALSE, quote = TRUE)
ignored_dir <- file.path(temporary, "ignored")
qc_dir <- file.path(temporary, "qc")
dir.create(ignored_dir)
dir.create(qc_dir)
gene_result <- helper_env$aggregate_entrez(
    finite, list(mapping = list(
        probe_map_path = map_path,
        probe_id_column = "probe_id",
        entrez_id_column = "entrez_gene_id",
        eligible_column = "eligible_probe",
        expected_probe_rows = 3L,
        expected_eligible_probe_count = 3L,
        expected_entrez_count = 2L
    )), root, ignored_dir, qc_dir)
gene_matrix <- readRDS(gene_result$matrix_path)
stop_if(!identical(rownames(gene_matrix), c("10", "11")), "Entrez feature order changed")
stop_if(abs(gene_matrix["10", "GSM1"] - 0.5) > 1e-12 ||
        abs(gene_matrix["10", "GSM2"] - 4.5) > 1e-12,
        "per-array median aggregation is incorrect")
stop_if(any(!is.finite(gene_matrix)), "synthetic Entrez matrix is not finite")

cat("GSE84161 normalization helper regression checks passed\n")
