#!/usr/bin/env Rscript

# GSE84161 all-array Affymetrix RMA and predeclared assay QC.
# --check-only validates sources and the isolated runtime without reading CEL
# intensities. --execute is refused until the JSON plan is explicitly frozen.

options(stringsAsFactors = FALSE, warn = 1)

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

first_field <- function(x, names, default = NULL) {
    if (!is.null(x)) for (nm in names)
        if (!is.null(x[[nm]]) && length(x[[nm]]) > 0L) return(x[[nm]])
    default
}
scalar_chr <- function(x, default = NULL) {
    if (is.null(x) || length(x) == 0L) default else as.character(x[[1L]])
}
scalar_int <- function(x, default = NULL) {
    if (is.null(x) || length(x) == 0L) return(default)
    y <- suppressWarnings(as.integer(x[[1L]])); if (is.na(y)) default else y
}
scalar_bool <- function(x, default = NULL) {
    if (is.null(x) || length(x) == 0L) return(default)
    if (is.logical(x)) return(isTRUE(x[[1L]]))
    y <- tolower(trimws(as.character(x[[1L]])))
    if (y %in% c("true", "t", "1", "yes")) TRUE
    else if (y %in% c("false", "f", "0", "no")) FALSE else default
}
stop_usage <- function(status = 2L) {
    cat("Usage: normalize_gse84161.R --plan PLAN.json --check-only|--execute\n")
    quit(save = "no", status = status, runLast = FALSE)
}
parse_args <- function(argv) {
    if (!length(argv)) stop_usage()
    out <- list(plan = NULL, mode = NULL); i <- 1L
    while (i <= length(argv)) {
        a <- argv[[i]]
        if (a %in% c("-h", "--help")) stop_usage(0L)
        if (a == "--plan") {
            if (i == length(argv)) stopf("--plan requires a path")
            i <- i + 1L; out$plan <- argv[[i]]
        } else if (a %in% c("--check-only", "--execute")) {
            if (!is.null(out$mode)) stopf("choose exactly one execution mode")
            out$mode <- if (a == "--execute") "execute" else "check-only"
        } else stopf("unknown argument: %s", a)
        i <- i + 1L
    }
    if (is.null(out$plan) || is.null(out$mode)) stopf("--plan and one execution mode are required")
    out
}

sha256_file <- function(path) {
    if (!file.exists(path)) stopf("file does not exist: %s", path)
    if (!nzchar(Sys.which("sha256sum"))) stopf("sha256sum executable is required")
    out <- system2("sha256sum", c(path), stdout = TRUE, stderr = TRUE)
    if (!length(out)) stopf("sha256sum returned no output for %s", path)
    token <- strsplit(trimws(out[[1L]]), "[[:space:]]+")[[1L]][[1L]]
    if (!grepl("^[[:xdigit:]]{64}$", token)) stopf("invalid sha256sum output for %s", path)
    tolower(token)
}
file_bytes <- function(path) {
    n <- file.info(path)$size
    if (is.na(n)) stopf("cannot stat file: %s", path)
    as.numeric(n)
}
resolve_path <- function(path, repo_root) {
    path <- scalar_chr(path)
    if (is.null(path) || !nzchar(path)) return(NULL)
    if (grepl("^/", path)) normalizePath(path, mustWork = FALSE)
    else normalizePath(file.path(repo_root, path), mustWork = FALSE)
}
read_plan <- function(path) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stopf("jsonlite is required")
    tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE),
             error = function(e) stopf("cannot parse plan: %s", e$message))
}
write_json <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE,
                         null = "null", digits = 17)
}

qsummary <- function(x, probs) {
    q <- stats::quantile(as.numeric(x), probs = probs, type = 7,
                         names = FALSE, na.rm = FALSE)
    as.list(stats::setNames(as.numeric(q), c("q0", "q01", "q25", "q50", "q75", "q99", "q100")))
}
column_summary <- function(x, sample_ids, probs, raw = FALSE) {
    if (!is.matrix(x)) x <- as.matrix(x)
    rows <- lapply(seq_len(ncol(x)), function(i) {
        v <- as.numeric(x[, i]); missing <- sum(is.na(v))
        z <- c(sample_accession = sample_ids[[i]], n_total = length(v),
               n_missing = missing,
               n_nonfinite = sum(!is.finite(v) & !is.na(v)),
               n_nonfinite_including_missing = sum(!is.finite(v)), qsummary(v, probs))
        if (raw) z <- c(z, n_negative = sum(v < 0, na.rm = TRUE),
                        n_zero = sum(v == 0, na.rm = TRUE))
        z
    })
    out <- as.data.frame(do.call(rbind, rows), check.names = FALSE,
                          stringsAsFactors = FALSE)
    out$sample_accession <- as.character(out$sample_accession)
    numeric_names <- setdiff(names(out), "sample_accession")
    out[numeric_names] <- lapply(out[numeric_names], function(y) as.numeric(as.character(y)))
    out
}

validate_runtime <- function(plan, repo_root) {
    runtime <- first_field(plan, c("runtime"), list())
    expected_r <- scalar_chr(first_field(runtime, c("r_version", "R_version")))
    actual_r <- as.character(getRversion())
    if (!is.null(expected_r) && !identical(expected_r, actual_r))
        stopf("R version mismatch: expected %s, running %s", expected_r, actual_r)
    lock_path <- resolve_path(first_field(runtime, c("lock_path")), repo_root)
    lock_sha <- scalar_chr(first_field(runtime, c("lock_sha256")))
    lock <- list(path = lock_path, expected_sha256 = lock_sha,
                 actual_sha256 = NULL, bytes = NULL, status = "not_pinned")
    if (!is.null(lock_path)) {
        if (!file.exists(lock_path)) stopf("runtime lock does not exist: %s", lock_path)
        lock$bytes <- file_bytes(lock_path); lock$actual_sha256 <- sha256_file(lock_path)
        lock$status <- "verified"
        if (!is.null(lock_sha) && !identical(lock$actual_sha256, tolower(lock_sha)))
            stopf("runtime lock SHA-256 mismatch: expected %s, got %s", lock_sha, lock$actual_sha256)
    }
    entries <- first_field(runtime, c("packages"), list()); versions <- list()
    for (entry in entries) {
        pkg <- scalar_chr(first_field(entry, c("r_package", "package", "name")))
        if (is.null(pkg) || !nzchar(pkg)) stopf("runtime package entry has no R package")
        if (!requireNamespace(pkg, quietly = TRUE)) stopf("required R package is unavailable: %s", pkg)
        actual <- as.character(utils::packageVersion(pkg)); expected <- scalar_chr(first_field(entry, c("version")))
        if (!is.null(expected) && !identical(actual, expected))
            stopf("R package %s is %s; expected %s", pkg, actual, expected)
        versions[[pkg]] <- actual
    }
    list(R_version = actual_r, platform = R.version$platform,
         bioconductor_release = scalar_chr(first_field(runtime, c("bioconductor_release"))),
         packages = versions, lock = lock)
}

validate_method <- function(plan) {
    method <- first_field(plan, c("normalization", "rma"), list())
    bools <- c("background", "normalize", "all_arrays_together")
    for (nm in bools) if (!isTRUE(scalar_bool(method[[nm]]))) stopf("method %s must be true", nm)
    if (!identical(scalar_int(method$bgversion), 2L)) stopf("method bgversion must be 2")
    if (!identical(scalar_int(method$expected_sample_count), 30L)) stopf("method expected_sample_count must be 30")
    if (!identical(scalar_int(method$expected_probe_rows), 54675L)) stopf("method expected_probe_rows must be 54675")
    if (!identical(scalar_int(method$quantile_type), 7L)) stopf("method quantile_type must be 7")
    if (!identical(scalar_chr(method$cdf_package), "hgu133plus2cdf")) stopf("method CDF must be hgu133plus2cdf")
    if (isTRUE(scalar_bool(method$destructive, FALSE))) stopf("method destructive must be false")
    probs <- as.numeric(unlist(method$quantile_probabilities %||% list(), use.names = FALSE))
    if (!identical(probs, c(0, .01, .25, .5, .75, .99, 1))) stopf("quantile probabilities are not the frozen set")
    list(background = TRUE, normalize = TRUE, bgversion = 2L,
         destructive = scalar_bool(method$destructive, FALSE), cdf_package = "hgu133plus2cdf",
         all_arrays_together = TRUE, expected_sample_count = 30L,
         expected_probe_rows = 54675L, quantile_probabilities = probs, quantile_type = 7L)
}

load_samples <- function(plan, repo_root, expected_count) {
    cfg <- first_field(plan, c("sample_metadata", "samples"), list())
    path <- resolve_path(first_field(cfg, c("path", "metadata_path")), repo_root)
    if (is.null(path) || !file.exists(path)) stopf("sample metadata is missing")
    actual_sha <- sha256_file(path); expected_sha <- scalar_chr(first_field(cfg, c("sha256", "metadata_sha256")))
    if (!is.null(expected_sha) && !identical(actual_sha, tolower(expected_sha)))
        stopf("sample metadata SHA-256 mismatch")
    tab <- tryCatch(utils::read.csv(path, check.names = FALSE, colClasses = "character",
                                    na.strings = character()),
                    error = function(e) stopf("cannot read sample metadata: %s", e$message))
    acc_col <- scalar_chr(first_field(cfg, c("accession_column", "sample_column")), "sample_accession")
    ord_col <- scalar_chr(first_field(cfg, c("order_column")), "source_order")
    if (!all(c(acc_col, ord_col) %in% names(tab))) stopf("sample metadata lacks accession/order columns")
    acc <- as.character(tab[[acc_col]]); ord <- suppressWarnings(as.integer(tab[[ord_col]]))
    if (length(acc) != expected_count || anyNA(acc) || any(!nzchar(acc)) || anyDuplicated(acc))
        stopf("sample metadata must contain %d unique accessions", expected_count)
    if (anyNA(ord) || anyDuplicated(ord)) stopf("sample metadata order is missing or duplicated")
    idx <- order(ord, method = "radix"); acc <- acc[idx]; tab <- tab[idx, , drop = FALSE]; ord <- ord[idx]
    if (!identical(ord, seq_len(expected_count))) stopf("sample order must be the complete 1..%d sequence", expected_count)
    planned <- first_field(plan, c("sample_order"), first_field(cfg, c("sample_order")))
    if (!is.null(planned) && !identical(as.character(unlist(planned, use.names = FALSE)), acc))
        stopf("pinned sample order differs from metadata")
    list(path = path, sha256 = actual_sha, table = tab, accession = acc, order = ord,
         accession_column = acc_col, order_column = ord_col,
         archive_member_column = scalar_chr(first_field(cfg, c("archive_member_column"))))
}

validate_sources <- function(plan, repo_root, sample_info) {
    cfg <- first_field(plan, c("source_archive", "raw_archive"), list())
    path <- resolve_path(first_field(cfg, c("path", "archive_path")), repo_root)
    if (is.null(path) || !file.exists(path)) stopf("source archive is missing")
    sha <- sha256_file(path); expected_sha <- scalar_chr(first_field(cfg, c("sha256", "archive_sha256")))
    if (!is.null(expected_sha) && !identical(sha, tolower(expected_sha))) stopf("source archive SHA-256 mismatch")
    bytes <- file_bytes(path); expected_bytes <- scalar_int(first_field(cfg, c("bytes", "archive_bytes")))
    if (!is.null(expected_bytes) && !identical(bytes, as.numeric(expected_bytes))) stopf("source archive byte count mismatch")
    inv_path <- resolve_path(first_field(cfg, c("inventory_path", "raw_inventory_path")), repo_root)
    inv <- NULL; inv_sha <- NULL; inv_idx <- NULL
    if (!is.null(inv_path)) {
        if (!file.exists(inv_path)) stopf("raw inventory is missing")
        inv_sha <- sha256_file(inv_path); inv_expected <- scalar_chr(first_field(cfg, c("inventory_sha256", "raw_inventory_sha256")))
        if (!is.null(inv_expected) && !identical(inv_sha, tolower(inv_expected))) stopf("raw inventory SHA-256 mismatch")
        inv <- utils::read.csv(inv_path, check.names = FALSE, colClasses = "character", na.strings = character())
        if (!all(c("gsm", "archive_member") %in% names(inv)) || anyDuplicated(inv$gsm) || anyDuplicated(inv$archive_member))
            stopf("raw inventory must have unique gsm and archive_member columns")
        inv_idx <- match(sample_info$accession, inv$gsm); if (anyNA(inv_idx)) stopf("raw inventory lacks a sample")
    }
    meta_member <- NULL
    if (!is.null(sample_info$archive_member_column) && sample_info$archive_member_column %in% names(sample_info$table))
        meta_member <- as.character(sample_info$table[[sample_info$archive_member_column]])
    expected_members <- if (!is.null(inv)) as.character(inv$archive_member[inv_idx]) else meta_member
    if (is.null(expected_members)) {
        cfg_samples <- first_field(plan, c("sample_metadata", "samples"), list()); col <- scalar_chr(first_field(cfg_samples, c("cel_url_column")), "cel_url")
        if (!(col %in% names(sample_info$table))) stopf("pin raw inventory or archive members")
        expected_members <- sub("^.*[/]", "", as.character(sample_info$table[[col]]))
    }
    if (length(expected_members) != length(sample_info$accession) || anyNA(expected_members) || anyDuplicated(expected_members) ||
        any(!grepl("[.]CEL[.]gz$", expected_members, ignore.case = TRUE))) stopf("invalid expected CEL member list")
    archive_members <- tryCatch(as.character(utils::untar(path, list = TRUE)),
                                error = function(e) stopf("cannot list archive: %s", e$message))
    archive_members <- archive_members[nzchar(archive_members)]
    expected_member_count <- scalar_int(first_field(cfg, c("expected_members", "member_count")), length(expected_members))
    if (!identical(length(expected_members), expected_member_count)) stopf("expected archive member count is not %d", expected_member_count)
    if (anyDuplicated(archive_members) || !setequal(archive_members, expected_members) || length(archive_members) != length(expected_members))
        stopf("source archive members do not match the 30 metadata CEL files")
    cel_dir <- resolve_path(first_field(cfg, c("cel_directory", "extracted_cel_directory")), repo_root)
    cel_paths <- file.path(cel_dir, expected_members); present <- file.exists(cel_paths); hashes <- list()
    if (all(present)) for (i in seq_along(cel_paths)) {
        if (!is.null(inv) && "compressed_bytes" %in% names(inv)) {
            n <- suppressWarnings(as.numeric(inv$compressed_bytes[inv_idx[[i]]]))
            if (!is.na(n) && !identical(file_bytes(cel_paths[[i]]), n)) stopf("extracted CEL byte mismatch: %s", sample_info$accession[[i]])
        }
        got <- sha256_file(cel_paths[[i]])
        if (!is.null(inv) && "compressed_sha256" %in% names(inv)) {
            want <- tolower(as.character(inv$compressed_sha256[inv_idx[[i]]]))
            if (nzchar(want) && !identical(got, want)) stopf("extracted CEL SHA-256 mismatch: %s", sample_info$accession[[i]])
        }
        hashes[[sample_info$accession[[i]]]] <- list(path = cel_paths[[i]], bytes = file_bytes(cel_paths[[i]]), sha256 = got)
    }
    list(archive_path = path, archive_url = scalar_chr(first_field(cfg, c("url", "source_url"))),
         archive_bytes = bytes, archive_sha256 = sha, archive_members = expected_members,
         archive_member_count = length(archive_members), inventory_path = inv_path, inventory_sha256 = inv_sha,
         cel_directory = cel_dir, extracted_paths = cel_paths, extracted_present = present,
         extracted_all_present = all(present), extracted_hashes = hashes)
}

execution_gate <- function(plan) {
    freeze <- first_field(plan, c("execution_freeze", "normalization_execution"), list())
    status <- scalar_chr(first_field(plan, c("normalization_run_status", "execution_status")), NULL)
    if (is.null(status)) status <- scalar_chr(first_field(freeze, c("status")), NULL)
    allowed <- scalar_bool(first_field(plan, c("normalization_run_allowed", "execution_allowed")), NULL)
    if (is.null(allowed)) allowed <- scalar_bool(first_field(freeze, c("normalization_run_allowed", "allowed")), FALSE)
    list(status = status, allowed = isTRUE(allowed), frozen = identical(status, "frozen"))
}
out_path <- function(outputs, key, default, repo_root) resolve_path(first_field(outputs, c(key), default), repo_root)

validate_input_hashes <- function(plan, repo_root) {
    entries <- first_field(plan, c("input_sha256", "inputs_sha256"), NULL)
    if (is.null(entries) || !length(entries)) return(list(status = "not_pinned", checked = list()))
    checked <- list()
    for (key in names(entries)) {
        entry <- entries[[key]]
        path_text <- key; expected <- entry
        if (is.list(entry)) {
            path_text <- scalar_chr(first_field(entry, c("path", "local_path")), key)
            expected <- first_field(entry, c("sha256", "expected_sha256"))
        }
        expected <- scalar_chr(expected)
        if (is.null(expected)) stopf("input hash entry has no SHA-256: %s", key)
        path <- resolve_path(path_text, repo_root)
        got <- sha256_file(path)
        if (!identical(got, tolower(expected)))
            stopf("input SHA-256 mismatch for %s: expected %s, got %s", path_text, expected, got)
        checked[[path_text]] <- list(path = path, sha256 = got, bytes = file_bytes(path))
    }
    list(status = "verified", checked = checked)
}

aggregate_entrez <- function(expr, plan, repo_root, ignored_dir, qc_dir) {
    mapping <- first_field(plan, c("mapping"), NULL)
    if (is.null(mapping)) return(list(status = "not_run_probe_map_not_pinned"))
    map_path <- resolve_path(first_field(mapping, c("probe_map_path", "path")), repo_root)
    if (is.null(map_path)) return(list(status = "not_run_probe_map_not_pinned"))
    if (!file.exists(map_path)) stopf("pinned probe map does not exist: %s", map_path)
    map_sha <- sha256_file(map_path); expected_sha <- scalar_chr(first_field(mapping, c("sha256", "probe_map_sha256")))
    if (!is.null(expected_sha) && !identical(map_sha, tolower(expected_sha))) stopf("probe map SHA-256 mismatch")
    probe_col <- scalar_chr(first_field(mapping, c("probe_id_column", "probe_set_column")), "probe_set_id")
    entrez_col <- scalar_chr(first_field(mapping, c("entrez_id_column", "entrez_column")), "entrez_gene_id")
    eligible_col <- scalar_chr(first_field(mapping, c("eligible_column")), "eligible_single_entrez")
    tab <- utils::read.csv(map_path, check.names = FALSE, colClasses = "character", na.strings = character())
    if (!all(c(probe_col, entrez_col, eligible_col) %in% names(tab)))
        stopf("probe map must contain %s, %s and %s", probe_col, entrez_col, eligible_col)
    probe_ids <- as.character(tab[[probe_col]]); entrez_ids <- trimws(as.character(tab[[entrez_col]]))
    eligible_raw <- tab[[eligible_col]]
    eligible <- if (is.logical(eligible_raw)) eligible_raw else {
        y <- tolower(trimws(as.character(eligible_raw))); z <- rep(NA, length(y))
        z[y %in% c("true", "t", "1", "yes")] <- TRUE; z[y %in% c("false", "f", "0", "no")] <- FALSE
        if (anyNA(z)) stopf("probe map eligible column is not boolean"); z
    }
    if (anyNA(probe_ids) || any(!nzchar(probe_ids)) || anyDuplicated(probe_ids)) stopf("probe map IDs are missing or duplicated")
    if (any(eligible & (is.na(entrez_ids) | !nzchar(entrez_ids) | grepl("///", entrez_ids, fixed = TRUE))))
        stopf("eligible probe map rows must contain one nonblank Entrez ID")
    expected_rows <- scalar_int(first_field(mapping, c("expected_probe_rows")), nrow(expr))
    if (nrow(tab) != expected_rows) stopf("probe map has %d rows; expected %d", nrow(tab), expected_rows)
    if (!setequal(probe_ids, rownames(expr))) stopf("probe map IDs do not match normalized probe IDs")
    idx <- match(rownames(expr), probe_ids); eligible <- eligible[idx]; entrez_ids <- entrez_ids[idx]
    genes <- sort(unique(entrez_ids[eligible]), method = "radix")
    expected_genes <- scalar_int(first_field(mapping, c("expected_entrez_count", "expected_gene_count")), 20486L)
    if (length(genes) != expected_genes) stopf("probe map represents %d Entrez genes; expected %d", length(genes), expected_genes)
    expected_eligible <- scalar_int(first_field(mapping, c("expected_eligible_probe_count")), NULL)
    if (!is.null(expected_eligible) && sum(eligible) != expected_eligible)
        stopf("probe map has %d eligible probes; expected %d", sum(eligible), expected_eligible)
    gene_matrix <- matrix(NA_real_, nrow = length(genes), ncol = ncol(expr),
                          dimnames = list(genes, colnames(expr)))
    probe_counts <- integer(length(genes))
    for (i in seq_along(genes)) {
        probes <- which(eligible & entrez_ids == genes[[i]])
        probe_counts[[i]] <- length(probes)
        gene_matrix[i, ] <- apply(expr[probes, , drop = FALSE], 2L, stats::median, na.rm = FALSE)
    }
    if (any(!is.finite(gene_matrix))) stopf("Entrez median matrix contains non-finite values")
    matrix_path <- file.path(ignored_dir, "normalized_entrez_median.rds")
    saveRDS(gene_matrix, matrix_path, compress = "xz")
    counts_path <- file.path(qc_dir, "entrez-probe-counts.csv")
    utils::write.csv(data.frame(entrez_gene_id = genes, contributing_probe_count = probe_counts),
                     counts_path, row.names = FALSE, quote = TRUE)
    list(status = "completed", probe_map_path = map_path, probe_map_sha256 = map_sha,
         eligible_probe_count = sum(eligible), entrez_gene_count = length(genes), dimensions = dim(gene_matrix),
         matrix_path = matrix_path, matrix_sha256 = sha256_file(matrix_path),
         contributing_probe_counts_path = counts_path, contributing_probe_counts_sha256 = sha256_file(counts_path),
         feature_order = "lexicographic Entrez Gene ID",
         aggregation = "median of all eligible single-Entrez normalized probe sets per array",
         ensembl_crosswalk_used = FALSE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
plan_path <- normalizePath(args$plan, mustWork = TRUE)
plan <- read_plan(plan_path)
plan_dir <- dirname(plan_path)
repo_hint <- scalar_chr(first_field(plan, c("repository_root", "repo_root")), "..")
repo_root <- if (grepl("^/", repo_hint)) normalizePath(repo_hint, mustWork = TRUE) else
    normalizePath(file.path(plan_dir, repo_hint), mustWork = TRUE)
runtime_info <- validate_runtime(plan, repo_root)
method_info <- validate_method(plan)
input_hash_info <- validate_input_hashes(plan, repo_root)
sample_info <- load_samples(plan, repo_root, method_info$expected_sample_count)
source_info <- validate_sources(plan, repo_root, sample_info)
gate <- execution_gate(plan)
check_result <- list(
    checked_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), mode = args$mode,
    plan_path = plan_path, plan_sha256 = sha256_file(plan_path), repository_root = repo_root,
    runtime = runtime_info, method = method_info, input_hashes = input_hash_info,
    samples = list(count = length(sample_info$accession), order_column = sample_info$order_column,
                   accessions = sample_info$accession, metadata_path = sample_info$path,
                   metadata_sha256 = sample_info$sha256),
    source = source_info, execution_gate = gate,
    expression_values_read = FALSE, expression_normalization_run = FALSE)
if (identical(args$mode, "check-only")) {
    cat(jsonlite::toJSON(check_result, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 17), "\n")
    quit(save = "no", status = 0L, runLast = FALSE)
}
if (!gate$frozen || !gate$allowed)
    stopf("refusing --execute: plan must set normalization_run_status=\"frozen\" and normalization_run_allowed=true")
if (!source_info$extracted_all_present)
    stopf("refusing --execute: all 30 extracted CEL files must be present")

outputs <- first_field(plan, c("outputs"), list())
ignored_dir <- out_path(outputs, "ignored_output_dir", "work/gse84161-normalization/output", repo_root)
qc_dir <- out_path(outputs, "versioned_qc_dir", "data/derived/gse84161-array-qc", repo_root)
dir.create(ignored_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
sample_ids <- sample_info$accession; cel_paths <- source_info$extracted_paths
names(cel_paths) <- sample_ids; probs <- method_info$quantile_probabilities
started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)

suppressPackageStartupMessages({
    if (!requireNamespace("affy", quietly = TRUE)) stopf("affy is unavailable")
    if (!requireNamespace("Biobase", quietly = TRUE)) stopf("Biobase is unavailable")
    if (!requireNamespace(method_info$cdf_package, quietly = TRUE)) stopf("CDF package is unavailable")
})
raw_ab <- affy::ReadAffy(filenames = unname(cel_paths), celfile.path = "",
                         sampleNames = sample_ids, cdfname = method_info$cdf_package, verbose = TRUE)
if (!identical(as.character(Biobase::sampleNames(raw_ab)), sample_ids)) stopf("ReadAffy sample order mismatch")
pm_matrix <- as.matrix(affy::pm(raw_ab)); if (ncol(pm_matrix) != length(sample_ids)) stopf("wrong PM sample count")
if (!is.null(colnames(pm_matrix)) && !identical(as.character(colnames(pm_matrix)), sample_ids)) stopf("PM sample order mismatch")
colnames(pm_matrix) <- sample_ids
pm_missing <- sum(is.na(pm_matrix)); pm_nonfinite <- sum(!is.finite(pm_matrix) & !is.na(pm_matrix)); pm_negative <- sum(pm_matrix < 0, na.rm = TRUE)
if (pm_missing > 0L || pm_nonfinite > 0L || pm_negative > 0L)
    stopf("raw CDF-defined PM entries failed finite nonnegative QC (missing=%d, nonfinite=%d, negative=%d)", pm_missing, pm_nonfinite, pm_negative)
raw_pm <- column_summary(pm_matrix, sample_ids, probs, raw = TRUE)
raw_pm_path <- out_path(outputs, "raw_pm_summary_csv", "data/derived/gse84161-array-qc/raw-pm-summary.csv", repo_root)
utils::write.csv(raw_pm, raw_pm_path, row.names = FALSE, quote = TRUE)

expr_eset <- affy::rma(raw_ab, background = TRUE, normalize = TRUE, bgversion = 2,
                       destructive = method_info$destructive, verbose = TRUE)
expr <- as.matrix(Biobase::exprs(expr_eset))
if (!identical(dim(expr), c(method_info$expected_probe_rows, length(sample_ids))))
    stopf("normalized matrix dimensions are %s; expected %d x %d", paste(dim(expr), collapse = " x "), method_info$expected_probe_rows, length(sample_ids))
if (is.null(rownames(expr)) || anyNA(rownames(expr)) || any(!nzchar(rownames(expr))) || anyDuplicated(rownames(expr))) stopf("probe IDs are missing or duplicated")
if (is.null(colnames(expr)) || !identical(as.character(colnames(expr)), sample_ids)) stopf("normalized sample order mismatch")
if (any(!is.finite(expr))) stopf("normalized matrix contains non-finite values")
normalized <- column_summary(expr, sample_ids, probs, raw = FALSE)
normalized_path <- out_path(outputs, "normalized_probe_summary_csv", "data/derived/gse84161-array-qc/normalized-probe-summary.csv", repo_root)
utils::write.csv(normalized, normalized_path, row.names = FALSE, quote = TRUE)
probe_medians <- apply(expr, 1L, stats::median, na.rm = FALSE); if (any(!is.finite(probe_medians))) stopf("RMA probe medians are non-finite")
rle_matrix <- sweep(expr, 1L, probe_medians, FUN = "-"); if (any(!is.finite(rle_matrix))) stopf("RLE matrix contains non-finite values")
rle <- column_summary(rle_matrix, sample_ids, probs, raw = FALSE); rle$rle_iqr <- rle$q75 - rle$q25
rle_path <- out_path(outputs, "rle_summary_csv", "data/derived/gse84161-array-qc/rle-summary.csv", repo_root)
utils::write.csv(rle, rle_path, row.names = FALSE, quote = TRUE)
sample_order_path <- out_path(outputs, "sample_order_csv", "data/derived/gse84161-array-qc/sample-order.csv", repo_root)
utils::write.csv(data.frame(source_order = sample_info$order, sample_accession = sample_ids, archive_member = source_info$archive_members),
                 sample_order_path, row.names = FALSE, quote = TRUE)

signal_plot_path <- out_path(outputs, "signal_distribution_plot", "data/derived/gse84161-array-qc/signal-distributions.pdf", repo_root)
rle_plot_path <- out_path(outputs, "rle_plot", "data/derived/gse84161-array-qc/rle-boxplots.pdf", repo_root)
dir.create(dirname(signal_plot_path), recursive = TRUE, showWarnings = FALSE)
grDevices::pdf(signal_plot_path, width = 13, height = 9); old_par <- graphics::par(no.readonly = TRUE)
graphics::par(mfrow = c(2, 2), mar = c(9, 4, 3, 1))
graphics::boxplot(log2(pm_matrix + 1), las = 2, main = "Raw CDF PM log2(PM + 1)", ylab = "log2(PM + 1)", cex.axis = .65)
graphics::boxplot(expr, las = 2, main = "RMA probe-set signal", ylab = "log2 RMA", cex.axis = .65)
graphics::hist(as.numeric(pm_matrix), breaks = 100, main = "Raw PM signal distribution", xlab = "PM intensity", col = "grey80", border = "white")
graphics::hist(as.numeric(expr), breaks = 100, main = "RMA signal distribution", xlab = "log2 RMA", col = "grey80", border = "white")
graphics::par(old_par); grDevices::dev.off()
grDevices::pdf(rle_plot_path, width = 13, height = 7); graphics::par(mar = c(9, 4, 3, 1))
graphics::boxplot(rle_matrix, las = 2, main = "Relative log expression", ylab = "RMA - per-probe across-array median", cex.axis = .65)
grDevices::dev.off()

probe_rds <- file.path(ignored_dir, "normalized_probe_rma.rds")
expression_set_rds <- file.path(ignored_dir, "normalized_probe_rma_expression_set.rds")
saveRDS(expr, probe_rds, compress = "xz"); saveRDS(expr_eset, expression_set_rds, compress = "xz")
mapping_info <- aggregate_entrez(expr, plan, repo_root, ignored_dir, qc_dir)
output_files <- unique(c(raw_pm_path, normalized_path, rle_path, sample_order_path, signal_plot_path, rle_plot_path, probe_rds, expression_set_rds))
if (identical(mapping_info$status, "completed")) output_files <- c(output_files, mapping_info$matrix_path, mapping_info$contributing_probe_counts_path)
output_manifest <- lapply(output_files, function(path) list(path = path, bytes = file_bytes(path), sha256 = sha256_file(path)))
names(output_manifest) <- basename(output_files)
audit <- utils::modifyList(check_result, list(
    started_at_utc = started_at, finished_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    mode = "execute", expression_values_read = TRUE, expression_normalization_run = TRUE,
    raw_pm = list(definition = "CDF-defined perfect-match entries returned by affy::pm",
                  dimensions = dim(pm_matrix), missing_count = pm_missing,
                  nonfinite_count_excluding_missing = pm_nonfinite, negative_count = pm_negative,
                  zero_count = sum(pm_matrix == 0), hard_stop_passed = TRUE,
                  quantile_type = 7L, quantile_probabilities = probs,
                  summary_path = raw_pm_path, summary_sha256 = sha256_file(raw_pm_path)),
    normalized_probe_rma = list(dimensions = dim(expr), finite = all(is.finite(expr)),
                                quantity = "log2 RMA probe-set expression", cdf_package = method_info$cdf_package,
                                background = TRUE, normalize = TRUE, bgversion = 2L,
                                destructive = method_info$destructive, summary_path = normalized_path,
                                summary_sha256 = sha256_file(normalized_path), probe_matrix_path = probe_rds,
                                probe_matrix_sha256 = sha256_file(probe_rds), expression_set_path = expression_set_rds,
                                expression_set_sha256 = sha256_file(expression_set_rds)),
    rle = list(definition = "full-probe RMA matrix minus each probe's median across all 30 arrays",
               dimensions = dim(rle_matrix), finite = all(is.finite(rle_matrix)), summary_path = rle_path,
               summary_sha256 = sha256_file(rle_path), plot_path = rle_plot_path,
               plot_sha256 = sha256_file(rle_plot_path)),
    plots = list(path = signal_plot_path, sha256 = sha256_file(signal_plot_path)),
    outputs = output_manifest, mapping = mapping_info,
    automatic_signal_based_exclusions = FALSE, program_scores_computed = FALSE,
    treatment_or_donor_contrasts_computed = FALSE
))
audit_path <- out_path(outputs, "audit_json", "data/derived/gse84161-array-qc/normalization-audit.json", repo_root)
write_json(audit, audit_path)
cat(jsonlite::toJSON(audit, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 17), "\n")
