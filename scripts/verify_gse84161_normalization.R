#!/usr/bin/env Rscript

# Independent structural and numerical check for the frozen GSE84161 array run.
# This script deliberately does not source the primary normalisation runner.
# It reads only the frozen plan, completed outputs and versioned QC summaries.

options(stringsAsFactors = FALSE, warn = 1)

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

scalar_chr <- function(x, default = NULL) {
    if (is.null(x) || !length(x)) return(default)
    y <- as.character(unlist(x, use.names = FALSE))
    if (!length(y) || is.na(y[[1L]]) || !nzchar(y[[1L]])) default else y[[1L]]
}

scalar_int <- function(x, default = NULL) {
    y <- suppressWarnings(as.integer(scalar_chr(x, NA_character_)))
    if (is.na(y)) default else y
}

field <- function(x, name, default = NULL) {
    if (!is.null(x) && !is.null(x[[name]])) x[[name]] else default
}

resolve_path <- function(path, root) {
    path <- scalar_chr(path)
    if (is.null(path)) return(NULL)
    normalizePath(if (grepl("^/", path)) path else file.path(root, path), mustWork = FALSE)
}

sha256_file <- function(path) {
    if (!file.exists(path)) stopf("file does not exist: %s", path)
    exe <- Sys.which("sha256sum")
    if (!nzchar(exe)) stopf("sha256sum executable is required")
    out <- suppressWarnings(system2(exe, c("--", path), stdout = TRUE, stderr = TRUE))
    if (!length(out)) stopf("sha256sum returned no output for %s", path)
    token <- strsplit(trimws(out[[1L]]), "[[:space:]]+")[[1L]][[1L]]
    if (!grepl("^[[:xdigit:]]{64}$", token)) stopf("invalid SHA-256 output for %s", path)
    tolower(token)
}

file_record <- function(path) {
    if (!file.exists(path)) stopf("missing output: %s", path)
    info <- file.info(path)
    list(path = path, bytes = as.numeric(info$size), sha256 = sha256_file(path))
}

read_plan <- function(path) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) stopf("jsonlite is required")
    tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE),
             error = function(e) stopf("cannot parse plan: %s", e$message))
}

read_json <- function(path) {
    tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE),
             error = function(e) stopf("cannot parse JSON %s: %s", path, e$message))
}

read_csv_chars <- function(path) {
    tryCatch(utils::read.csv(path, check.names = FALSE, colClasses = "character",
                             na.strings = character()),
             error = function(e) stopf("cannot read CSV %s: %s", path, e$message))
}

read_gpl_probe_ids <- function(path) {
    lines <- tryCatch(readLines(path, warn = FALSE),
                      error = function(e) stopf("cannot read GPL SOFT %s: %s", path, e$message))
    lines <- sub("\\r$", "", lines)
    begin <- which(lines == "!platform_table_begin")
    end <- which(lines == "!platform_table_end")
    if (length(begin) != 1L || length(end) != 1L || end <= begin + 1L)
        stopf("GPL SOFT has no unique platform table")
    header <- strsplit(lines[[begin + 1L]], "\\t", fixed = FALSE)[[1L]]
    if (!identical(header[[1L]], "ID")) stopf("GPL platform table does not start with ID")
    rows <- lines[(begin + 2L):(end - 1L)]
    ids <- sub("\\t.*$", "", rows)
    if (length(ids) != 54675L || any(!nzchar(ids)) || anyDuplicated(ids))
        stopf("GPL platform table has %d non-unique probe IDs; expected 54675", length(ids))
    ids
}

read_stable_rds <- function(path, label) {
    before <- file_record(path)
    value <- tryCatch(readRDS(path), error = function(e) stopf("cannot read %s: %s", label, e$message))
    after <- file_record(path)
    if (!identical(before$sha256, after$sha256) || !identical(before$bytes, after$bytes))
        stopf("%s changed while being read", label)
    list(value = value, artifact = after)
}

read_stable_csv <- function(path, label) {
    before <- file_record(path)
    value <- read_csv_chars(path)
    after <- file_record(path)
    if (!identical(before$sha256, after$sha256) || !identical(before$bytes, after$bytes))
        stopf("%s changed while being read", label)
    list(value = value, artifact = after)
}

check_input_hashes <- function(plan, root) {
    entries <- field(plan, "input_sha256")
    if (is.null(entries) || !length(entries) || is.null(names(entries)))
        stopf("the execution plan has no named input_sha256 lock")
    checked <- list()
    for (key in names(entries)) {
        expected <- tolower(scalar_chr(entries[[key]], ""))
        if (!grepl("^[[:xdigit:]]{64}$", expected)) stopf("invalid pinned SHA-256 for %s", key)
        path <- resolve_path(key, root)
        got <- sha256_file(path)
        if (!identical(got, expected)) stopf("pinned input SHA-256 mismatch for %s", key)
        checked[[key]] <- list(path = path, bytes = file_record(path)$bytes, sha256 = got)
    }
    checked
}

load_samples <- function(plan, root, expected_n) {
    cfg <- field(plan, "sample_metadata", list())
    path <- resolve_path(field(cfg, "path", "data/derived/gse84161-samples.csv"), root)
    expected_sha <- scalar_chr(field(cfg, "sha256"))
    got_sha <- sha256_file(path)
    if (!is.null(expected_sha) && !identical(got_sha, tolower(expected_sha)))
        stopf("sample metadata SHA-256 mismatch")
    tab <- read_csv_chars(path)
    acc_col <- scalar_chr(field(cfg, "accession_column"), "sample_accession")
    order_col <- scalar_chr(field(cfg, "order_column"), "source_order")
    if (!all(c(acc_col, order_col) %in% names(tab))) stopf("sample metadata lacks accession/order columns")
    acc <- as.character(tab[[acc_col]])
    ord <- suppressWarnings(as.integer(tab[[order_col]]))
    if (length(acc) != expected_n || anyNA(acc) || any(!nzchar(acc)) || anyDuplicated(acc))
        stopf("sample metadata does not contain %d unique accessions", expected_n)
    if (anyNA(ord) || anyDuplicated(ord)) stopf("sample metadata order is missing or duplicated")
    ix <- order(ord, method = "radix")
    acc <- acc[ix]; ord <- ord[ix]
    if (!identical(ord, seq_len(expected_n))) stopf("sample metadata order is not 1..%d", expected_n)
    planned <- field(plan, "sample_order")
    if (!is.null(planned) && !identical(as.character(unlist(planned, use.names = FALSE)), acc))
        stopf("plan sample order differs from sample metadata")
    list(path = path, sha256 = got_sha, accession = acc, order = ord)
}

check_matrix <- function(x, label, nrow_expected, sample_ids, row_order = NULL) {
    if (!is.matrix(x) || !is.numeric(x) || !identical(typeof(x), "double"))
        stopf("%s is not a numeric double matrix", label)
    if (!identical(dim(x), c(nrow_expected, length(sample_ids))))
        stopf("%s dimensions are %s; expected %d x %d", label,
              paste(dim(x), collapse = " x "), nrow_expected, length(sample_ids))
    ids <- rownames(x); cols <- colnames(x)
    if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
        stopf("%s row IDs are missing or duplicated", label)
    if (is.null(cols) || !identical(as.character(cols), sample_ids))
        stopf("%s sample columns do not match the frozen sample order", label)
    if (any(!is.finite(x))) stopf("%s contains non-finite values", label)
    if (!is.null(row_order) && !identical(as.character(ids), as.character(row_order)))
        stopf("%s row IDs are not in the frozen order", label)
    list(dimensions = dim(x), storage_mode = typeof(x), finite = TRUE,
         row_count = length(ids), first_row = ids[[1L]], last_row = ids[[length(ids)]],
         first_column = cols[[1L]], last_column = cols[[length(cols)]])
}

expected_summary <- function(x, sample_ids, probs, rle = FALSE) {
    out <- data.frame(sample_accession = sample_ids,
                      n_total = rep(nrow(x), ncol(x)),
                      n_missing = integer(ncol(x)),
                      n_nonfinite = integer(ncol(x)),
                      n_nonfinite_including_missing = integer(ncol(x)),
                      stringsAsFactors = FALSE, check.names = FALSE)
    qnames <- c("q0", "q01", "q25", "q50", "q75", "q99", "q100")
    for (j in seq_along(qnames)) out[[qnames[[j]]]] <- vapply(seq_len(ncol(x)), function(i) {
        as.numeric(stats::quantile(x[, i], probs = probs[[j]], type = 7,
                                   names = FALSE, na.rm = FALSE))
    }, numeric(1))
    if (rle) out$rle_iqr <- out$q75 - out$q25
    out
}

compare_summary <- function(path, label, expected, fields) {
    got <- read_stable_csv(path, label)
    tab <- got$value
    if (nrow(tab) != nrow(expected)) stopf("%s has %d rows; expected %d", label, nrow(tab), nrow(expected))
    if (!"sample_accession" %in% names(tab) ||
        !identical(as.character(tab$sample_accession), as.character(expected$sample_accession)))
        stopf("%s sample order differs from the frozen order", label)
    diffs <- list(); overall_max <- 0; failures <- 0L
    for (nm in fields) {
        if (!(nm %in% names(tab))) stopf("%s lacks column %s", label, nm)
        actual <- suppressWarnings(as.numeric(tab[[nm]])); wanted <- as.numeric(expected[[nm]])
        if (anyNA(actual)) stopf("%s column %s contains non-numeric values", label, nm)
        d <- abs(actual - wanted)
        m <- if (length(d)) max(d) else 0
        tol <- 1e-12 + 1e-12 * pmax(abs(actual), abs(wanted))
        bad <- !is.finite(d) | d > tol
        diffs[[nm]] <- list(max_abs_diff = m, exact_match_count = sum(actual == wanted),
                            row_count = length(actual), tolerant_mismatch_count = sum(bad))
        overall_max <- max(overall_max, m); failures <- failures + sum(bad)
    }
    if (failures) stopf("%s differs from independent recomputation in %d values", label, failures)
    list(artifact = got$artifact, fields = diffs, max_abs_diff = overall_max,
         tolerant_mismatch_count = failures)
}

main <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    script_arg <- commandArgs()[grepl("^--file=", commandArgs())]
    script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else "scripts/verify_gse84161_normalization.R"
    root <- normalizePath(file.path(dirname(normalizePath(script_path, mustWork = TRUE)), ".."), mustWork = TRUE)
    plan_arg <- "config/gse84161-normalization-plan.json"
    report_arg <- "reports/gse84161-normalization-independent-check.json"
    i <- 1L
    while (i <= length(argv)) {
        if (argv[[i]] == "--plan" && i < length(argv)) { i <- i + 1L; plan_arg <- argv[[i]]
        } else if (argv[[i]] == "--report" && i < length(argv)) { i <- i + 1L; report_arg <- argv[[i]]
        } else if (argv[[i]] %in% c("-h", "--help")) {
            cat("Usage: verify_gse84161_normalization.R [--plan PLAN.json] [--report REPORT.json]\n")
            return(invisible(NULL))
        } else stopf("unknown or incomplete argument: %s", argv[[i]])
        i <- i + 1L
    }
    plan_path <- resolve_path(plan_arg, root); report_path <- resolve_path(report_arg, root)
    plan_artifact <- file_record(plan_path); plan <- read_plan(plan_path)
    if (!identical(scalar_chr(field(plan, "status")), "frozen")) stopf("normalization plan is not frozen")
    input_hashes <- check_input_hashes(plan, root)
    method <- field(plan, "normalization", list())
    n_samples <- scalar_int(field(method, "expected_sample_count"), 30L)
    n_probes <- scalar_int(field(method, "expected_probe_rows"), 54675L)
    probs <- as.numeric(unlist(field(method, "quantile_probabilities", list()), use.names = FALSE))
    if (!identical(n_samples, 30L) || !identical(n_probes, 54675L) ||
        !identical(scalar_int(field(method, "quantile_type")), 7L) ||
        !identical(probs, c(0, .01, .25, .5, .75, .99, 1)))
        stopf("plan normalization dimensions or quantile specification is not the frozen method")
    samples <- load_samples(plan, root, n_samples)

    gpl_keys <- names(field(plan, "input_sha256"))[grepl("GSE84161-family-soft[.]txt$",
                                                            names(field(plan, "input_sha256")))]
    if (length(gpl_keys) != 1L) stopf("the plan does not pin one GSE84161 GPL SOFT file")
    gpl_path <- resolve_path(gpl_keys[[1L]], root)
    gpl_probe_ids <- read_gpl_probe_ids(gpl_path)

    map_cfg <- field(plan, "mapping", list())
    map_path <- resolve_path(field(map_cfg, "probe_map_path"), root)
    map_artifact <- file_record(map_path)
    map_expected_sha <- scalar_chr(field(map_cfg, "probe_map_sha256"))
    if (!is.null(map_expected_sha) && !identical(map_artifact$sha256, tolower(map_expected_sha)))
        stopf("probe map SHA-256 mismatch")
    map <- read_csv_chars(map_path)
    probe_col <- scalar_chr(field(map_cfg, "probe_id_column"), "probe_id")
    entrez_col <- scalar_chr(field(map_cfg, "entrez_id_column"), "entrez_gene_id")
    eligible_col <- scalar_chr(field(map_cfg, "eligible_column"), "eligible_probe")
    if (!all(c(probe_col, entrez_col, eligible_col) %in% names(map))) stopf("probe map lacks required columns")
    map_probe_ids <- as.character(map[[probe_col]])
    if (nrow(map) != n_probes || anyNA(map_probe_ids) || any(!nzchar(map_probe_ids)) || anyDuplicated(map_probe_ids))
        stopf("probe map does not contain %d unique probe IDs", n_probes)
    eligible_text <- tolower(trimws(as.character(map[[eligible_col]])))
    if (any(!eligible_text %in% c("true", "false"))) stopf("probe map eligible column is not boolean")
    eligible <- eligible_text == "true"
    map_entrez <- trimws(as.character(map[[entrez_col]]))
    if (any(eligible & (is.na(map_entrez) | !nzchar(map_entrez) | grepl("///", map_entrez, fixed = TRUE))))
        stopf("eligible probe map rows do not contain one Entrez ID")
    eligible_count <- sum(eligible)
    gene_ids <- sort(unique(map_entrez[eligible]), method = "radix")
    expected_eligible <- scalar_int(field(map_cfg, "expected_eligible_probe_count"), 41834L)
    expected_genes <- scalar_int(field(map_cfg, "expected_entrez_count"), 20486L)
    if (!identical(eligible_count, expected_eligible) || !identical(length(gene_ids), expected_genes))
        stopf("probe map accounting is %d eligible probes and %d Entrez IDs; expected %d and %d",
              eligible_count, length(gene_ids), expected_eligible, expected_genes)

    outputs <- field(plan, "outputs", list())
    out_dir <- resolve_path(field(outputs, "ignored_output_dir", "work/gse84161-normalization/output"), root)
    qc_dir <- resolve_path(field(outputs, "versioned_qc_dir", "data/derived/gse84161-array-qc"), root)
    probe_path <- file.path(out_dir, "normalized_probe_rma.rds")
    entrez_path <- file.path(out_dir, "normalized_entrez_median.rds")
    probe_read <- read_stable_rds(probe_path, "normalized probe RDS")
    entrez_read <- read_stable_rds(entrez_path, "normalized Entrez RDS")
    probe <- probe_read$value; entrez <- entrez_read$value
    probe_info <- check_matrix(probe, "normalized probe matrix", n_probes, samples$accession)
    # The complete GPL/CDF row order is validated by the exact probe-ID join below;
    # map row order is deliberately not used to index the expression matrix.
    if (!setequal(rownames(probe), map_probe_ids)) stopf("normalized probe IDs do not match the complete probe map")
    if (!setequal(rownames(probe), gpl_probe_ids)) stopf("normalized probe IDs do not match the pinned GPL probe IDs")
    probe_order_matches_map <- identical(as.character(rownames(probe)), as.character(map_probe_ids))
    probe_order_matches_gpl <- identical(as.character(rownames(probe)), as.character(gpl_probe_ids))
    probe_order_mismatch_count <- sum(as.character(rownames(probe)) != as.character(gpl_probe_ids))
    probe_info$probe_id_set_matches_probe_map <- TRUE
    probe_info$probe_id_set_matches_gpl <- TRUE
    probe_info$probe_row_order_matches_probe_map <- probe_order_matches_map
    probe_info$probe_row_order_matches_gpl <- probe_order_matches_gpl
    probe_info$probe_row_order_mismatch_count_vs_gpl <- probe_order_mismatch_count
    entrez_info <- check_matrix(entrez, "normalized Entrez matrix", expected_genes, samples$accession, gene_ids)

    eligible_idx <- match(map_probe_ids[eligible], rownames(probe))
    if (anyNA(eligible_idx) || anyDuplicated(eligible_idx)) stopf("explicit probe-ID join is incomplete or duplicated")
    eligible_values <- probe[eligible_idx, , drop = FALSE]
    groups <- split(seq_len(nrow(eligible_values)), map_entrez[eligible], drop = TRUE)
    group_names <- sort(names(groups), method = "radix")
    recomputed_entrez <- t(vapply(group_names, function(g) {
        apply(eligible_values[groups[[g]], , drop = FALSE], 2L, stats::median, na.rm = FALSE)
    }, numeric(ncol(probe))))
    rownames(recomputed_entrez) <- group_names
    colnames(recomputed_entrez) <- colnames(probe)
    expected_entrez <- entrez[match(group_names, rownames(entrez)), , drop = FALSE]
    delta <- recomputed_entrez - expected_entrez
    abs_delta <- abs(delta)
    scale_delta <- pmax(pmax(abs(recomputed_entrez), abs(expected_entrez)), 1e-300)
    tol_delta <- 1e-12 + 1e-12 * scale_delta
    tolerant_bad <- !is.finite(abs_delta) | abs_delta > tol_delta
    if (any(tolerant_bad)) stopf("independent Entrez median recomputation differs in %d cells", sum(tolerant_bad))
    median_check <- list(total_cells = length(delta), exact_match_count = sum(delta == 0),
                         exact_mismatch_count = sum(delta != 0), tolerant_mismatch_count = sum(tolerant_bad),
                         max_abs_diff = max(abs_delta), max_relative_diff = max(abs_delta / scale_delta),
                         absolute_tolerance = 1e-12, relative_tolerance = 1e-12,
                         explicit_probe_id_join = TRUE, aggregation = "median of all eligible single-Entrez probes per array")

    normalized_expected <- expected_summary(probe, samples$accession, probs, rle = FALSE)
    normalized_qc_path <- resolve_path(field(outputs, "normalized_probe_summary_csv", file.path(qc_dir, "normalized-probe-summary.csv")), root)
    normalized_fields <- c("n_total", "n_missing", "n_nonfinite", "n_nonfinite_including_missing",
                           "q0", "q01", "q25", "q50", "q75", "q99", "q100")
    normalized_qc <- compare_summary(normalized_qc_path, "normalized probe QC", normalized_expected, normalized_fields)

    probe_medians <- apply(probe, 1L, stats::median, na.rm = FALSE)
    rle <- sweep(probe, 1L, probe_medians, FUN = "-")
    if (any(!is.finite(rle))) stopf("independent RLE contains non-finite values")
    rle_expected <- expected_summary(rle, samples$accession, probs, rle = TRUE)
    rle_qc_path <- resolve_path(field(outputs, "rle_summary_csv", file.path(qc_dir, "rle-summary.csv")), root)
    rle_fields <- c("n_total", "n_missing", "n_nonfinite", "n_nonfinite_including_missing",
                    "q0", "q01", "q25", "q50", "q75", "q99", "q100", "rle_iqr")
    rle_qc <- compare_summary(rle_qc_path, "RLE QC", rle_expected, rle_fields)

    order_path <- resolve_path(field(outputs, "sample_order_csv", file.path(qc_dir, "sample-order.csv")), root)
    order_qc <- read_stable_csv(order_path, "sample-order QC")
    if (nrow(order_qc$value) != n_samples || !all(c("source_order", "sample_accession") %in% names(order_qc$value)))
        stopf("sample-order QC has the wrong shape")
    order_got <- suppressWarnings(as.integer(order_qc$value$source_order))
    if (anyNA(order_got) || !identical(order_got, samples$order) ||
        !identical(as.character(order_qc$value$sample_accession), samples$accession))
        stopf("sample-order QC differs from pinned metadata order")

    audit_path <- resolve_path(field(outputs, "audit_json", "reports/gse84161-normalization-audit.json"), root)
    audit_artifact <- file_record(audit_path); audit <- read_json(audit_path)
    if (!identical(scalar_chr(field(audit, "mode")), "execute") ||
        !identical(scalar_chr(field(audit, "expression_values_read")), "TRUE") ||
        !identical(scalar_chr(field(audit, "expression_normalization_run")), "TRUE") ||
        !identical(scalar_chr(field(audit, "program_scores_computed")), "FALSE") ||
        !identical(scalar_chr(field(audit, "treatment_or_donor_contrasts_computed")), "FALSE"))
        stopf("primary normalization audit does not record a completed expression run")
    mapping_audit <- field(audit, "mapping", list())
    if (!identical(scalar_chr(field(mapping_audit, "status")), "completed"))
        stopf("primary normalization audit does not record completed Entrez aggregation")

    result <- list(
        verifier_id = "gse84161-normalization-independent-check-v1",
        status = "passed",
        checked_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        plan = list(path = plan_path, bytes = plan_artifact$bytes, sha256 = plan_artifact$sha256,
                    plan_id = scalar_chr(field(plan, "plan_id")), status = scalar_chr(field(plan, "status")),
                    input_sha256 = input_hashes),
        sources = list(gpl = list(path = gpl_path, sha256 = input_hashes[[gpl_keys[[1L]]]]$sha256,
                                 probe_count = length(gpl_probe_ids)),
                       sample_metadata = list(path = samples$path, sha256 = samples$sha256,
                                               sample_count = n_samples, sample_order = samples$accession),
                       probe_map = list(path = map_artifact$path, bytes = map_artifact$bytes,
                                        sha256 = map_artifact$sha256, probe_count = nrow(map),
                                        eligible_probe_count = eligible_count, entrez_gene_count = length(gene_ids))),
        primary_audit = audit_artifact,
        matrices = list(
            normalized_probe_rma = c(probe_read$artifact, probe_info),
            normalized_entrez_median = c(entrez_read$artifact, entrez_info)
        ),
        checks = list(
            probe_id_join = list(full_probe_set = TRUE, explicit_join = TRUE,
                                 joined_eligible_probe_count = length(eligible_idx),
                                 probe_id_set_matches_probe_map = TRUE,
                                 probe_id_set_matches_gpl = TRUE,
                                 probe_row_order_matches_probe_map = probe_order_matches_map,
                                 probe_row_order_matches_gpl = probe_order_matches_gpl,
                                 probe_row_order_mismatch_count_vs_gpl = probe_order_mismatch_count),
            entrez_median_recomputation = median_check,
            normalized_probe_qc = normalized_qc,
            rle_qc = rle_qc,
            sample_order_qc = order_qc$artifact
        ),
        limits = list(raw_cel_values_read = FALSE, program_scores_computed = FALSE,
                      treatment_or_donor_contrasts_computed = FALSE,
                      values_used = "completed normalized RMA matrices only for structural, quantile, RLE and fixed-probe median checks")
    )
    dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(result, report_path, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 17)
    cat(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 17), "\n")
    cat(sprintf("report_sha256=%s\n", sha256_file(report_path)))
    invisible(result)
}

main()
