# Reproduce the GSE84161 array runtime

Run these commands from the repository root. The runtime is isolated at
`work/gse84161-normalization/runtime`; the existing coloc runtime is not
modified. The base Linux-64 package lock is
`requirements-gse84161-array-linux-64.lock` (17,270 bytes,
SHA-256 `f04873f928f31d5ebc0f634788ca78e9f2d823eee0e47d800500c3e78477064b`).
The effective normalization plan must also pin the supplemental
[`preprocessCore` override manifest](../config/gse84161-preprocesscore-thread-override.json)
(SHA-256 `bdedf9e39f649146157906e7278df12e18b3b31d446119d8fb55c236b8555259`).

## Validate inputs

The archive and small source inventories are checked before any CEL is opened:

```bash
sha256sum \
  requirements-gse84161-array-linux-64.lock \
  data/raw/gse84161/GSE84161_RAW.tar \
  data/derived/gse84161-raw-arrays.csv \
  data/derived/gse84161-samples.csv \
  data/derived/gse84161-probe-map.csv.gz
tar -tf data/raw/gse84161/GSE84161_RAW.tar | wc -l
mkdir -p work/gse84161-normalization/cel
tar -xf data/raw/gse84161/GSE84161_RAW.tar -C work/gse84161-normalization/cel
```

The pinned raw archive is 146,739,200 bytes with SHA-256
`6d700c2e2cdcee9ef4eabb654b4c61f5d32952aa6f07bb00d11900f12ec78fff`. The
runner rechecks the archive member set and every extracted file against the
sample and array inventories before reading probe values.

## Create the base environment

The first environment was bootstrapped with the exact package specifications
recorded in its Conda history:

```bash
work/coloc-runtime/micromamba create -y \
  -p work/gse84161-normalization/runtime \
  -c conda-forge -c bioconda --override-channels \
  r-base=4.4.3 \
  bioconductor-affy=1.84.0=r44h3df3fcb_1 \
  bioconductor-hgu133plus2cdf=2.18.0=r44hdfd78af_14 \
  r-jsonlite=2.0.0
```

The explicit lock was exported after that bootstrap. For a fresh Linux-64
prefix, the lock-based form is:

```bash
work/coloc-runtime/micromamba create -y \
  -p <new-prefix> --file requirements-gse84161-array-linux-64.lock
```

A local dry run resolved all 146 locked packages with zero downloads because
the package cache was already populated. It did not test a clean host. The
lock records exact artifact URLs and hashes, so a clean replay needs access to
those channels and will fail on a changed artifact rather than silently
accepting it. The actual bootstrap did not use the lock because it was created
after the environment; no package-resolution fallback was needed.

The `hgu133plus2cdf` Conda package has a post-link hook. The hook calls
`installBiocDataPackage.sh hgu133plus2cdf-2.18.0`, reads the local
`dataURLs.json`, downloads
`https://bioconductor.org/packages/3.22/data/annotation/src/contrib/hgu133plus2cdf_2.18.0.tar.gz`
or its listed mirrors, checks MD5
`284fef2f0b777d7b53451538ddd53de3`, and installs the source package into the
R library. In this runtime the hook completed; the installed
`hgu133plus2cdf.rda` is 4,350,398 bytes with SHA-256
`30c2e18b6df8151eb22d6d5ca2470e6de68a4317d5d32c2ab1c7ad02687bb0ce`.
The lock verifies the 11,659-byte Conda hook package (SHA-256
`c5224648389dcb6a760a6a225068d57a949b0e24f86b386652d5831e844a0dfc`), but it
does not contain the post-link source tarball. A fresh or offline install
therefore needs the source tarball in its available cache, or network access
for the hook; no manual post-link invocation was used here.

## Apply the same-version runtime repair

The original Conda `preprocessCore` artifact is version 1.68.0, 201,647 bytes,
SHA-256
`f4a708cbea277d6ea11f53b70299057e4635dd9a106c68e68c90cca32a4b4a72`. Its
threaded shared object failed with `pthread_create()` error 22 during
background correction, including on a synthetic matrix. The isolated runtime
was repaired from the official Bioconductor 3.20 source:

```bash
mkdir -p work/gse84161-normalization/source work/gse84161-normalization/build
curl -L \
  https://bioconductor.org/packages/3.20/bioc/src/contrib/preprocessCore_1.68.0.tar.gz \
  -o work/gse84161-normalization/source/preprocessCore_1.68.0.tar.gz
sha256sum work/gse84161-normalization/source/preprocessCore_1.68.0.tar.gz
tar -xzf work/gse84161-normalization/source/preprocessCore_1.68.0.tar.gz \
  -C work/gse84161-normalization/build
runtime_bin="$PWD/work/gse84161-normalization/runtime/bin"
PATH="$runtime_bin:$PATH" "$runtime_bin/R" CMD INSTALL \
  --library="$PWD/work/gse84161-normalization/runtime/lib/R/library" \
  --configure-args=--disable-threading \
  "$PWD/work/gse84161-normalization/build/preprocessCore"
```

The source tarball is 127,397 bytes with SHA-256
`50d9b6c2c1e4758dd19cc11470951d7e0d60f5ba315701835bc13dc9cbf3fb62`. The
runtime compiler directory must be on `PATH`; otherwise configure cannot find
`x86_64-conda-linux-gnu-cc`. The build log records
`configure: Disabling threading for preprocessCore`. The resulting package is
still version 1.68.0. The complete installed-library inventory and build-log
hashes are pinned in the override manifest. This source rebuild is the only
fallback required for the observed host-threading failure; it leaves the base
Conda lock and base runtime manifest unchanged.

Verify the replacement before using the study arrays:

```bash
work/gse84161-normalization/runtime/bin/Rscript --vanilla -e '
library(preprocessCore)
x <- matrix(seq_len(600), nrow=200L, ncol=3L)
bg <- rma.background.correct(x)
q <- normalize.quantiles(bg)
stopifnot(identical(dim(bg), c(200L, 3L)), identical(dim(q), c(200L, 3L)),
          all(is.finite(bg)), all(is.finite(q)),
          identical(as.character(packageVersion("preprocessCore")), "1.68.0"))
cat("toy_background_quantile_ok\\n")
'
```

This smoke test passed. It verifies the repaired background and quantile
functions only; it does not read or normalize GSE84161 arrays.

## Run the frozen plan

After the plan has been frozen with the source, runtime and override hashes,
run the provenance gate and then the single study execution:

```bash
work/gse84161-normalization/runtime/bin/Rscript \
  scripts/normalize_gse84161.R \
  --plan config/gse84161-normalization-plan.json --check-only
work/gse84161-normalization/runtime/bin/Rscript \
  scripts/normalize_gse84161.R \
  --plan config/gse84161-normalization-plan.json --execute
```

The runner uses all 30 arrays in saved GSM order, `affy::rma` with background
correction, quantile normalization, `bgversion=2`, and the full CDF. It writes
the full probe and Entrez matrices to ignored workspace storage and small raw
PM, normalized-signal and RLE outputs under
`data/derived/gse84161-array-qc/`, with the completed run audit at
`reports/gse84161-normalization-audit.json`. Full vector PDFs are kept ignored;
their compact PNG renders are versioned. Do not reuse a prior normalized matrix when
the plan's input or override hashes change.

## Reproduction limits

The lock is a Linux-64 explicit environment record; the package cache and
post-link CDF source were available locally for this run. A clean-machine
replay has not been independently completed. The Conda CDF package's post-link
source is versioned by its URL and MD5 in the local hook metadata, while the
installed `.rda` SHA above verifies the result actually used. The
`preprocessCore` replacement is a local same-version source build, so its
compiler flags and installed file inventory must be checked in addition to the
base package lock. The toy check and identity/QC gates establish executable,
finite numeric inputs; they do not establish donor correspondence, remove
technical confounding, or support treatment-effect interpretation.
