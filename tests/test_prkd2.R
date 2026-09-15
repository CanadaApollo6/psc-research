source('scripts/run_prkd2_signals.R')

# A real last-column variant must retain the same role under any SNP order.
ids <- paste0('v',seq_len(200))
a <- matrix(c(rep(0,199),35),nrow=1,dimnames=list(NULL,ids))
b <- a
one <- coloc::coloc.bf_bf(cbind(a,null=0),cbind(b,null=0),p12=1e-6,overlap.min=.9)
order <- c(200,seq_len(199))
two <- coloc::coloc.bf_bf(cbind(a[,order,drop=FALSE],null=0),cbind(b[,order,drop=FALSE],null=0),p12=1e-6,overlap.min=.9)
cols <- paste0('PP.H',0:4,'.abf')
stopifnot(one$summary$PP.H4.abf>.99,isTRUE(all.equal(as.numeric(one$summary[,..cols]),as.numeric(two$summary[,..cols]),tolerance=1e-12)))

# A missing dominant variant fails the full-vector gate despite many common rows.
missing <- suppressWarnings(coloc::coloc.bf_bf(cbind(a,null=0),cbind(b[,1:199,drop=FALSE],null=0),p12=1e-6,overlap.min=.9))
stopifnot(is.na(missing$summary$PP.H4.abf),missing$summary$nsnps==199)

# Exercise both sides of the 90% threshold independently of signal significance.
for(fraction in c(.89,.91)) {
  w <- c(rep(fraction/199,199),1-fraction)
  full <- matrix(log(w)+35,nrow=1,dimnames=list(NULL,ids))
  kept <- full[,1:199,drop=FALSE]
  result <- suppressWarnings(coloc::coloc.bf_bf(cbind(full,null=0),cbind(kept,null=0),p12=1e-6,overlap.min=.9))
  stopifnot(is.na(result$summary$PP.H4.abf)==(fraction<.9))
}
cat('Explicit null, variant-order invariance, missing dominant allele, and both coverage-threshold sides passed\n')
