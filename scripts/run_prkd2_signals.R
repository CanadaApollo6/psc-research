# Explicit null-model columns for the pinned coloc overlap calculation.
# Selected before any PRKD2 component comparison. Real variants remain intact.
source('scripts/run_prkd2.R')

prkd2_compare <- function(plan) {
  cs <- read.csv(gzfile('data/derived/prkd2-published-credible-sets.csv.gz'),stringsAsFactors=FALSE)
  rows <- list();details <- list();notes <- list()
  priors <- c(plan$methods$priors$p12_primary,unlist(plan$methods$priors$p12_sensitivity))
  for(ds in vapply(plan$datasets,function(d)d$dataset_id,character(1))) {
    selected <- cs[cs$dataset_id==ds,,drop=FALSE]
    ids <- sort(unique(as.integer(sub('^ENSG00000105287_L','',selected$cs_id))))
    if(!length(ids)||any(is.na(ids)))stop('No usable published component indices')
    qb <- prkd2_published(ds,ids)
    for(n in c(11386,3504,14890)) {
      record <- read_json(paste0('reports/prkd2-GWAS-fit-N',n,'.json'),simplifyVector=TRUE)
      if(record$status!='converged'||is.null(record$n_credible_sets)||record$n_credible_sets==0) {
        for(j in ids)for(prior in priors)rows[[length(rows)+1]] <- data.frame(dataset_id=ds,n_approximation=n,
             gwas_component=NA_integer_,qtl_component=j,p12=prior,status='no_usable_GWAS_signal')
        next
      }
      fit <- readRDS(paste0('data/raw/prkd2-ld/GWAS-PRKD2-N',n,'.rds'))
      for(i in fit$sets$cs_index)for(j in seq_along(ids)) {
        gb <- fit$lbf_variable[i,,drop=FALSE]; mb <- qb$bf[j,,drop=FALSE]
        common <- intersect(colnames(gb),colnames(mb))
        gm <- bf_mass(gb[1,],colnames(gb)%in%common); qm <- bf_mass(mb[1,],colnames(mb)%in%common)
        gate <- length(common)>=100 && gm>=.9 && qm>=.9
        for(prior in priors) {
          captured <- character()
          result <- withCallingHandlers(coloc::coloc.bf_bf(cbind(gb,null=0),cbind(mb,null=0),p1=plan$methods$priors$p1,p2=plan$methods$priors$p2,
                 p12=prior,overlap.min=.9,trim_by_posterior=TRUE),
                 warning=function(w){captured<<-c(captured,conditionMessage(w));invokeRestart('muffleWarning')})
          if(is.null(result$summary)||nrow(result$summary)!=1)stop('Unexpected single-component result shape')
          tab <- as.data.frame(result$summary)
          rows[[length(rows)+1]] <- cbind(data.frame(dataset_id=ds,n_approximation=n,gwas_component=i,qtl_component=ids[j],p12=prior,
              status=if(gate)'external_LD_sample_scope_sensitivity' else 'insufficient_signal_overlap',
              full_gwas_component_variants=ncol(gb),full_qtl_component_variants=ncol(mb),
              gwas_signal_mass_retained=gm,qtl_signal_mass_retained=qm,coverage_gate=gate),tab)
          if(length(captured))notes[[length(notes)+1]] <- list(dataset_id=ds,n_approximation=n,gwas_component=i,qtl_component=ids[j],p12=prior,warnings=unique(captured))
          if(gate && n==11386 && prior==plan$methods$priors$p12_primary) {
            if(is.null(result$results))stop('Passing component lacks variant posteriors')
            details[[length(details)+1]] <- cbind(dataset_id=ds,gwas_component=i,qtl_component=ids[j],as.data.frame(result$results))
          }
        }
      }
    }
  }
  write.csv(data.table::rbindlist(rows,fill=TRUE),'reports/prkd2-signal-comparison.csv',row.names=FALSE,na='')
  con <- gzfile('reports/prkd2-signal-variant-posteriors.csv.gz','wt')
  if(length(details))write.csv(data.table::rbindlist(details,fill=TRUE),con,row.names=FALSE,na='') else writeLines('dataset_id,gwas_component,qtl_component,snp,SNP.PP.H4.abf',con)
  close(con)
  write_json(notes,'reports/prkd2-signal-warnings.json')
  cat(jsonlite::toJSON(list(stage='signal-comparison',rows=length(rows),status='complete'),auto_unbox=TRUE),'\n')
}

if(sys.nframe()==0)prkd2_main()
