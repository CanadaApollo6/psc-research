# Separately frozen complete-allele PRKD2 extension; identical statistical settings.
source('scripts/run_coloc_platform.R')

prkd2_read <- function(id, suffix='PRKD2') {
  read.csv(gzfile(file.path('data/derived/prkd2-expanded-inputs',paste0(id,'-',suffix,'.csv.gz'))),
           stringsAsFactors=FALSE,check.names=FALSE)
}

prkd2_baseline <- function(plan) {
  g <- prkd2_read('GWAS'); gf <- prkd2_read('GWAS','full-evidence')
  bg <- platform_logbf(g); bgf <- platform_logbf(gf)
  rows <- list(); detail <- list(); notes <- list()
  priors <- c(plan$methods$priors$p12_primary,unlist(plan$methods$priors$p12_sensitivity))
  for(ds in vapply(plan$datasets,function(d)d$dataset_id,character(1))) {
    q <- prkd2_read(ds); qf <- prkd2_read(ds,'full-evidence')
    if(anyDuplicated(g$snp)||anyDuplicated(q$snp)||anyDuplicated(qf$source_variant))stop('Duplicate trait identity')
    common <- intersect(g$snp,q$snp); gi <- match(common,g$snp); qi <- match(common,q$snp)
    qfi <- match(q$source_variant[qi],qf$source_variant)
    if(length(unique(qf$N))!=1)stop('Variable N prevents frozen sdY sensitivity')
    for(mode in c('unit','estimated')) {
      warnings <- character()
      sd_y <- if(mode=='unit') 1 else withCallingHandlers(coloc:::sdY.est(qf$varbeta,qf$MAF,unique(qf$N)),
                   warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
      bq <- coloc:::approx.bf.estimates(q$beta/sqrt(q$varbeta),q$varbeta,'quant',sdY=sd_y)$lABF
      bqf <- coloc:::approx.bf.estimates(qf$beta/sqrt(qf$varbeta),qf$varbeta,'quant',sdY=sd_y)$lABF
      gm <- bf_mass(bgf,gf$source_row %in% g$source_row[gi])
      qm <- bf_mass(bqf,seq_len(nrow(qf)) %in% qfi)
      gate <- length(common)>=100 && gm>=.9 && qm>=.9
      for(prior in priors) {
        posterior <- coloc:::combine.abf(bg[gi],bq[qi],plan$methods$priors$p1,plan$methods$priors$p2,prior,quiet=TRUE)
        rows[[length(rows)+1]] <- as.data.frame(c(list(dataset_id=ds,sd_mode=mode,p12=prior,
          status=if(gate)'single_causal_baseline' else 'coverage_failed',nsnps=length(common),
          full_gwas_rows=nrow(gf),eligible_gwas_snvs=nrow(g),full_qtl_variants=nrow(qf),eligible_qtl_snvs=nrow(q),
          full_gwas_bf_mass_retained=gm,full_qtl_bf_mass_retained=qm,
          eligible_gwas_bf_mass_retained=bf_mass(bg,seq_len(nrow(g)) %in% gi),
          eligible_qtl_bf_mass_retained=bf_mass(bq,seq_len(nrow(q)) %in% qi),
          coverage_gate=gate,sdY_used=sd_y,gwas_top_snp=g$snp[which.max(bg)],qtl_top_snp=q$snp[which.max(bq)],
          gwas_min_p=min(g$pvalue),qtl_min_p=min(q$pvalue)),as.list(posterior)),stringsAsFactors=FALSE)
        if(mode=='unit'&&prior==plan$methods$priors$p12_primary) {
          joint <- bg[gi]+bq[qi]
          detail[[length(detail)+1]] <- data.frame(dataset_id=ds,snp=common,position=g$position[gi],
            gwas_pvalue=g$pvalue[gi],qtl_pvalue=q$pvalue[qi],gwas_logbf=bg[gi],qtl_logbf=bq[qi],
            conditional_H4_snp_probability=exp(joint-logsum(joint)))
        }
      }
      if(length(warnings))notes[[length(notes)+1]] <- list(dataset_id=ds,sd_mode=mode,warnings=unique(warnings))
    }
    cat(jsonlite::toJSON(list(stage='baseline',dataset=ds,status='complete'),auto_unbox=TRUE),'\n')
  }
  write.csv(data.table::rbindlist(rows),'reports/prkd2-expanded-baseline.csv',row.names=FALSE)
  con <- gzfile('reports/prkd2-expanded-baseline-variant-posteriors.csv.gz','wt')
  write.csv(data.table::rbindlist(detail),con,row.names=FALSE);close(con)
  write_json(notes,'reports/prkd2-expanded-baseline-warnings.json')
}

prkd2_fit <- function(plan) {
  trait <- prkd2_read('GWAS','PRKD2-ld'); m <- nrow(trait)
  prefix <- 'data/raw/prkd2-expanded-ld/GWAS-PRKD2'
  con <- file(paste0(prefix,'.f64'),'rb'); R <- matrix(readBin(con,'double',n=m*m,size=8,endian='little'),m,m);close(con)
  if(any(!is.finite(R))||max(abs(R-t(R)))>1e-10||max(abs(diag(R)-1))>1e-10)stop('Invalid signed LD')
  dimnames(R) <- list(trait$snp,trait$snp)
  z <- qnorm(trait$pvalue/2,lower.tail=FALSE)*trait$canonical_sign
  cat(jsonlite::toJSON(list(stage='LD-eigendecomposition',variants=m),auto_unbox=TRUE),'\n')
  e <- eigen(R,symmetric=TRUE)
  if(min(e$values)< -1e-8)stop('Non-PSD reference LD')
  attr(R,'eigen') <- e
  for(n in c(11386,3504,14890)) {
    warnings <- character();set.seed(plan$methods$susie$seed)
    mismatch <- susieR::estimate_s_rss(z,R,n=n)
    diagnostic <- withCallingHandlers(susieR::kriging_rss(z,R,n=n,s=mismatch),
      warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    diagnostic$conditional_dist$snp <- trait$snp
    write.csv(diagnostic$conditional_dist,paste0('reports/prkd2-expanded-LD-diagnostics-N',n,'.csv'),row.names=FALSE)
    cat(jsonlite::toJSON(list(stage='SuSiE',n_approximation=n,mismatch=mismatch),auto_unbox=TRUE),'\n')
    fit <- tryCatch(withCallingHandlers(susieR::susie_rss(z=z,R=R,n=n,L=plan$methods$susie$L,
        coverage=plan$methods$susie$coverage,min_abs_corr=plan$methods$susie$min_abs_corr,
        estimate_residual_variance=FALSE,max_iter=plan$methods$susie$max_iter),
        warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
        error=function(e)list(error=conditionMessage(e)))
    record <- list(status=if(!is.null(fit$error))'fit_error' else if(isTRUE(fit$converged))'converged' else 'not_converged',
       scalar_n_approximation=n,actual_variant_N_range=range(trait$N),reference_donors=503,in_sample_LD=FALSE,
       variants=m,min_eigenvalue=min(e$values),reference_rank=sum(e$values>1e-8),estimated_ld_mismatch_s=mismatch,
       flagged_ld_diagnostic_variants=sum(diagnostic$conditional_dist$logLR>2 & abs(diagnostic$conditional_dist$z)>2),
       warnings=unique(warnings),error=fit$error)
    if(is.null(fit$error)) {
      colnames(fit$lbf_variable) <- trait$snp;colnames(fit$alpha) <- trait$snp;fit$snp <- trait$snp;names(fit$pip) <- trait$snp
      record$n_credible_sets <- length(fit$sets$cs);record$iterations <- fit$niter
      record$credible_sets <- lapply(seq_along(fit$sets$cs),function(i){idx<-fit$sets$cs_index[i];list(component=idx,
        size=length(fit$sets$cs[[i]]),snps=trait$snp[fit$sets$cs[[i]]],purity=as.list(fit$sets$purity[i,,drop=FALSE]),log_bf=fit$lbf[idx])})
      saveRDS(fit,paste0(prefix,'-N',n,'.rds'),version=2,compress='gzip')
      write.csv(data.frame(snp=trait$snp,pip=fit$pip),paste0('reports/prkd2-expanded-GWAS-pips-N',n,'.csv'),row.names=FALSE)
      con <- gzfile(paste0('data/derived/prkd2-expanded-GWAS-component-bfs-N',n,'.csv.gz'),'wt')
      write.csv(data.frame(snp=trait$snp,t(fit$lbf_variable),check.names=FALSE),con,row.names=FALSE);close(con)
    }
    write_json(record,paste0('reports/prkd2-expanded-GWAS-fit-N',n,'.json'))
    cat(jsonlite::toJSON(record[c('scalar_n_approximation','status','n_credible_sets')],auto_unbox=TRUE),'\n')
  }
}

prkd2_published <- function(ds, indices) {
  rows <- read.csv(gzfile(paste0('data/derived/prkd2-expanded-published-bfs/',ds,'.csv.gz')),
                   stringsAsFactors=FALSE,check.names=FALSE)
  if(!all(rows$molecular_trait_id=='ENSG00000105287')||anyDuplicated(rows$snp))stop('Invalid mapped full published vectors')
  mat <- t(as.matrix(rows[,paste0('lbf_variable',indices),drop=FALSE]))
  if(any(is.na(mat)|mat==Inf))stop('Invalid component BF')
  colnames(mat) <- rows$snp;rownames(mat) <- as.character(indices)
  list(bf=mat,rows=rows)
}

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
      record <- read_json(paste0('reports/prkd2-expanded-GWAS-fit-N',n,'.json'),simplifyVector=TRUE)
      if(record$status!='converged'||is.null(record$n_credible_sets)||record$n_credible_sets==0) {
        for(j in ids)for(prior in priors)rows[[length(rows)+1]] <- data.frame(dataset_id=ds,n_approximation=n,
             gwas_component=NA_integer_,qtl_component=j,p12=prior,status='no_usable_GWAS_signal')
        next
      }
      fit <- readRDS(paste0('data/raw/prkd2-expanded-ld/GWAS-PRKD2-N',n,'.rds'))
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
  write.csv(data.table::rbindlist(rows,fill=TRUE),'reports/prkd2-expanded-signal-comparison.csv',row.names=FALSE,na='')
  con <- gzfile('reports/prkd2-expanded-signal-variant-posteriors.csv.gz','wt')
  if(length(details))write.csv(data.table::rbindlist(details,fill=TRUE),con,row.names=FALSE,na='') else writeLines('dataset_id,gwas_component,qtl_component,snp,SNP.PP.H4.abf',con)
  close(con)
  write_json(notes,'reports/prkd2-expanded-signal-warnings.json')
  cat(jsonlite::toJSON(list(stage='signal-comparison',rows=length(rows),status='complete'),auto_unbox=TRUE),'\n')
}

prkd2_main <- function() {
  if(as.character(packageVersion('coloc'))!='5.2.3'||as.character(packageVersion('susieR'))!='0.14.2')stop('Pinned runtime changed')
  args <- commandArgs(trailingOnly=TRUE)
  plan <- read_json('config/prkd2-plan.json',simplifyVector=FALSE)
  if(args[1]=='baseline')prkd2_baseline(plan)
  else if(args[1]=='fit')prkd2_fit(plan)
  else if(args[1]=='compare')prkd2_compare(plan)
  else stop('Expected baseline, fit or compare')
}
if(sys.nframe()==0)prkd2_main()
