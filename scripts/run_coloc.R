# Official coloc/susieR implementations with explicit full-region coverage audits.
suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(susieR))
suppressPackageStartupMessages(library(jsonlite))

read_trait <- function(id, gene, ld=FALSE) {
  read.csv(gzfile(file.path('data/derived/coloc-inputs', paste0(id,'-',gene,if(ld) '-ld' else '', '.csv.gz'))),
           stringsAsFactors=FALSE, check.names=FALSE)
}
write_json <- function(x, path) {
  jsonlite::write_json(x,path,pretty=TRUE,auto_unbox=TRUE,digits=NA,null='null',na='null')
}
logsum <- function(x) {
  m <- max(x)
  m + log(sum(exp(x-m)))
}
bf_mass <- function(lbf, keep) {
  if (!any(keep)) return(0)
  exp(logsum(lbf[keep])-logsum(lbf))
}

baseline_pair <- function(gwas, qtl, plan, p12, sd_mode='unit') {
  if(anyDuplicated(gwas$snp) || anyDuplicated(qtl$snp)) stop('Duplicate aligned variants')
  if(sd_mode=='unit') sd_y <- 1 else {
    if(length(unique(qtl$N))!=1) stop('Variable QTL sample sizes prevent the fixed sdY sensitivity')
    sd_y <- coloc:::sdY.est(qtl$varbeta,qtl$MAF,unique(qtl$N))
  }
  dg <- list(snp=gwas$snp,position=gwas$position,pvalues=gwas$pvalue,MAF=gwas$MAF,
             N=plan$gwas$n,s=plan$gwas$cases/plan$gwas$n,type='cc')
  dq <- list(snp=qtl$snp,position=qtl$position,beta=qtl$beta,varbeta=qtl$varbeta,
             MAF=qtl$MAF,sdY=sd_y,type='quant')
  common <- intersect(gwas$snp,qtl$snp)
  bg <- coloc:::approx.bf.p(gwas$pvalue,gwas$MAF,'cc',plan$gwas$n,plan$gwas$cases/plan$gwas$n)$lABF
  bq <- coloc:::approx.bf.estimates(qtl$beta/sqrt(qtl$varbeta),qtl$varbeta,'quant',sdY=sd_y)$lABF
  output <- NULL
  invisible(capture.output(output <- coloc::coloc.abf(dg,dq,p1=plan$methods$priors$p1,p2=plan$methods$priors$p2,p12=p12)))
  summary <- as.list(output$summary)
  summary$gwas_full_variants <- nrow(gwas)
  summary$qtl_full_variants <- nrow(qtl)
  summary$gwas_bf_mass_retained <- bf_mass(bg,gwas$snp %in% common)
  summary$qtl_bf_mass_retained <- bf_mass(bq,qtl$snp %in% common)
  summary$coverage_gate <- length(common)>=100 && summary$gwas_bf_mass_retained>=0.9 && summary$qtl_bf_mass_retained>=0.9
  summary$sdY_used <- sd_y
  summary$p12 <- p12
  summary$sd_mode <- sd_mode
  summary$gwas_min_p <- min(gwas$pvalue)
  summary$qtl_min_p <- min(qtl$pvalue)
  summary$gwas_top_snp <- gwas$snp[which.max(bg)]
  summary$qtl_top_snp <- qtl$snp[which.max(bq)]
  list(summary=summary,results=output$results)
}

run_baseline <- function(plan) {
  summaries <- list(); variants <- list(); warnings <- list()
  priors <- c(plan$methods$priors$p12_primary,unlist(plan$methods$priors$p12_sensitivity))
  for(region in plan$regions) {
    gene <- region$gene_name
    gwas <- read_trait('GWAS',gene)
    for(dataset in plan$datasets) {
      ds <- dataset$dataset_id; qtl <- read_trait(ds,gene)
      for(mode in c('unit','estimated')) for(prior in priors) {
        messages <- character()
        result <- tryCatch(withCallingHandlers(baseline_pair(gwas,qtl,plan,prior,mode),
                            warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart('muffleWarning')}),
                           error=function(e)list(error=conditionMessage(e)))
        if(!is.null(result$error)) {
          summaries[[length(summaries)+1]] <- data.frame(gene_name=gene,dataset_id=ds,sd_mode=mode,p12=prior,status='unavailable',reason=result$error)
        } else {
          summaries[[length(summaries)+1]] <- as.data.frame(c(list(gene_name=gene,dataset_id=ds,status='complete',reason=''),result$summary),stringsAsFactors=FALSE)
          if(mode=='unit' && prior==plan$methods$priors$p12_primary) {
            variants[[length(variants)+1]] <- cbind(gene_name=gene,dataset_id=ds,as.data.frame(result$results))
          }
        }
        if(length(messages)) warnings[[length(warnings)+1]] <- list(gene_name=gene,dataset_id=ds,sd_mode=mode,p12=prior,messages=unique(messages))
      }
      cat(jsonlite::toJSON(list(stage='baseline',gene=gene,dataset=ds,status='complete'),auto_unbox=TRUE),'\n')
    }
  }
  all_summary <- data.table::rbindlist(summaries,fill=TRUE)
  write.csv(all_summary,'reports/coloc-baseline.csv',row.names=FALSE,na='')
  connection <- gzfile('reports/coloc-baseline-variant-posteriors.csv.gz','wt')
  write.csv(data.table::rbindlist(variants,fill=TRUE),connection,row.names=FALSE,na='');close(connection)
  write_json(warnings,'reports/coloc-baseline-warnings.json')
}

fit_reference <- function(plan,ds,gene) {
  trait <- read_trait(ds,gene,TRUE)
  output <- file.path('data/raw/coloc-ld',paste0(ds,'-',gene))
  dimension <- nrow(trait)
  con <- file(paste0(output,'.f64'),'rb')
  R <- matrix(readBin(con,'double',n=dimension*dimension,size=8,endian='little'),dimension,dimension)
  close(con)
  if(any(!is.finite(R)) || max(abs(R-t(R)))>1e-10 || max(abs(diag(R)-1))>1e-10) stop('Invalid signed LD')
  dimnames(R) <- list(trait$snp,trait$snp)
  if(ds=='GWAS') {
    z <- qnorm(trait$pvalue/2,lower.tail=FALSE)*trait$canonical_sign
    n <- plan$gwas$n
  } else {
    if(length(unique(trait$N))!=1) stop('Variable sample sizes are not compatible with this fixed RSS fit')
    z <- trait$beta/sqrt(trait$varbeta);n <- unique(trait$N)
  }
  messages <- character()
  cat(jsonlite::toJSON(list(stage='LD-eigendecomposition',gene=gene,dataset=ds,variants=dimension),auto_unbox=TRUE),'\n')
  eigen_R <- eigen(R,symmetric=TRUE)
  if(min(eigen_R$values)< -1e-8) stop('Reference LD is not positive semidefinite')
  attr(R,'eigen') <- eigen_R
  s <- susieR::estimate_s_rss(z,R,n=n)
  cat(jsonlite::toJSON(list(stage='LD-diagnostics',gene=gene,dataset=ds,estimated_mismatch=s),auto_unbox=TRUE),'\n')
  diagnostics <- withCallingHandlers(susieR::kriging_rss(z,R,n=n,s=s),
                    warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart('muffleWarning')})
  diagnostics$conditional_dist$snp <- trait$snp
  write.csv(diagnostics$conditional_dist,file.path('reports',paste0('coloc-LD-diagnostics-',ds,'-',gene,'.csv')),row.names=FALSE)
  set.seed(plan$methods$susie$seed)
  cat(jsonlite::toJSON(list(stage='SuSiE-fit',gene=gene,dataset=ds),auto_unbox=TRUE),'\n')
  fit <- tryCatch(withCallingHandlers(susieR::susie_rss(z=z,R=R,n=n,L=plan$methods$susie$L,
                        coverage=plan$methods$susie$coverage,min_abs_corr=plan$methods$susie$min_abs_corr,
                        estimate_residual_variance=FALSE,max_iter=plan$methods$susie$max_iter),
                    warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart('muffleWarning')}),
                    error=function(e)list(error=conditionMessage(e)))
  status <- if(!is.null(fit$error)) 'fit_error' else if(!isTRUE(fit$converged)) 'not_converged' else 'converged'
  record <- list(dataset_id=ds,gene_name=gene,status=status,n_variants=dimension,n=n,reference_donors=503,
                 in_sample_LD=FALSE,estimated_ld_mismatch_s=s,min_eigenvalue=min(eigen_R$values),
                 reference_rank=sum(eigen_R$values>1e-8),warnings=unique(messages),error=fit$error,
                 diagnostic_flagged_variants=sum(diagnostics$conditional_dist$logLR>2 & abs(diagnostics$conditional_dist$z)>2))
  if(is.null(fit$error)) {
    colnames(fit$lbf_variable) <- trait$snp
    colnames(fit$alpha) <- trait$snp
    fit$snp <- trait$snp
    names(fit$pip) <- trait$snp
    record$n_credible_sets <- length(fit$sets$cs)
    record$iterations <- fit$niter
    record$credible_sets <- lapply(seq_along(fit$sets$cs),function(i){
      idx <- fit$sets$cs_index[i]
      list(component=idx,snps=trait$snp[fit$sets$cs[[i]]],size=length(fit$sets$cs[[i]]),
           purity=as.list(fit$sets$purity[i,,drop=FALSE]),log_bf=fit$lbf[idx])
    })
    saveRDS(fit,paste0(output,'.rds'),version=2,compress='gzip')
    write.csv(data.frame(snp=trait$snp,pip=fit$pip),file.path('reports',paste0('coloc-reference-pips-',ds,'-',gene,'.csv')),row.names=FALSE)
  }
  write_json(record,file.path('reports',paste0('coloc-reference-fit-',ds,'-',gene,'.json')))
  cat(jsonlite::toJSON(record[c('dataset_id','gene_name','status','n_credible_sets','estimated_ld_mismatch_s')],auto_unbox=TRUE),'\n')
}

compare_reference <- function(plan) {
  summaries <- list();all_details <- list()
  priors <- c(plan$methods$priors$p12_primary,unlist(plan$methods$priors$p12_sensitivity))
  for(region in plan$regions) {
    gene <- region$gene_name
    records <- lapply(c('GWAS','QTD000031'),function(ds)read_json(file.path('reports',paste0('coloc-reference-fit-',ds,'-',gene,'.json')),simplifyVector=TRUE))
    if(any(vapply(records,function(x)x$status!='converged' || is.null(x$n_credible_sets) || x$n_credible_sets==0,logical(1)))) {
      summaries[[length(summaries)+1]] <- data.frame(gene_name=gene,dataset_id='QTD000031',status='unavailable',reason='Fit did not converge or no credible set in one trait')
      next
    }
    fits <- lapply(c('GWAS','QTD000031'),function(ds)readRDS(file.path('data/raw/coloc-ld',paste0(ds,'-',gene,'.rds'))))
    for(prior in priors) {
      result <- coloc::coloc.susie(fits[[1]],fits[[2]],p1=plan$methods$priors$p1,p2=plan$methods$priors$p2,
                                    p12=prior,overlap.min=0.9,trim_by_posterior=TRUE)
      if(is.null(result$summary) || nrow(result$summary)==0) {
        summaries[[length(summaries)+1]] <- data.frame(gene_name=gene,dataset_id='QTD000031',p12=prior,status='unavailable',reason='No signal pair meets 90% posterior-mass overlap')
      } else {
        tab <- cbind(gene_name=gene,dataset_id='QTD000031',p12=prior,status='external_LD_sensitivity',as.data.frame(result$summary))
        summaries[[length(summaries)+1]] <- tab
        if(prior==plan$methods$priors$p12_primary && !is.null(result$results)) all_details[[length(all_details)+1]] <- cbind(gene_name=gene,as.data.frame(result$results))
      }
    }
  }
  write.csv(data.table::rbindlist(summaries,fill=TRUE),'reports/coloc-reference-comparison.csv',row.names=FALSE,na='')
  if(length(all_details)) {
    con<-gzfile('reports/coloc-reference-variant-posteriors.csv.gz','wt')
    write.csv(data.table::rbindlist(all_details,fill=TRUE),con,row.names=FALSE,na='');close(con)
  }
}

main <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if(as.character(packageVersion('coloc'))!='5.2.3' || as.character(packageVersion('susieR'))!='0.14.2') stop('Pinned package version changed')
  plan <- read_json('config/coloc-plan.json',simplifyVector=FALSE)
  if(args[1]=='baseline') run_baseline(plan)
  else if(args[1]=='fit') fit_reference(plan,args[2],args[3])
  else if(args[1]=='compare') compare_reference(plan)
  else stop('Expected baseline, fit or compare')
}
if(sys.nframe()==0) main()
