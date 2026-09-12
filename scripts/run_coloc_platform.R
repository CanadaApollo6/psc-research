# Follow-up with audited per-platform sample counts; original results stay intact.
source('scripts/run_coloc.R')

read_platform <- function(gene,ld=FALSE) {
  read.csv(gzfile(file.path('data/derived/coloc-platform-inputs',paste0('GWAS-',gene,if(ld) '-ld' else '', '.csv.gz'))),stringsAsFactors=FALSE)
}

platform_logbf <- function(gwas) {
  coloc:::approx.bf.p(gwas$pvalue,gwas$MAF,'cc',gwas$N,gwas$case_fraction)$lABF
}

platform_pair <- function(gwas,qtl,original,p12,mode='unit') {
  if(anyDuplicated(gwas$snp)||anyDuplicated(qtl$snp))stop('Duplicate aligned SNPs')
  sd_y <- if(mode=='unit') 1 else coloc:::sdY.est(qtl$varbeta,qtl$MAF,unique(qtl$N))
  bg <- platform_logbf(gwas)
  bq <- coloc:::approx.bf.estimates(qtl$beta/sqrt(qtl$varbeta),qtl$varbeta,'quant',sdY=sd_y)$lABF
  common <- intersect(gwas$snp,qtl$snp)
  if(length(common)<100)stop('Fewer than 100 common SNPs')
  gi<-match(common,gwas$snp);qi<-match(common,qtl$snp)
  posterior <- coloc:::combine.abf(bg[gi],bq[qi],original$methods$priors$p1,original$methods$priors$p2,p12,quiet=TRUE)
  joint <- bg[gi]+bq[qi]
  list(summary=c(list(nsnps=length(common)),as.list(posterior),list(gwas_full_variants=nrow(gwas),qtl_full_variants=nrow(qtl),
           gwas_bf_mass_retained=bf_mass(bg,gwas$snp %in% common),qtl_bf_mass_retained=bf_mass(bq,qtl$snp %in% common),
           coverage_gate=bf_mass(bg,gwas$snp %in% common)>=0.9 && bf_mass(bq,qtl$snp %in% common)>=0.9,
           sdY_used=sd_y,gwas_min_p=min(gwas$pvalue),qtl_min_p=min(qtl$pvalue),gwas_top_snp=gwas$snp[which.max(bg)],qtl_top_snp=qtl$snp[which.max(bq)])),
       variants=data.frame(snp=common,position=gwas$position[gi],gwas_pvalue=gwas$pvalue[gi],qtl_pvalue=qtl$pvalue[qi],
                           gwas_logbf=bg[gi],qtl_logbf=bq[qi],conditional_H4_snp_probability=exp(joint-logsum(joint))))
}

platform_baseline <- function(original,amendment) {
  summaries<-list();details<-list();messages<-list()
  priors<-c(amendment$priors_unchanged$p12_primary,unlist(amendment$priors_unchanged$p12_sensitivity))
  for(region in original$regions) {
    gene<-region$gene_name;gwas<-read_platform(gene)
    for(dataset in original$datasets) {
      ds<-dataset$dataset_id;qtl<-read_trait(ds,gene)
      for(mode in c('unit','estimated'))for(prior in priors) {
        captured<-character()
        result<-withCallingHandlers(platform_pair(gwas,qtl,original,prior,mode),warning=function(w){captured<<-c(captured,conditionMessage(w));invokeRestart('muffleWarning')})
        summaries[[length(summaries)+1]]<-as.data.frame(c(list(gene_name=gene,dataset_id=ds,p12=prior,sd_mode=mode),result$summary),stringsAsFactors=FALSE)
        if(mode=='unit'&&prior==amendment$priors_unchanged$p12_primary)details[[length(details)+1]]<-cbind(gene_name=gene,dataset_id=ds,result$variants)
        if(length(captured))messages[[length(messages)+1]]<-list(gene_name=gene,dataset_id=ds,sd_mode=mode,p12=prior,warnings=unique(captured))
      }
      cat(jsonlite::toJSON(list(stage='platform-baseline',gene=gene,dataset=ds,status='complete'),auto_unbox=TRUE),'\n')
    }
  }
  write.csv(data.table::rbindlist(summaries),'reports/coloc-platform-baseline.csv',row.names=FALSE)
  con<-gzfile('reports/coloc-platform-variant-posteriors.csv.gz','wt')
  write.csv(data.table::rbindlist(details),con,row.names=FALSE);close(con)
  write_json(messages,'reports/coloc-platform-warnings.json')
}

platform_fit <- function(original,amendment,gene) {
  trait<-read_platform(gene,TRUE);m<-nrow(trait)
  prefix<-file.path('data/raw/coloc-platform-ld',paste0('GWAS-',gene))
  con<-file(paste0(prefix,'.f64'),'rb');R<-matrix(readBin(con,'double',n=m*m,size=8,endian='little'),m,m);close(con)
  if(any(!is.finite(R))||max(abs(R-t(R)))>1e-10||max(abs(diag(R)-1))>1e-10)stop('Invalid signed LD')
  dimnames(R)<-list(trait$snp,trait$snp)
  z<-qnorm(trait$pvalue/2,lower.tail=FALSE)*trait$canonical_sign
  cat(jsonlite::toJSON(list(stage='platform-LD-eigendecomposition',gene=gene,variants=m),auto_unbox=TRUE),'\n')
  e<-eigen(R,symmetric=TRUE)
  if(min(e$values)< -1e-8)stop('Non-PSD reference LD')
  attr(R,'eigen')<-e
  for(n in c(11386,3504,14890)) {
    warnings<-character();set.seed(original$methods$susie$seed)
    mismatch<-susieR::estimate_s_rss(z,R,n=n)
    if(n==11386) {
      diagnostic<-withCallingHandlers(susieR::kriging_rss(z,R,n=n,s=mismatch),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
      diagnostic$conditional_dist$snp<-trait$snp
      write.csv(diagnostic$conditional_dist,file.path('reports',paste0('coloc-platform-LD-diagnostics-',gene,'.csv')),row.names=FALSE)
    }
    cat(jsonlite::toJSON(list(stage='platform-SuSiE',gene=gene,n_approximation=n,mismatch=mismatch),auto_unbox=TRUE),'\n')
    fit<-tryCatch(withCallingHandlers(susieR::susie_rss(z=z,R=R,n=n,L=original$methods$susie$L,
             coverage=original$methods$susie$coverage,min_abs_corr=original$methods$susie$min_abs_corr,
             estimate_residual_variance=FALSE,max_iter=original$methods$susie$max_iter),
             warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(error=conditionMessage(e)))
    record<-list(gene_name=gene,status=if(!is.null(fit$error))'fit_error' else if(isTRUE(fit$converged))'converged' else 'not_converged',
                 scalar_n_approximation=n,actual_variant_N_range=range(trait$N),reference_donors=503,in_sample_LD=FALSE,
                 variants=m,min_eigenvalue=min(e$values),estimated_ld_mismatch_s=mismatch,warnings=unique(warnings),error=fit$error)
    if(is.null(fit$error)) {
      colnames(fit$lbf_variable)<-trait$snp;colnames(fit$alpha)<-trait$snp;fit$snp<-trait$snp;names(fit$pip)<-trait$snp
      record$n_credible_sets<-length(fit$sets$cs);record$iterations<-fit$niter
      record$credible_sets<-lapply(seq_along(fit$sets$cs),function(i){idx<-fit$sets$cs_index[i];list(component=idx,size=length(fit$sets$cs[[i]]),
              snps=trait$snp[fit$sets$cs[[i]]],purity=as.list(fit$sets$purity[i,,drop=FALSE]),log_bf=fit$lbf[idx])})
      saveRDS(fit,paste0(prefix,'-N',n,'.rds'),version=2,compress='gzip')
      write.csv(data.frame(snp=trait$snp,pip=fit$pip),file.path('reports',paste0('coloc-platform-pips-',gene,'-N',n,'.csv')),row.names=FALSE)
    }
    write_json(record,file.path('reports',paste0('coloc-platform-fit-',gene,'-N',n,'.json')))
    cat(jsonlite::toJSON(record[c('gene_name','scalar_n_approximation','status','n_credible_sets')],auto_unbox=TRUE),'\n')
  }
}

published_bfs <- function(ds,gene,indices) {
  rows<-read.delim(gzfile(file.path('data/derived/coloc-query-rows',paste0(ds,'-',gene,'-published-bfs.tsv.gz'))),stringsAsFactors=FALSE,check.names=FALSE)
  fields<-strsplit(rows$variant,'_',fixed=TRUE)
  keys<-vapply(fields,function(f){
    if(length(f)!=4)stop('Malformed published BF variant')
    if(nchar(f[3])==1&&nchar(f[4])==1&&f[3]%in%c('A','C','G','T')&&f[4]%in%c('A','C','G','T'))paste(c(f[1],f[2],sort(f[3:4])),collapse=':')
    else paste0('unmatched_nonSNV:',paste(f,collapse='_'))
  },character(1))
  if(anyDuplicated(keys))stop('Duplicate published BF identity')
  mat<-t(as.matrix(rows[,paste0('lbf_variable',indices),drop=FALSE]))
  if(any(is.na(mat)|mat==Inf))stop('Invalid published signal BF')
  colnames(mat)<-keys;rownames(mat)<-paste0('published_component_',indices)
  list(bf=mat,rows=rows)
}

compare_published <- function(original,amendment) {
  summaries<-list();mass_rows<-list();details<-list();coverage_rows<-list()
  priors<-c(amendment$priors_unchanged$p12_primary,unlist(amendment$priors_unchanged$p12_sensitivity))
  for(region in original$regions)for(dataset in original$datasets) {
    gene<-region$gene_name;ds<-dataset$dataset_id
    spec<-Filter(function(x)x$dataset_id==ds&&x$gene_id==region$gene_id,amendment$published_qtl_signal_sources)
    if(!length(spec)) {
      summaries[[length(summaries)+1]]<-data.frame(gene_name=gene,dataset_id=ds,status='no_published_QTL_credible_set',reason='No published 95% SuSiE credible set for this gene/context; not a null RNA effect')
      next
    }
    qb<-published_bfs(ds,gene,unlist(spec[[1]]$cs_indices))
    for(n in c(11386,3504,14890)) {
      fit_record<-read_json(file.path('reports',paste0('coloc-platform-fit-',gene,'-N',n,'.json')),simplifyVector=TRUE)
      if(fit_record$status!='converged'||is.null(fit_record$n_credible_sets)||fit_record$n_credible_sets==0) {
        summaries[[length(summaries)+1]]<-data.frame(gene_name=gene,dataset_id=ds,n_approximation=n,status='no_usable_GWAS_signal');next
      }
      fit<-readRDS(file.path('data/raw/coloc-platform-ld',paste0('GWAS-',gene,'-N',n,'.rds')))
      idx<-fit$sets$cs_index;gb<-fit$lbf_variable[idx,,drop=FALSE]
      common<-intersect(colnames(gb),colnames(qb$bf))
      for(i in seq_len(nrow(gb)))for(j in seq_len(nrow(qb$bf)))mass_rows[[length(mass_rows)+1]]<-data.frame(gene_name=gene,dataset_id=ds,n_approximation=n,
           gwas_component=idx[i],qtl_component=unlist(spec[[1]]$cs_indices)[j],common_variants=length(common),
           gwas_signal_mass_retained=bf_mass(gb[i,],colnames(gb)%in%common),qtl_signal_mass_retained=bf_mass(qb$bf[j,],colnames(qb$bf)%in%common))
      for(prior in priors) {
        result<-coloc::coloc.bf_bf(gb,qb$bf,p1=original$methods$priors$p1,p2=original$methods$priors$p2,p12=prior,
                                    overlap.min=0.9,trim_by_posterior=TRUE)
        if(is.null(result$summary)||nrow(result$summary)==0) {
          summaries[[length(summaries)+1]]<-data.frame(gene_name=gene,dataset_id=ds,n_approximation=n,p12=prior,status='insufficient_signal_overlap')
        } else {
          tab<-as.data.frame(result$summary)
          tab$gwas_component<-idx[tab$idx1];tab$qtl_component<-unlist(spec[[1]]$cs_indices)[tab$idx2]
          summaries[[length(summaries)+1]]<-cbind(gene_name=gene,dataset_id=ds,n_approximation=n,p12=prior,status='reference_LD_and_sample_heterogeneity_sensitivity',tab)
          if(n==11386&&prior==amendment$priors_unchanged$p12_primary&&!is.null(result$results))details[[length(details)+1]]<-cbind(gene_name=gene,dataset_id=ds,as.data.frame(result$results))
        }
      }
    }
  }
  write.csv(data.table::rbindlist(summaries,fill=TRUE),'reports/coloc-published-signal-comparison.csv',row.names=FALSE,na='')
  write.csv(data.table::rbindlist(mass_rows,fill=TRUE),'reports/coloc-published-signal-overlap.csv',row.names=FALSE,na='')
  if(length(details)) {
    con<-gzfile('reports/coloc-published-signal-variant-posteriors.csv.gz','wt')
    write.csv(data.table::rbindlist(details,fill=TRUE),con,row.names=FALSE,na='');close(con)
  }
}

platform_main <- function() {
  if(as.character(packageVersion('coloc'))!='5.2.3'||as.character(packageVersion('susieR'))!='0.14.2')stop('Pinned implementation changed')
  args<-commandArgs(trailingOnly=TRUE)
  original<-read_json('config/coloc-plan.json',simplifyVector=FALSE)
  amendment<-read_json('config/coloc-platform-plan.json',simplifyVector=FALSE)
  if(args[1]=='baseline')platform_baseline(original,amendment)
  else if(args[1]=='fit')platform_fit(original,amendment,args[2])
  else if(args[1]=='compare')compare_published(original,amendment)
  else stop('Expected baseline, fit or compare')
}
if(sys.nframe()==0)platform_main()
