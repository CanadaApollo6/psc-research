"""Audit and summarize the frozen PSC/CD4 colocalisation results without refitting."""

import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[1]/'work/matplotlib-coloc'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import norm

from coloc_common import ROOT, save_json, sha


LABELS = {
    'QTD000031': 'BLUEPRINT naive CD4', 'QTD000439': 'DICE Tfh memory',
    'QTD000444': 'DICE Th17 memory', 'QTD000449': 'DICE Th1 memory',
    'QTD000454': 'DICE Th2 memory', 'QTD000459': 'DICE Th1-17 memory',
    'QTD000464': 'DICE Treg 464*', 'QTD000469': 'DICE Treg 469*',
    'QTD000479': 'DICE resting CD4', 'QTD000484': 'DICE activated CD4',
}


def read(name):
    return json.loads((ROOT/name).read_text())


def check_posterior(table):
    columns = [f'PP.H{i}.abf' for i in range(5)]
    measured = table[columns].dropna()
    if not ((measured >= 0) & (measured <= 1)).all().all():
        raise ValueError('Invalid hypothesis probabilities')
    if not np.allclose(measured.sum(axis=1), 1, atol=1e-10, rtol=0):
        raise ValueError('Hypothesis probabilities do not sum to one')
    return len(measured)


def curate_signals(raw, overlap):
    """Account for every signal pair, including rows suppressed by coloc's gate."""
    rows=[]
    for row in raw[raw.status.eq('no_published_QTL_credible_set')].to_dict('records'):
        rows.append({k:row[k] for k in ['gene_name','dataset_id','status','reason']})
    for record in overlap.to_dict('records'):
        gate=(record['common_variants']>=100 and record['gwas_signal_mass_retained']>=.9 and record['qtl_signal_mass_retained']>=.9)
        for prior in [1e-6,1e-7,1e-5]:
            entry={**record,'p12':prior,'coverage_gate':gate}
            candidate=raw
            for key in ['gene_name','dataset_id','n_approximation','gwas_component','qtl_component']:
                candidate=candidate[candidate[key].eq(record[key])]
            candidate=candidate[candidate.p12.eq(prior)]
            if not gate:
                entry.update(status='insufficient_signal_overlap',reason='Below fixed 90% signal-weight coverage; no hypothesis posterior is interpretable')
            else:
                if len(candidate)!=1 or check_posterior(candidate)!=1:
                    raise ValueError('Passing signal coverage lacks one valid posterior vector')
                entry.update(status='distinct_signals_under_reference_LD_sensitivity' if candidate['PP.H3.abf'].iloc[0]>.9 else 'reference_LD_sensitivity_result')
                entry.update({f'PP.H{i}.abf':float(candidate[f'PP.H{i}.abf'].iloc[0]) for i in range(5)})
            rows.append(entry)
    return pd.DataFrame(rows)


def gwas_logbf(table):
    """Coloc's case-control P-value ABF, cross-checked against recorded R output."""
    f, n, s = table.MAF.to_numpy(), table.N.to_numpy(), table.case_fraction.to_numpy()
    variance = 1 / (2*n*f*(1-f)*s*(1-s))
    shrinkage = .2**2 / (.2**2 + variance)
    z = norm.isf(table.pvalue.to_numpy()/2)
    return .5*(np.log1p(-shrinkage) + shrinkage*z*z)


def published_audit(plan):
    cs = pd.read_csv(ROOT/'data/derived/coloc-published-credible-sets.csv.gz')
    results = []
    for spec in plan['published_qtl_signal_sources']:
        ds = spec['dataset_id']
        bf = pd.read_csv(ROOT/f'data/derived/coloc-query-rows/{ds}-PFKFB3-published-bfs.tsv.gz', sep='\t')
        if bf.variant.duplicated().any():
            raise ValueError('Duplicate full published BF variant')
        values = bf[[f'lbf_variable{i}' for i in range(1,11)]].to_numpy()
        weights = np.exp(values-logsumexp(values, axis=0))
        pip = 1-np.prod(1-weights, axis=1)
        # The published CSV repeats variants for rsID aliases. Audit equality
        # first, then count distinct alleles; retain the original file unchanged.
        source = cs[cs.dataset_id.eq(ds)]
        identity = ['dataset_id', 'cs_id', 'variant']
        for _, aliases in source.groupby(identity):
            if len(aliases.drop(columns='rsid').drop_duplicates()) != 1:
                raise ValueError('Conflicting credible-set alias rows')
        source = source.drop_duplicates(identity)
        joined = source.merge(bf[['variant']].assign(recomputed_pip=pip), on='variant', validate='many_to_one')
        if len(joined) != len(source) or not np.allclose(joined.pip, joined.recomputed_pip, atol=1e-10, rtol=0):
            raise ValueError('Complete BF vectors do not reproduce published credible-set PIPs')
        for component in spec['cs_indices']:
            subset = joined[joined.cs_id.eq(spec['gene_id']+f'_L{component}')]
            if subset.cs_size.nunique() != 1 or int(subset.cs_size.iloc[0]) != len(subset):
                raise ValueError('Distinct credible-set size differs from published size')
            alleles = bf.variant.str.split('_')
            is_snv = alleles.apply(lambda a: len(a)==4 and len(a[2])==len(a[3])==1 and a[2] in 'ACGT' and a[3] in 'ACGT')
            top = np.argsort(-weights[:,component-1])[:10]
            results.append({
                'dataset_id': ds, 'qtl_component': component, 'complete_gene_variants': len(bf),
                'published_distinct_cs_variants': len(subset),
                'published_source_cs_rows_including_aliases': len(cs[cs.dataset_id.eq(ds)&cs.cs_id.eq(spec['gene_id']+f'_L{component}')]),
                'maximum_pip_reconstruction_difference': float((subset.pip-subset.recomputed_pip).abs().max()),
                'non_snv_signal_bf_mass': float(weights[~is_snv,component-1].sum()),
                'top_signal_variants': [{'source_variant': bf.variant.iloc[i], 'component_bf_weight': float(weights[i,component-1])} for i in top],
            })
    return results


def save_figure(fig, name):
    folder = ROOT/'reports/figures'
    folder.mkdir(exist_ok=True)
    fig.savefig(folder/(name+'.png'), dpi=170, facecolor='white', metadata={'Software':'PSC public-data research'})
    svg=folder/(name+'.svg')
    fig.savefig(svg, facecolor='white', metadata={'Date':None, 'Creator':'PSC public-data research'})
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)


def regional_figure(gene, contexts, anchor, region):
    gwas = pd.read_csv(ROOT/f'data/derived/coloc-platform-inputs/GWAS-{gene}.csv.gz')
    fig, axes = plt.subplots(1+len(contexts), 1, figsize=(11,.8+2.15*(1+len(contexts))), sharex=True)
    fig.subplots_adjust(left=.1, right=.97, top=.90, bottom=.16, hspace=.45)
    fig.suptitle(f'{gene}: PSC and RNA association patterns', x=.1, ha='left', fontsize=17, weight='bold')
    for i, ax in enumerate(axes):
        table = gwas if i==0 else pd.read_csv(ROOT/f'data/derived/coloc-inputs/{contexts[i-1]}-{gene}.csv.gz')
        if i==0:
            ax.scatter(table.position/1e6, -np.log10(table.pvalue), s=8, color='#273f66', alpha=.55, linewidths=0)
            title = 'PSC discovery GWAS · all three genotyping platforms'
        else:
            shared = table.snp.isin(gwas.snp)
            for mask, color, label in [(~shared,'#c4c9d0','SNP absent from aligned GWAS'),(shared,'#167a7c','SNP shared with aligned GWAS')]:
                ax.scatter(table.loc[mask,'position']/1e6, -np.log10(table.loc[mask,'pvalue']), s=8, color=color, alpha=.65, linewidths=0, label=label)
            title = LABELS[contexts[i-1]]+' · '+contexts[i-1]
        ax.axvline(anchor/1e6, color='#b25b25', lw=1.1, ls='--')
        ax.set_title(title, loc='left', fontsize=11)
        ax.set_ylabel('−log₁₀ P')
        ax.set_ylim(bottom=0)
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='y', alpha=.15)
    axes[-1].legend(loc='upper right', frameon=False, fontsize=8, markerscale=2)
    axes[-1].set_xlim(region['start_grch38_0based']/1e6, region['end_grch38_exclusive']/1e6)
    axes[-1].set_xlabel(f"Chromosome {region['chromosome'][3:]} · GRCh38 position (Mb)")
    fig.text(.1,.025,'Dashed line: previously screened variant. Full fixed TSS ±1 Mb window; eligible SNPs only.\nDifferent peaks are descriptive; omitted indels and missing variants limit colocalisation. *Treg naive/memory label unresolved.',fontsize=9,color='#4c5561')
    save_figure(fig,'coloc-'+gene.lower()+'-regions')


def coverage_figure(primary):
    fig, axes = plt.subplots(1,2,figsize=(12,6.4),sharey=True)
    fig.subplots_adjust(left=.20,right=.98,top=.82,bottom=.14,wspace=.16)
    fig.suptitle('Does each comparison retain the association evidence?', x=.03,y=.97,ha='left',fontsize=17,weight='bold')
    fig.text(.03,.895,'Fraction of single-signal Bayes-factor weight retained on shared SNPs.\nThis is a coverage check within the eligible SNPs; it excludes indel evidence.',fontsize=11)
    for ax, gene in zip(axes,['PFKFB3','BCL2L11']):
        part=primary[primary.gene_name.eq(gene)].set_index('dataset_id').loc[list(LABELS)]
        y=np.arange(len(part))
        ax.barh(y-.17,part.gwas_bf_mass_retained,height=.3,label='PSC',color='#273f66')
        ax.barh(y+.17,part.qtl_bf_mass_retained,height=.3,label='RNA',color='#167a7c')
        ax.axvline(.9,ls='--',color='#b25b25',lw=1.2)
        ax.set_title(gene,loc='left',weight='bold'); ax.set_xlim(0,1.04)
        ax.set_yticks(y,list(LABELS.values()))
        ax.set_xticks([0,.25,.5,.75,1],['0%','25%','50%','75%','100%'])
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='x',alpha=.12)
    axes[0].set_ylim(len(LABELS)-.5,-.5)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper right',bbox_to_anchor=(.98,.91),ncol=2,frameon=False)
    fig.text(.03,.03,'Dashed line: fixed 90% threshold for both traits. *Treg naive/memory label unresolved.\nPublished signal-wise coverage is a separate check and retains omitted indels in its denominator.',fontsize=9,color='#4c5561')
    save_figure(fig,'coloc-coverage')


def main():
    matplotlib.rcParams.update({'font.family':'DejaVu Sans','svg.hashsalt':'psc-coloc-v1'})
    paths=set()
    recorded_pins={}
    for lock in ['config/coloc-plan.json','config/coloc-platform-plan.json','config/coloc-analysis-inputs.json','config/coloc-platform-analysis-inputs.json']:
        record=read(lock)
        for name,digest in record.get('file_sha256',record.get('frozen_file_sha256',{})).items():
            target=ROOT/name
            if target.exists() and sha(target.read_bytes())!=digest:
                raise ValueError('Frozen input changed: '+name)
            if not target.exists() and not name.startswith('data/raw/'):
                raise ValueError('Missing required aggregate input: '+name)
            recorded_pins[name]=digest
        paths.add(lock)
    plan=read('config/coloc-plan.json');amendment=read('config/coloc-platform-plan.json')
    old=pd.read_csv(ROOT/'reports/coloc-baseline.csv')
    baseline=pd.read_csv(ROOT/'reports/coloc-platform-baseline.csv')
    signal=pd.read_csv(ROOT/'reports/coloc-published-signal-comparison.csv')
    overlap=pd.read_csv(ROOT/'reports/coloc-published-signal-overlap.csv')
    old_count,base_count,signal_count=map(check_posterior,[old,baseline,signal])
    curated=curate_signals(signal,overlap)
    curated.to_csv(ROOT/'reports/coloc-signal-status.csv',index=False)
    primary=baseline[baseline.sd_mode.eq('unit')&baseline.p12.eq(1e-6)].copy()
    if len(primary)!=20 or primary.duplicated(['gene_name','dataset_id']).any():
        raise ValueError('Incomplete or duplicated primary comparison family')
    full_range=baseline.groupby(['gene_name','dataset_id'])['PP.H4.abf'].agg(H4_min_all_sensitivities='min',H4_max_all_sensitivities='max')
    primary=primary.merge(full_range,on=['gene_name','dataset_id'],validate='one_to_one')
    primary.insert(2,'context',primary.dataset_id.map(LABELS))
    primary.to_csv(ROOT/'reports/coloc-context-summary.csv',index=False)
    missing=[];max_logbf_error=0.
    recorded=pd.read_csv(ROOT/'reports/coloc-platform-variant-posteriors.csv.gz')
    for gene in ['PFKFB3','BCL2L11']:
        gwas=pd.read_csv(ROOT/f'data/derived/coloc-platform-inputs/GWAS-{gene}.csv.gz')
        logbf=gwas_logbf(gwas);gwas['single_signal_bf_weight']=np.exp(logbf-logsumexp(logbf))
        check=recorded[recorded.gene_name.eq(gene)].merge(gwas[['snp']].assign(recomputed_logbf=logbf),on='snp',validate='many_to_one')
        error=float((check.gwas_logbf-check.recomputed_logbf).abs().max())
        if error>1e-10: raise ValueError('Python coverage calculation differs from official R implementation')
        max_logbf_error=max(max_logbf_error,error)
        for ds in LABELS:
            qtl=pd.read_csv(ROOT/f'data/derived/coloc-inputs/{ds}-{gene}.csv.gz')
            absent=gwas[~gwas.snp.isin(qtl.snp)].nlargest(20,'single_signal_bf_weight').copy()
            absent.insert(0,'dataset_id',ds);absent.insert(0,'gene_name',gene)
            missing.extend(absent.to_dict('records'))
    pd.DataFrame(missing).to_csv(ROOT/'reports/coloc-missing-gwas-variants.csv',index=False)
    published=published_audit(amendment)
    for region,contexts,anchor_id in [(plan['regions'][0],['QTD000031','QTD000469'],'rs7923054'),(plan['regions'][1],['QTD000031','QTD000439','QTD000444','QTD000454'],'rs72837826')]:
        gene=region['gene_name'];table=pd.read_csv(ROOT/f'data/derived/coloc-platform-inputs/GWAS-{gene}.csv.gz')
        anchor=table.loc[table.source_id.eq(anchor_id),'position'].item()
        regional_figure(gene,contexts,anchor,region)
    coverage_figure(primary)
    for pattern in ['config/coloc-*-queries.json','reports/coloc-*run-provenance.json']:
        paths.update(str(p.relative_to(ROOT)) for p in ROOT.glob(pattern))
    for name in ['reports/coloc-baseline.csv','reports/coloc-platform-baseline.csv','reports/coloc-platform-variant-posteriors.csv.gz','reports/coloc-published-signal-comparison.csv','reports/coloc-published-signal-overlap.csv','data/derived/coloc-published-credible-sets.csv.gz','scripts/summarize_coloc.py','scripts/verify_coloc_artifacts.py','scripts/replay_coloc_baseline.py','requirements-coloc-python.txt']:
        paths.add(name)
    for query in read('config/coloc-published-bf-queries.json')['queries']:
        if not query['archive_eof_verified'] or sha((ROOT/query['file']).read_bytes())!=query['sha256']:
            raise ValueError('Incomplete or changed full published BF evidence')
        paths.add(query['file'])
    fits={f'{gene}-N{n}':read(f'reports/coloc-platform-fit-{gene}-N{n}.json') for gene in ['PFKFB3','BCL2L11'] for n in [11386,3504,14890]}
    overlap['coverage_gate']=overlap.gwas_signal_mass_retained.ge(.9)&overlap.qtl_signal_mass_retained.ge(.9)&overlap.common_variants.ge(100)
    summary={
        'analysis':'PSC/CD4 regional colocalisation, original pass plus explicit platform-count amendment',
        'date':'2026-09-12','planned_gene_context_pairs':20,'completed_original_baseline_rows':old_count,
        'completed_amended_baseline_rows':base_count,'signal_rows_with_complete_hypothesis_probabilities':signal_count,
        'primary_snv_coverage_passes':int(primary.coverage_gate.sum()),
        'published_signal_audits':published,'published_signal_overlap':overlap.to_dict('records'),
        'signal_status_counts':curated.status.value_counts().to_dict(),
        'raw_signal_output_note':'The frozen coloc output omits BLUEPRINT component 1 after its coverage filter and returns incomplete NA vectors for DICE Treg. The separately derived signal-status table enumerates all 27 signal/prior/N combinations plus 18 gene-contexts without published credible sets; failed coverage has no posterior.',
        'gwas_fits':{k:{field:v[field] for field in ['status','scalar_n_approximation','variants','estimated_ld_mismatch_s','n_credible_sets']} for k,v in fits.items()},
        'maximum_python_vs_R_logbf_difference':max_logbf_error,
        'missing_variant_weight_definition':'Normalized single-signal ABF weight among eligible GWAS SNPs; not a clinical or model-free causal probability. Table retains top 20 absent SNPs per gene/context.',
        'interpretation':{
            'PFKFB3':'Eligible-SNP baseline favors distinct PSC/RNA signals in both earlier positive contexts. The strongest published RNA components fail full signal coverage because indels carry substantial evidence. Any interpretable secondary-signal result is restricted to that component and to approximate external-LD GWAS fits.',
            'BCL2L11':'Unresolved: all contexts fail the joint baseline evidence-coverage gate and no published 95% QTL credible set was found in the selected contexts.',
            'clinical':'No validated mechanism, treatment target, PSC progression inference or clinical recommendation.'},
        'source_sha256':dict(sorted({**recorded_pins,**{name:sha((ROOT/name).read_bytes()) for name in paths}}.items())),
        'generated_sha256':{str(p.relative_to(ROOT)):sha(p.read_bytes()) for p in sorted([ROOT/'reports/coloc-context-summary.csv',ROOT/'reports/coloc-missing-gwas-variants.csv',ROOT/'reports/coloc-signal-status.csv',*ROOT.glob('reports/figures/coloc-*')])},
    }
    save_json(ROOT/'reports/coloc-summary.json',summary)
    print(json.dumps({k:summary[k] for k in ['planned_gene_context_pairs','completed_amended_baseline_rows','signal_rows_with_complete_hypothesis_probabilities','primary_snv_coverage_passes','maximum_python_vs_R_logbf_difference']}))


if __name__=='__main__':
    main()
