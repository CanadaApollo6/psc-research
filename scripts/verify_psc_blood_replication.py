#!/usr/bin/env python3
"""Independent PSC-blood verification. Never import the analysis scripts.

The default preparation mode reads source hashes, sample headers and annotations
only. Effect calculations require --run plus a matching execution-freeze hash.
Patient-level matrices and samples remain in ignored work/.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
KNOWN = ('UBASH3A','ETS2','PRKD2','PFKFB3','BCL2L11','IFIT1','G0S2')
PLATES = ('PS449','PS450','PS451','PS452','PS453')
ARRAY_GROUPS = {'primary sclerosing cholangitis':'PSC','control':'Control',
                'ulcerative colitis':'UC','primary biliary cholangitis':'PBC',
                "Crohn's disease":'CD'}
CONTRASTS = ('PSC_vs_Control','PSC_vs_UC','PSC_vs_PBC','PSC_vs_CD')
REPORT = ROOT/'reports/psc-blood-independent-verification.json'


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def op_text(path):
    return gzip.open(path,'rt',newline='') if str(path).endswith('.gz') else open(path,newline='')


def stamp(): return datetime.now(timezone.utc).isoformat()


class Checks:
    def __init__(self): self.rows=[]
    def exact(self,name,observed,expected):
        passed=observed==expected
        if isinstance(passed,np.ndarray): passed=bool(passed.all())
        def compact(value):
            if isinstance(value,(list,tuple,dict)) and len(value)>30:
                encoded=json.dumps(value,sort_keys=True,default=str,separators=(',',':')).encode()
                return {'items':len(value),'sha256':hashlib.sha256(encoded).hexdigest()}
            return value
        self.rows.append({'check':name,'passed':bool(passed),'observed':compact(observed),'expected':compact(expected)})
        return bool(passed)
    def truth(self,name,value,detail=None):
        self.rows.append({'check':name,'passed':bool(value),'detail':detail})
        return bool(value)
    def numeric(self,name,actual,expected,rtol=1e-9,atol=1e-10):
        a=np.asarray(actual,dtype=float); b=np.asarray(expected,dtype=float)
        shape=a.shape==b.shape
        passed=shape and bool(np.allclose(a,b,rtol=rtol,atol=atol,equal_nan=True))
        if shape:
            finite=np.isfinite(a)&np.isfinite(b)
            diff=np.abs(a[finite]-b[finite])
            maxerr=float(diff.max()) if diff.size else None
            na_agree=bool(np.array_equal(np.isnan(a),np.isnan(b)))
        else: maxerr=None; na_agree=False
        self.rows.append({'check':name,'passed':passed,'values':int(a.size),
                          'shape_matches':shape,'missingness_matches':na_agree,
                          'max_abs_error':maxerr,'rtol':rtol,'atol':atol})
        return passed
    def require(self,name,value,detail=None):
        if not self.truth(name,value,detail): raise ValueError(name)
    @property
    def passed(self): return all(r['passed'] for r in self.rows)


def parse_soft(path):
    records={}; current=None
    with op_text(path) as f:
        for line in f:
            if line.startswith('^SAMPLE = '):
                gsm=line.split(' = ',1)[1].strip()
                if gsm in records: raise ValueError('duplicate GSM metadata')
                current={'gsm':gsm,'characteristics':{}}
                records[gsm]=current
            elif current is not None and line.startswith('!Sample_') and ' = ' in line:
                key,value=line.rstrip('\r\n').split(' = ',1)
                if key=='!Sample_characteristics_ch1':
                    k,v=value.split(':',1)
                    if k.strip() in current['characteristics']: raise ValueError('duplicate sample characteristic')
                    current['characteristics'][k.strip()]=v.strip()
                else:
                    current.setdefault(key,[]).append(value)
    return records


def one(record,key):
    vals=record.get(key,[])
    if len(vals)!=1: raise ValueError(f'{key}: expected one metadata field')
    return vals[0]


def count_header(path):
    with op_text(path) as f: return next(csv.reader(f))


def series_header(path):
    meta=[]
    with op_text(path) as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):
                header=next(csv.reader([next(f)],delimiter='\t'))
                return header,meta
            if line.startswith('!Sample_'): meta.append(next(csv.reader([line],delimiter='\t')))
    raise ValueError('no series matrix table')


def series_matrix(path,nprobe,selected):
    all_ids=[]; selected=np.asarray(selected,dtype=int)
    adult=np.empty((nprobe,len(selected)),dtype=float)
    with op_text(path) as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'): break
        header=next(csv.reader([next(f)],delimiter='\t'))
        for row in csv.reader(f,delimiter='\t'):
            if row and row[0].startswith('!series_matrix_table_end'): break
            if not row: continue
            i=len(all_ids)
            if i>=nprobe: raise ValueError('unexpected extra probe')
            values=np.array(row[1:],dtype=float)
            if len(values)!=len(header)-1 or not np.all(np.isfinite(values)&(values>0)):
                raise ValueError('nonpositive/nonfinite/width source array failure')
            all_ids.append(row[0]); adult[i]=np.log2(values[selected])
    if len(all_ids)!=nprobe or len(set(all_ids))!=nprobe: raise ValueError('array probe count/uniqueness')
    return all_ids,adult


def parse_entrez(field):
    tokens=re.split(r'///|[;,|]',str(field))
    tokens=[s.strip() for s in tokens]
    missing={'','---','na','n/a','null','none','nan'}
    if any(t.lower() in missing for t in tokens): return None,'missing_or_mixed_missing'
    if not all(re.fullmatch(r'[1-9][0-9]*',t) for t in tokens): return None,'invalid_token'
    unique=set(tokens)
    if len(unique)!=1: return None,'multiple_entrez'
    return next(iter(unique)),'single_exact_entrez'


def read_platform(path,probe_col,entrez_col,symbol_col):
    rows=[]; header=None
    with op_text(path) as f:
        for line in f:
            if line.startswith(('!','#','^')): continue
            values=next(csv.reader([line],delimiter='\t'))
            if header is None:
                if probe_col in values and entrez_col in values:
                    header=values
                continue
            if not values or all(v=='' for v in values): continue
            if len(values)!=len(header): raise ValueError('platform row width')
            row=dict(zip(header,values)); entrez,reason=parse_entrez(row[entrez_col])
            rows.append({'probe':row[probe_col],'entrez':entrez,'mapping_reason':reason,
                         'entrez_original':row[entrez_col],'symbol':row[symbol_col]})
    if header is None or not rows: raise ValueError('platform annotation table not found')
    if len({r['probe'] for r in rows})!=len(rows): raise ValueError('duplicate platform probes')
    return rows


def crosswalk(path,checks):
    # Reparse ORIGINAL source columns, not the historical derived eligibility flags.
    ens_all=set(); e2n=defaultdict(set); n2e=defaultdict(set); symbols=defaultdict(set)
    count=0; invalid_ens=0; invalid_entrez=0
    with op_text(path) as f:
        for row in csv.DictReader(f):
            count+=1
            rawens=row['ensembl_gene_id_original']; rawn=row['entrez_gene_id_original']
            ensg=re.sub(r'\.\d+$','',rawens.strip())
            if not re.fullmatch('ENSG[0-9]{11}',ensg): invalid_ens+=1; continue
            ens_all.add(ensg)
            if row['hgnc_symbol_original'].strip(): symbols[row['hgnc_symbol_original'].strip()].add(ensg)
            if not rawn.strip(): continue
            if not re.fullmatch('[1-9][0-9]*',rawn.strip()): invalid_entrez+=1; continue
            entrez=rawn.strip(); e2n[ensg].add(entrez); n2e[entrez].add(ensg)
    eligible={e:next(iter(ids)) for e,ids in e2n.items() if len(ids)==1 and len(n2e[next(iter(ids))])==1}
    checks.exact('crosswalk_full_relationship_rows',count,91748)
    checks.exact('crosswalk_full_ensembl_universe',len(ens_all),86411)
    checks.exact('crosswalk_invalid_ensembl_ids',invalid_ens,0)
    checks.exact('crosswalk_invalid_nonblank_entrez_ids',invalid_entrez,0)
    return eligible,symbols,{'rows':count,'ensembl_ids':len(ens_all),'reciprocal_pairs':len(eligible)},ens_all


def corr(p,method='fdr_bh',restore_missing=False):
    p=np.asarray(p,dtype=float); missing=~np.isfinite(p)
    q=multipletests(np.where(missing,1.,p),method=method)[1]
    if restore_missing: q[missing]=np.nan
    return q


def numeric_column(frame,names):
    for name in names:
        if name in frame: return pd.to_numeric(frame[name],errors='coerce').to_numpy()
    raise ValueError('missing required result column: '+repr(names))


def bool_column(frame,name):
    values=frame[name].astype(str).str.lower()
    if not values.isin(['true','false','1','0']).all(): raise ValueError('nonboolean '+name)
    return values.isin(['true','1']).to_numpy()


def ordered_csv(path,key,ids):
    frame=pd.read_csv(path,dtype={key:str})
    if frame[key].duplicated().any(): raise ValueError('duplicate output key '+key)
    frame=frame.set_index(key)
    if set(frame.index)!=set(ids): raise ValueError('output gene family mismatch: '+str(path))
    return frame.loc[list(ids)]


def verify_sources_and_design(source_manifest,array_plan,rna_plan,checks):
    pins=list(source_manifest['source_files'])
    for plan in (array_plan,rna_plan):
        entries=plan.get('sources',{}) if 'sources' in plan else plan.get('rnaseq',{})
        for value in entries.values():
            if isinstance(value,dict) and 'path' in value and 'sha256' in value: pins.append(value)
            if isinstance(value,list): pins.extend(v for v in value if isinstance(v,dict) and 'path' in v and 'sha256' in v)
    unique={e['path']:e for e in pins}
    for path,e in unique.items():
        checks.require('source_exists:'+path,(ROOT/path).exists())
        checks.exact('source_sha256:'+path,sha(ROOT/path),e['sha256'])
        if 'bytes' in e: checks.exact('source_size:'+path,(ROOT/path).stat().st_size,e['bytes'])
    cross=source_manifest['reused_gene_crosswalk']
    checks.exact('complete_crosswalk_sha256',sha(ROOT/cross['path']),cross['sha256'])
    amap=ROOT/array_plan['sources']['sample_metadata']['path']
    rmap=ROOT/rna_plan['rnaseq']['metadata']['path']
    array_meta=parse_soft(amap); rna_meta=parse_soft(rmap)
    checks.exact('array_all_source_sample_records',len(array_meta),370)
    checks.exact('RNA_all_source_sample_records',len(rna_meta),1035)
    # Complete joins across all three deposited RNA count headers, not just cases.
    titles={one(r,'!Sample_title'):r for r in rna_meta.values()}
    checks.exact('RNA_unique_metadata_titles',len(titles),len(rna_meta))
    headers=[]; header_receipts=[]
    for pin in source_manifest['source_files']:
        if 'GSE177044_' in pin['path'] and pin['path'].endswith('.csv.gz'):
            hdr=count_header(ROOT/pin['path']); names=hdr[2:]
            checks.exact('RNA_header_identifier_columns:'+Path(pin['path']).name,hdr[:2],['Geneid','gene_name'])
            checks.exact('RNA_header_unique_samples:'+Path(pin['path']).name,len(set(names)),len(names))
            headers.extend(names)
            header_receipts.append({'path':pin['path'],'samples':len(names),
                                    'header_sha256':hashlib.sha256(('\t'.join(hdr)+'\n').encode()).hexdigest()})
    checks.exact('RNA_all_header_samples_unique',len(set(headers)),len(headers))
    checks.exact('RNA_all_headers_exact_metadata_join',sorted(headers),sorted(titles))
    primary_header=count_header(ROOT/rna_plan['rnaseq']['counts'][0]['path'])
    selected_names=primary_header[2:]
    selected=[]
    for name in selected_names:
        r=titles[name]; c=r['characteristics']
        plate=c.get('sequencing plate',c.get('sequencing plate number (used for batch correction)'))
        selected.append({'sample_id':name,'gsm':r['gsm'],'disease':c['disease'],
                         'age':float(c['age']),'sex':c['Sex'],'plate':plate})
    rn=pd.DataFrame(selected).set_index('sample_id')
    checks.exact('RNA_Norway_count',len(rn),297)
    checks.exact('RNA_Norway_plates',sorted(rn['plate'].unique()),list(PLATES))
    checks.exact('RNA_Norway_diagnosis_counts',dict(rn['disease'].value_counts()),{'PSCUC':177,'Control':77,'PSC':43})
    eligible_meta={one(r,'!Sample_title') for r in rna_meta.values()
                   if r['characteristics'].get('sequencing plate') in PLATES}
    checks.exact('RNA_plate_selection_matches_primary_header',sorted(eligible_meta),sorted(selected_names))
    age=rn['age'].to_numpy()-rn['age'].mean()
    X=np.column_stack([np.ones(len(rn)),rn['disease'].ne('Control').astype(float),age,
                       rn['sex'].eq('male').astype(float),
                       *[rn['plate'].eq(p).astype(float) for p in PLATES[1:]]]).astype(float)
    xcols=['Intercept','disease_PSC','age_centered','sex_male',*[f'plate_{p}' for p in PLATES[1:]]]
    checks.exact('RNA_design_rank',int(np.linalg.matrix_rank(X)),8)
    row_counts=Counter(map(tuple,X)); cooks_samples=np.array([row_counts[tuple(v)]>=3 for v in X])
    header,series_meta=series_header(ROOT/array_plan['sources']['series_matrix']['path'])
    gsm_order=header[1:]
    checks.exact('array_header_unique_GSMs',len(set(gsm_order)),370)
    checks.exact('array_complete_header_to_metadata_join',sorted(gsm_order),sorted(array_meta))
    agegroups=[]; groups=[]; adult_indices=[]
    for i,gsm in enumerate(gsm_order):
        c=array_meta[gsm]['characteristics']; agegroups.append(c['age group'])
        groups.append(ARRAY_GROUPS[c['condition']])
        if c['age group']=='adult': adult_indices.append(i)
    counts=Counter(groups[i] for i in adult_indices)
    checks.exact('array_adult_counts',dict(counts),{'PSC':45,'Control':47,'UC':45,'PBC':90,'CD':48})
    checks.exact('array_child_count',sum(v=='child' for v in agegroups),95)
    # Series-matrix characteristic arrays independently confirm the quick SOFT join.
    series_gsm=None; characteristics=[]
    for row in series_meta:
        key=row[0].split(' = ',1)[0]
        if key=='!Sample_geo_accession': series_gsm=row[1:]
        if key=='!Sample_characteristics_ch1': characteristics.append(row[1:])
    checks.exact('array_series_metadata_GSM_order',series_gsm,gsm_order)
    for characteristic,expected in [('age group',agegroups),('condition',[array_meta[g]['characteristics']['condition'] for g in gsm_order])]:
        candidates=[r for r in characteristics if all(v.startswith(characteristic+':') for v in r)]
        checks.exact('array_series_'+characteristic+'_row_count',len(candidates),1)
        if candidates:
            checks.exact('array_series_'+characteristic+'_exact_join',[v.split(':',1)[1].strip() for v in candidates[0]],expected)
    annotation=array_plan['annotation']
    platform=read_platform(ROOT/array_plan['sources']['platform_annotation']['path'],
                           annotation['probe_column'],annotation['entrez_column'],annotation['symbol_column'])
    checks.exact('array_platform_rows',len(platform),array_plan['expected_annotation_rows'])
    e2n,symbols,cross_summary,ens_all=crosswalk(ROOT/cross['path'],checks)
    return {'rna_meta':rn,'X':X,'xcols':xcols,'cooks_filter_samples':cooks_samples,
            'array_gsms':gsm_order,'adult_indices':adult_indices,
            'adult_gsms':[gsm_order[i] for i in adult_indices],
            'adult_groups':np.array([groups[i] for i in adult_indices]),
            'platform':platform,'e2n':e2n,'symbols':symbols,'crosswalk_gene_universe':ens_all,
            'summary':{'rna_headers':header_receipts,'rna_design_columns':xcols,
                       'rna_design_sha256':hashlib.sha256(X.astype('<f8').tobytes()).hexdigest(),
                       'rna_sample_ids_sha256':hashlib.sha256(('\n'.join(selected_names)+'\n').encode()).hexdigest(),
                       'array_GSM_header_sha256':hashlib.sha256(('\n'.join(gsm_order)+'\n').encode()).hexdigest(),
                       'array_adult_GSM_order_sha256':hashlib.sha256(('\n'.join([gsm_order[i] for i in adult_indices])+'\n').encode()).hexdigest(),
                       'cooks_filter_eligible_sample_count':int(cooks_samples.sum()),
                       'cooks_filter_threshold_F99':float(stats.f.ppf(.99,X.shape[1],len(rn)-X.shape[1])),
                       'crosswalk':cross_summary}}


def verify_array(plan,design,output_dir,checks,workdir):
    probe_ids,probe_log=series_matrix(ROOT/plan['sources']['series_matrix']['path'],
                                     plan['expected_probe_count'],design['adult_indices'])
    ann={r['probe']:r for r in design['platform']}
    group=defaultdict(list)
    for i,probe in enumerate(probe_ids):
        if probe in ann and ann[probe]['entrez'] is not None: group[ann[probe]['entrez']].append(i)
    genes=sorted(group,key=lambda s:(int(s),s))
    matrix=np.array([np.median(probe_log[group[g]],axis=0) for g in genes])
    del probe_log
    checks.truth('array_all_gene_medians_finite',np.isfinite(matrix).all())
    effects=pd.read_csv(output_dir/'array-gene-effects.csv.gz',dtype={'entrez_id':str,'entrez_gene_id':str})
    idcol='entrez_id' if 'entrez_id' in effects else 'entrez_gene_id'
    checks.exact('array_effect_contrast_set',sorted(effects['contrast'].unique()),sorted(CONTRASTS))
    checks.exact('array_effect_total_rows',len(effects),len(genes)*4)
    genecol='entrez_id'
    universe=pd.read_csv(output_dir/'array-gene-universe.csv.gz',dtype=str,keep_default_na=False)
    if genecol not in universe: genecol='entrez_gene_id'
    checks.exact('array_full_gene_universe',sorted(universe[genecol]),sorted(genes))
    results={}
    label=design['adult_groups']; psc=matrix[:,label=='PSC']; n1=psc.shape[1]
    for contrast in CONTRASTS:
        control=contrast.split('_vs_',1)[1]; ctl=matrix[:,label==control]; n0=ctl.shape[1]
        # scipy independently computes Welch's statistic, P and Satterthwaite df.
        test=stats.ttest_ind(psc,ctl,axis=1,equal_var=False)
        estimate=psc.mean(axis=1)-ctl.mean(axis=1)
        v1=psc.var(axis=1,ddof=1); v0=ctl.var(axis=1,ddof=1)
        se=np.sqrt(v1/n1+v0/n0)
        t=np.asarray(test.statistic).copy(); p=np.asarray(test.pvalue).copy(); df=np.asarray(test.df).copy()
        unavailable=(v1==0)&(v0==0)
        t[unavailable]=np.nan; p[unavailable]=np.nan; df[unavailable]=np.nan
        q=corr(p,restore_missing=True)
        observed=effects[effects['contrast']==contrast].copy()
        checks.require('array_unique_ids:'+contrast,not observed[idcol].duplicated().any())
        checks.require('array_family_exact:'+contrast,set(observed[idcol])==set(genes))
        observed=observed.set_index(idcol).loc[genes]
        for key,expected,aliases in [
            ('effect',estimate,['estimate','log2_intensity_difference']),
            ('Welch_SE',se,['welch_se','se']),('Welch_t',t,['welch_t','welch_statistic','statistic']),
            ('Welch_df',df,['welch_df']),('Welch_P',p,['welch_p','pvalue']),
            ('Welch_BH_full_family',q,['welch_bh_q','padj'])]:
            checks.numeric('array_'+key+':'+contrast,numeric_column(observed,aliases),expected)
        for side,sig in [('low',-1),('high',1)]:
            expected=estimate+sig*stats.t.ppf(.975,df)*se
            checks.numeric('array_Welch_CI_'+side+':'+contrast,numeric_column(observed,['welch_ci95_'+side]),expected)
        # Independent HC3 closed form for two-group OLS, then optional statsmodels subset below.
        hcse=np.sqrt(v1/(n1-1)+v0/(n0-1)); z=np.divide(estimate,hcse,out=np.full_like(estimate,np.nan),where=hcse>0)
        hp=2*stats.norm.sf(np.abs(z)); hq=corr(hp,restore_missing=True)
        for name,expected in [('hc3_se',hcse),('hc3_z',z),('hc3_p',hp),('hc3_bh_q',hq)]:
            checks.numeric('array_'+name+':'+contrast,numeric_column(observed,[name]),expected)
        results[contrast]=pd.DataFrame({'entrez_id':genes,'effect':estimate,'pvalue':p,'qvalue':q}).set_index('entrez_id')
    # Original platform/matrix feature sets, not only gene members, remain auditable.
    feature=pd.read_csv(output_dir/'array-feature-join.csv.gz',dtype=str,keep_default_na=False)
    probe_col=next((c for c in ['probe_id','probe','ID'] if c in feature),None)
    if probe_col:
        checks.exact('array_complete_annotation_matrix_union',sorted(feature[probe_col]),sorted(set(probe_ids)|set(ann)))
    summary={'all_genes':len(genes),'all_contrasts':4,'scipy_Welch_comparisons':len(genes)*4,
             'eligible_probes':sum(map(len,group.values())),
             'probe_order_sha256':hashlib.sha256(('\n'.join(probe_ids)+'\n').encode()).hexdigest(),
             'gene_order_sha256':hashlib.sha256(('\n'.join(genes)+'\n').encode()).hexdigest(),
             'adult_gene_matrix_sha256':hashlib.sha256(matrix.astype('<f8').tobytes()).hexdigest()}
    # No patient matrix is exported from the independent verifier.
    return results,summary


def verify_rna(plan,design,output_dir,fit_work,checks,workdir):
    p=plan['rnaseq']
    raw=pd.read_csv(ROOT/p['counts'][0]['path'],keep_default_na=False)
    ids=raw['Geneid'].astype(str).tolist(); names=raw['gene_name'].astype(str).tolist()
    checks.exact('RNA_full_source_gene_rows',len(ids),p['expected_feature_count'])
    checks.exact('RNA_unique_source_gene_ids',len(set(ids)),len(ids))
    checks.exact('RNA_count_column_order',list(raw.columns[2:]),design['rna_meta'].index.tolist())
    Y=raw.iloc[:,2:].to_numpy(dtype=np.float64)
    checks.require('RNA_source_finite_nonnegative_integer',bool(np.all(np.isfinite(Y)&(Y>=0)&(Y==np.floor(Y)))))
    libs=Y.sum(axis=0)
    checks.require('RNA_full_unfiltered_library_totals_positive',bool(np.all(libs>0)))
    eligible=(Y>=p['filter']['min_count']).sum(axis=1)>=p['filter']['min_samples']
    X=design['X']; n=len(libs)
    normalized_ids=[re.sub(r'\.\d+$','',g) for g in ids]
    checks.exact('RNA_version_removal_no_collisions',len(set(normalized_ids)),len(ids))
    observed=pd.read_csv(output_dir/'rnaseq-primary.csv',dtype={'source_gene_id':str,'source_gene_name':str,'gene_id':str,'ensembl_gene_id_versionless':str})
    checks.require('RNA_output_unique_source_ids',not observed['source_gene_id'].duplicated().any())
    checks.require('RNA_full_output_gene_family',set(observed['source_gene_id'])==set(ids))
    observed=observed.set_index('source_gene_id').loc[ids]
    checks.exact('RNA_reported_versionless_ids',observed['ensembl_gene_id_versionless'].tolist(),normalized_ids)
    checks.exact('RNA_count_eligibility_all_genes',bool_column(observed,'eligible').tolist(),eligible.tolist())
    checks.exact('RNA_source_gene_names_preserved',observed['source_gene_name'].fillna('').tolist(),names)
    # Exact source-to-normalization and model-design joins, no inferred sample renames.
    norm=pd.read_csv(fit_work/'sample-normalization.csv',dtype={'sample_id':str}).set_index('sample_id')
    checks.require('RNA_normalization_sample_set',set(norm.index)==set(design['rna_meta'].index))
    norm=norm.loc[design['rna_meta'].index]
    d=pd.read_csv(fit_work/'design-matrix.csv',dtype={'sample_id':str}).set_index('sample_id')
    checks.exact('RNA_reported_design_column_names',list(d.columns),design['xcols'])
    checks.require('RNA_reported_design_sample_set',set(d.index)==set(design['rna_meta'].index))
    checks.numeric('RNA_complete_design_matrix',d.loc[design['rna_meta'].index].to_numpy(),X,rtol=0,atol=1e-12)
    checks.numeric('RNA_unfiltered_full_library_sums',numeric_column(norm,['library_size_all_source_features']),libs,rtol=0,atol=0)
    checks.numeric('RNA_count_eligible_library_sums',numeric_column(norm,['library_size_all_fitted_features']),Y[eligible].sum(axis=0),rtol=0,atol=0)
    checks.exact('RNA_Cook_filter_sample_eligibility',bool_column(norm,'cooks_filter_eligible_sample').tolist(),design['cooks_filter_samples'].tolist())
    allpositive=eligible & np.all(Y>0,axis=1)
    checks.require('RNA_ratio_has_all_positive_eligible_genes',bool(allpositive.any()))
    # SciPy's geometric mean, followed by the median of log ratios, is independent
    # of the analysis package's normalization implementation.
    geometric=stats.gmean(Y[allpositive],axis=1)
    sizef=np.exp(np.median(np.log(Y[allpositive]/geometric[:,None]),axis=0))
    checks.numeric('RNA_ratio_size_factors_all_samples',numeric_column(norm,['size_factor']),sizef,rtol=1e-10,atol=1e-11)
    eobs=observed.iloc[np.flatnonzero(eligible)]
    checks.numeric('RNA_baseMean_all_eligible_genes',numeric_column(eobs,['baseMean']),(Y[eligible]/sizef).mean(axis=1),rtol=1e-9,atol=1e-8)
    # All-row consistency checks of supplied NB results: not a second NB fit.
    lfc=numeric_column(eobs,['log2FoldChange']); se=numeric_column(eobs,['lfcSE'])
    z=numeric_column(eobs,['stat']); source_p=numeric_column(eobs,['pvalue'])
    expected_z=np.divide(lfc,se,out=np.full_like(lfc,np.nan),where=np.isfinite(se)&(se>0))
    finite_stat=np.isfinite(expected_z)
    checks.numeric('RNA_all_NB_Wald_stat_from_effect_SE',z[finite_stat],expected_z[finite_stat],rtol=1e-8,atol=1e-9)
    finite_p=np.isfinite(source_p)
    checks.numeric('RNA_all_available_NB_P_from_Wald_z',source_p[finite_p],2*stats.norm.sf(np.abs(expected_z[finite_p])),rtol=1e-8,atol=1e-12)
    checks.numeric('RNA_primary_BH_full_eligible_family',numeric_column(eobs,['padj']),corr(source_p),rtol=1e-9,atol=1e-12)
    cooks=bool_column(eobs,'cooks_outlier')
    status=eobs['model_status'].fillna('').astype(str)
    checks.truth('RNA_Cook_filtered_P_remains_missing',bool(np.all(~finite_p[cooks])))
    checks.truth('RNA_Cook_filtered_status_preserved',bool(status[cooks].str.contains('cooks_filtered_pvalue_unavailable').all()))
    checks.truth('RNA_available_Wald_has_tested_status',bool(status[finite_p].str.startswith('tested').all()))
    checks.truth('RNA_unavailable_Wald_status_retained',bool(status[~finite_p].str.contains('pvalue_unavailable').all()))
    checks.truth('RNA_ineligible_P_not_imputed',bool(observed.loc[~eligible,'pvalue'].isna().all()))
    checks.truth('RNA_ineligible_status_retained',bool(observed.loc[~eligible,'model_status'].str.contains('excluded_prespecified_count_filter').all()))
    max_cook=numeric_column(eobs,['max_cooks_distance_filter_eligible_samples'])
    cutoff=stats.f.ppf(.99,X.shape[1],X.shape[0]-X.shape[1])
    checks.truth('RNA_Cook_mask_implies_filtered_sample_distance_above_cutoff',bool(np.all(max_cook[cooks]>cutoff)))
    checks.truth('RNA_Cook_max_subset_le_all_samples',bool(np.all(max_cook<=numeric_column(eobs,['max_cooks_distance_all_samples'])+1e-12)))
    # The full library's three-higher-count safeguard is not reconstructible from
    # maxima alone. No assertion is made that every distance > cutoff is filtered.
    eligible_indices=np.flatnonzero(eligible).tolist()
    deterministic=sorted(eligible_indices,key=lambda i:hashlib.sha256(('psc-blood-independent-v1:'+ids[i]).encode()).hexdigest())[:20]
    known=[i for i,name in enumerate(names) if name in KNOWN]
    subset=sorted(set(deterministic+known),key=lambda i:ids[i])
    model_subset=[i for i in subset if eligible[i]]
    source_y=np.log2(1+Y[model_subset]/libs*1e6)
    primary={}; subgroup={}; loo={}
    diagnosis=design['rna_meta']['disease'].to_numpy()
    Xsub=np.column_stack([np.ones(n),(diagnosis=='PSC').astype(float),(diagnosis=='PSCUC').astype(float),X[:,2:]])
    source_plates=design['rna_meta']['plate'].to_numpy()
    for row_index,y in zip(model_subset,source_y):
        fit=sm.OLS(y,X).fit(cov_type='HC3',use_t=False)
        ci=fit.conf_int(alpha=.05)[1]
        primary[ids[row_index]]=[fit.params[1],fit.bse[1],fit.tvalues[1],fit.pvalues[1],ci[0],ci[1]]
        subfit=sm.OLS(y,Xsub).fit(cov_type='HC3',use_t=False)
        subgroup[ids[row_index]]={}
        for j,label in [(1,'ols_psc_alone_'),(2,'ols_pscuc_')]:
            ci=subfit.conf_int(alpha=.05)[j]
            subgroup[ids[row_index]][label]=[subfit.params[j],subfit.bse[j],subfit.tvalues[j],subfit.pvalues[j],ci[0],ci[1]]
        loo[ids[row_index]]={}
        for plate in PLATES:
            keep=source_plates!=plate
            # Rebuild plate indicators with lexical first remaining plate as reference.
            rem=[v for v in PLATES if v!=plate]
            xx=np.column_stack([np.ones(keep.sum()),(diagnosis[keep]!='Control').astype(float),
                                X[keep,2],X[keep,3],*[(source_plates[keep]==v).astype(float) for v in rem[1:]]])
            loo[ids[row_index]][plate]=sm.OLS(y[keep],xx).fit().params[1]
    columns=['effect_log2_cpm','hc3_se','normal_z','pvalue','ci95_lower','ci95_upper']
    subset_ids=list(primary)
    for j,suffix in enumerate(columns):
        checks.numeric('RNA_statsmodels_HC3_subset:'+suffix,
                       numeric_column(observed.loc[subset_ids],['ols_'+suffix]),[primary[g][j] for g in subset_ids],rtol=1e-7,atol=1e-9)
        for prefix in ['ols_psc_alone_','ols_pscuc_']:
            checks.numeric('RNA_statsmodels_subgroup_subset:'+prefix+suffix,
                           numeric_column(observed.loc[subset_ids],[prefix+suffix]),[subgroup[g][prefix][j] for g in subset_ids],rtol=1e-7,atol=1e-9)
    # Verify every reported OLS correction independently, without claiming all fits were refitted.
    for prefix in ['ols_','ols_psc_alone_','ols_pscuc_']:
        pp=numeric_column(eobs,[prefix+'pvalue'])
        checks.numeric('RNA_'+prefix+'BH_full_eligible_family',numeric_column(eobs,[prefix+'padj_bh']),corr(pp),rtol=1e-9,atol=1e-12)
    lout=pd.read_csv(output_dir/'rnaseq-leave-one-plate-out.csv',dtype=str)
    key=next((k for k in ['source_gene_id','gene_id','ensembl_gene_id_versionless'] if k in lout),None)
    if key is None: raise ValueError('leave-one-plate output missing gene ID')
    lout=lout.set_index(key)
    # In the current input IDs are unversioned; retain an explicit crosswalk otherwise.
    loo_keys=[g if g in lout.index else re.sub(r'\.\d+$','',g) for g in subset_ids]
    for plate in PLATES:
        checks.numeric('RNA_statsmodels_leave_one_plate_subset:'+plate,
                       numeric_column(lout.loc[loo_keys],['effect_without_'+plate]),[loo[g][plate] for g in subset_ids],rtol=1e-7,atol=1e-9)
    # The selected verification subset is recorded before inspecting its model agreement.
    subset_rows=[{'source_gene_id':ids[i],'source_gene_name':names[i],'eligible':bool(eligible[i]),
                  'reason':';'.join((['deterministic_SHA256_first20'] if i in deterministic else [])+(['prespecified_benchmark'] if i in known else []))}
                 for i in subset]
    pd.DataFrame(subset_rows).to_csv(workdir/'rnaseq-independent-subset.csv',index=False)
    diagnostics={}
    nonconverged_any=np.zeros(len(eobs),dtype=bool)
    for flag in ['genewise_converged','MAP_converged','LFC_converged']:
        vals=eobs[flag].fillna('missing').astype(str).str.lower()
        false=vals.eq('false').to_numpy(); nonconverged_any |= false
        diagnostics[flag]={'source_value_counts':dict(vals.value_counts()),
                           'false_with_finite_P':int((false&finite_p).sum()),
                           'false_with_primary_BH_le_0_05':int((false&(numeric_column(eobs,['padj'])<=.05)).sum())}
    diagnostics['any_false_flag_with_finite_P']=int((nonconverged_any&finite_p).sum())
    diagnostics['any_false_flag_genes']=int(nonconverged_any.sum())
    diagnostics['interpretation']='Finite Wald arithmetic is not evidence of fit convergence. Source convergence flags are retained; no gene excluded post-effect by verifier.'
    for flag in ['genewise_dispersions','fitted_dispersions','dispersions']:
        values=numeric_column(eobs,[flag])
        diagnostics[flag]={'nonfinite':int((~np.isfinite(values)).sum()),'nonpositive_finite':int(((values<=0)&np.isfinite(values)).sum())}
    log_files=[ROOT/'work/research-path/rnaseq-run.log',fit_work/'pydeseq2-fit.log']
    diagnostics['runtime_logs']={str(p.relative_to(ROOT)):{'sha256':sha(p),'warning_lines':[line for line in p.read_text().splitlines() if any(word in line.lower() for word in ['warning','overflow','converg','fallback'])]} for p in log_files if p.exists()}
    full_nonconverged=np.zeros(len(ids),dtype=bool); full_nonconverged[eligible]=nonconverged_any
    failed_flags_by_id={g:[flag for flag in ['genewise_converged','MAP_converged','LFC_converged'] if str(observed.at[g,flag]).lower()=='false'] for g in observed.index[full_nonconverged]}
    diagnostics['flagged_source_genes']=[{'source_gene_id':g,'source_gene_name':observed.at[g,'source_gene_name'],'false_flags':flags} for g,flags in failed_flags_by_id.items()]
    result=pd.DataFrame({'source_gene_id':ids,'gene_id':normalized_ids,'source_gene_name':names,
                         'eligible':eligible,'effect':numeric_column(observed,['log2FoldChange']),
                         'pvalue':numeric_column(observed,['pvalue']),'model_status':observed['model_status'].to_numpy(),
                         'any_source_convergence_flag_false':full_nonconverged,
                         'false_convergence_flags':[';'.join(failed_flags_by_id.get(g,[])) for g in ids]})
    summary={'source_genes':len(ids),'eligible_genes':int(eligible.sum()),'ratio_all_positive_genes':int(allpositive.sum()),
             'deterministic_genes':20,'benchmark_source_rows':len(known),'model_refit_subset_genes':len(model_subset),
             'subset':subset_rows,'NB_refit_performed':False,'NB_available_P_checks':int(finite_p.sum()),
             'cooks_filtered_genes':int(cooks.sum()),'unavailable_Wald_genes':int((~finite_p).sum()),
             'cooks_filter_eligible_samples':int(design['cooks_filter_samples'].sum()),
             'convergence_and_warnings':diagnostics,
             'source_counts_replaced':False,
             'full_source_gene_id_sha256':hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()}
    return result,summary


def verify_conjunctions(rna,array,design,output_dir,checks):
    e2n=design['e2n']; array_genes=set(array['PSC_vs_Control'].index)
    eligible=rna['eligible'].to_numpy(dtype=bool)
    mapped=[e2n.get(g) for g in rna['gene_id']]
    inarray=np.array([n in array_genes if n is not None else False for n in mapped])
    shared=eligible & np.array([n is not None for n in mapped]) & inarray
    family=rna.loc[shared].copy(); family['entrez_gene_id']=[mapped[i] for i in np.flatnonzero(shared)]
    family=family.sort_values('gene_id').set_index('gene_id')
    observed=pd.read_csv(output_dir/'cross-cohort-gene-results.csv.gz',dtype={'ensembl_gene_id':str,'entrez_gene_id':str})
    checks.require('joint_unique_Ensembl_gene_rows',not observed['ensembl_gene_id'].duplicated().any())
    checks.exact('joint_complete_eligible_shared_family',sorted(observed['ensembl_gene_id']),family.index.tolist())
    observed=observed.set_index('ensembl_gene_id').loc[family.index]
    checks.exact('joint_reciprocal_Entrez_identity',observed['entrez_gene_id'].tolist(),family['entrez_gene_id'].tolist())
    checks.numeric('joint_RNA_P_preserved',numeric_column(observed,['rna_pvalue']),family['pvalue'].to_numpy())
    checks.numeric('joint_RNA_effect_preserved',numeric_column(observed,['rna_effect']),family['effect'].to_numpy())
    arr={}
    for contrast,short in zip(CONTRASTS,['control','uc','pbc','cd']):
        rr=array[contrast].loc[family['entrez_gene_id']]
        arr[short]=rr
        checks.numeric('joint_array_'+short+'_P',numeric_column(observed,['array_'+short+'_pvalue']),rr['pvalue'].to_numpy())
        checks.numeric('joint_array_'+short+'_effect',numeric_column(observed,['array_'+short+'_effect']),rr['effect'].to_numpy())
    outcomes={}
    published_summary=json.loads((output_dir/'summary.json').read_text())
    checks.exact('joint_summary_shared_family',published_summary['shared_gene_family'],len(family))
    for name,components in [('replication',['control']),('uc_support',['control','uc']),('all_disease_support',['control','uc','pbc','cd'])]:
        pp=np.column_stack([family['pvalue'],*[arr[c]['pvalue'].to_numpy() for c in components]])
        ff=np.column_stack([family['effect'],*[arr[c]['effect'].to_numpy() for c in components]])
        available=np.all(np.isfinite(pp)&(pp>=0)&(pp<=1)&np.isfinite(ff),axis=1)
        same=(np.all(ff>0,axis=1)|np.all(ff<0,axis=1))&available
        value=np.ones(len(family)); value[same]=np.max(pp[same],axis=1)
        statuses=np.full(len(family),'discordant_or_zero_direction',dtype=object)
        statuses[~available]='unavailable_constituent'; statuses[same]='concordant'
        checks.numeric('joint_'+name+'_maxP_full_family',numeric_column(observed,[name+'_pvalue']),value,rtol=1e-9,atol=1e-12)
        checks.exact('joint_'+name+'_missingness_sign_status',observed[name+'_status'].tolist(),statuses.tolist())
        for method,suffix in [('fdr_bh','bh'),('fdr_by','by')]:
            checks.numeric('joint_'+name+'_'+suffix+'_full_family',numeric_column(observed,[name+'_qvalue_'+suffix]),corr(value,method),rtol=1e-9,atol=1e-12)
        passing=corr(value)<=.05; by_passing=corr(value,'fdr_by')<=.05
        summary_counts={'q_bh_le_0_05':int(passing.sum()),'q_by_le_0_05':int(by_passing.sum()),
                        'passing_up_in_psc':int((passing&(family['effect'].to_numpy()>0)).sum()),
                        'passing_down_in_psc':int((passing&(family['effect'].to_numpy()<0)).sum()),
                        'status_counts':dict(Counter(statuses))}
        checks.exact('joint_summary_counts:'+name,published_summary['families'][name],summary_counts)
        flagged=passing&family['any_source_convergence_flag_false'].to_numpy()
        outcomes[name]={'family_size':len(family),**summary_counts,
                        'corrections_checked':['statsmodels fdr_bh','statsmodels fdr_by'],
                        'BH_passing_with_any_RNA_convergence_flag_false':int(flagged.sum()),
                        'BH_passing_flagged_gene_details':[{'ensembl_gene_id':family.index[i],
                            'source_gene_name':family.iloc[i]['source_gene_name'],
                            'false_convergence_flags':family.iloc[i]['false_convergence_flags'],
                            'independent_conjunction_BH_q':float(corr(value)[i])} for i in np.flatnonzero(flagged)]}
    coverage=pd.read_csv(output_dir/'cross-platform-coverage.csv.gz',dtype={'ensembl_gene_id':str,'entrez_gene_id':str},keep_default_na=False)
    # gene_id is the preserved source key; ensembl_gene_id is the annotation join
    # key and is intentionally blank for genes absent from the snapshot.
    checks.require('joint_full_coverage_unique_source_ids',not coverage['gene_id'].duplicated().any())
    checks.exact('joint_full_RNA_coverage_set',sorted(coverage['gene_id']),sorted(rna['gene_id']))
    coverage=coverage.set_index('gene_id').loc[rna['gene_id']]
    expected_annotation_ids=[g if g in design['crosswalk_gene_universe'] else '' for g in rna['gene_id']]
    checks.exact('joint_missing_annotation_IDs_preserved_separately',coverage['ensembl_gene_id'].tolist(),expected_annotation_ids)
    checks.exact('joint_eligibility_all_source_genes',bool_column(coverage,'eligible_shared_family').tolist(),shared.tolist())
    states=['rna_low_count' if not eligible[i] else ('no_unambiguous_cross_platform_identity' if mapped[i] is None else ('not_measured_on_array' if not inarray[i] else 'eligible_shared_family')) for i in range(len(rna))]
    checks.exact('joint_complete_family_status',coverage['family_status'].tolist(),states)
    return {'complete_shared_family':len(family),'complete_RNA_coverage_rows':len(rna),
            'gene_family_sha256':hashlib.sha256(('\n'.join(family.index)+'\n').encode()).hexdigest(),
            'coverage_status_counts':dict(Counter(states)),'families':outcomes,
            'NB_component_not_refit':True}


def json_default(value):
    if isinstance(value,np.generic): return value.item()
    if isinstance(value,np.ndarray): return value.tolist()
    if isinstance(value,Path): return str(value)
    raise TypeError(type(value).__name__)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run',action='store_true',help='Only after parent authorization and freeze.')
    ap.add_argument('--assays-only',action='store_true',help='Run source/assay checks while the integration is still being assembled.')
    ap.add_argument('--sources',default='config/psc-blood-sources.json')
    ap.add_argument('--array-plan',default='config/psc-blood-array-plan.json')
    ap.add_argument('--rna-plan',default='config/psc-blood-rnaseq-plan.json')
    ap.add_argument('--execution-freeze')
    ap.add_argument('--execution-freeze-sha256')
    ap.add_argument('--array-dir')
    ap.add_argument('--rna-dir')
    ap.add_argument('--rna-work-dir')
    ap.add_argument('--replication-dir')
    ap.add_argument('--work-dir',default='work/psc-blood-independent-verification')
    args=ap.parse_args()
    checks=Checks(); workdir=ROOT/args.work_dir; workdir.mkdir(parents=True,exist_ok=True)
    prior=json.loads(REPORT.read_text()) if REPORT.exists() else {}
    report={k:v for k,v in prior.items() if k in ['verification_id','prepared_at_utc','scope','method_independence','rna_subset_rule','array_scope','rna_scope','identity_scope','numerical_tolerances','source_manifest_at_preparation_sha256','protocol_at_preparation_sha256']}
    report.update({'last_checked_at_utc':stamp(),'status':'running' if args.run else 'source_header_preparation',
                   'effects_inspected_or_computed':False,'NB_refit_performed':False,'source_patient_matrices_versioned':False,
                   'implementation_sha256':sha(Path(__file__)),'python':sys.version,
                   'verification_implementation_repairs':['Closed a missing bracket before numerical execution.','Adapted draft expected_gene_count key to frozen expected_feature_count; unchanged63677 criterion.','Joined coverage by preserved source gene_id; annotation-side ensembl_gene_id is intentionally blank outside the snapshot. No analysis output or numerical tolerance changed.'],
                   'crosswalk_limit':'Reparsed all original relationship fields from the complete pinned Ensembl116 export; not fresh upstream identity validation.',
                   'array_annotation_snapshot':'GEO GPL10558.annot.gz Aug2016 reannotation, as frozen; not a manufacturer-original annotation claim.',
                   'input_plan_sha256':{p:sha(ROOT/p) for p in [args.sources,args.array_plan,args.rna_plan]},'checks':checks.rows})
    try:
        sources=json.loads((ROOT/args.sources).read_text())
        array_plan=json.loads((ROOT/args.array_plan).read_text())
        rna_plan=json.loads((ROOT/args.rna_plan).read_text())
        if args.run:
            if not args.execution_freeze or not args.execution_freeze_sha256:
                raise ValueError('Numerical run requires explicit execution-freeze path and hash')
            checks.require('execution_freeze_hash',sha(ROOT/args.execution_freeze)==args.execution_freeze_sha256)
            report['execution_freeze']={'path':args.execution_freeze,'sha256':args.execution_freeze_sha256}
            freeze=json.loads((ROOT/args.execution_freeze).read_text())
            for pinned_path,pinned_hash in freeze['sha256'].items():
                checks.exact('execution_pinned_sha256:'+pinned_path,sha(ROOT/pinned_path),pinned_hash)
            for label,plan in [('array',array_plan),('RNA',rna_plan)]:
                checks.require(label+'_assay_plan_frozen',bool(plan.get('frozen_at_utc')) and plan.get('status')!='draft')
            if not all([args.array_dir,args.rna_dir,args.rna_work_dir,args.replication_dir]):
                raise ValueError('Numerical run requires all aggregate/output and RNA work directories')
        design=verify_sources_and_design(sources,array_plan,rna_plan,checks)
        report['source_design_verification']=design['summary']
        if args.run:
            if not checks.passed: raise ValueError('Source/design checks failed before any effect verification')
            report['effects_inspected_or_computed']=True
            array,report['array_verification']=verify_array(array_plan,design,ROOT/args.array_dir,checks,workdir)
            rna,report['rna_verification']=verify_rna(rna_plan,design,ROOT/args.rna_dir,ROOT/args.rna_work_dir,checks,workdir)
            if not args.assays_only:
                report['conjunction_verification']=verify_conjunctions(rna,array,design,ROOT/args.replication_dir,checks)
            report['status']=('assay_checks_passed_waiting_for_integration' if args.assays_only else 'passed') if checks.passed else 'failed'
            report['verified_output_sha256']={str(p.relative_to(ROOT)):sha(p) for directory in [args.array_dir,args.rna_dir,args.replication_dir] for p in sorted((ROOT/directory).glob('*')) if p.is_file()}
        else:
            report['status']='prepared_waiting_for_parent_execution_signal' if checks.passed else 'source_preparation_failed'
            report['pending']=['Parent start signal/execution freeze','All fitted result directories','No real gene/group effects have been inspected']
    except Exception as exc:
        report['status']='error'; report['error']=f'{type(exc).__name__}: {exc}'
        report['checks']=checks.rows
        REPORT.write_text(json.dumps(report,indent=2,default=json_default,allow_nan=False)+'\n')
        raise
    report['checks']=checks.rows
    report['checks_total']=len(checks.rows); report['checks_failed']=sum(not r['passed'] for r in checks.rows)
    report['completed_at_utc']=stamp()
    REPORT.write_text(json.dumps(report,indent=2,default=json_default,allow_nan=False)+'\n')
    print(json.dumps({'status':report['status'],'checks_total':report['checks_total'],'checks_failed':report['checks_failed'],
                      'effects_inspected_or_computed':report['effects_inspected_or_computed'],'report':str(REPORT.relative_to(ROOT))},indent=2))
    if not checks.passed: raise SystemExit(1)


if __name__=='__main__': main()
