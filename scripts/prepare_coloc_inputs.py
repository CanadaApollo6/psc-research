"""Audit every regional row and align complete traits before intersecting them."""

import csv
import gzip
import io
import json
import math
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from coloc_common import ROOT, ChainMap, complement, gzip_bytes, reciprocal_base, save_json, sha, verify_plan
from summarize_cd4_rna import gene_id, measured_numbers


def canonical_key(chrom, position, allele_a, allele_b):
    if len(allele_a) != 1 or len(allele_b) != 1 or allele_a not in 'ACGT' or allele_b not in 'ACGT' or allele_a == allele_b:
        raise ValueError('Not a two-allele SNV')
    low, high = sorted((allele_a, allele_b))
    return f'{chrom}:{int(position)}:{low}:{high}', low, high


def save_frame(frame, path):
    plain = frame.to_csv(index=False, lineterminator='\n').encode()
    path.write_bytes(gzip_bytes(plain) if path.name.endswith('.gz') else plain)


def load_queries(stage):
    path = ROOT / ('config/coloc-'+stage+'-queries.json')
    run = json.loads(path.read_text())
    if run['plan_sha256'] != sha((ROOT/'config/coloc-plan.json').read_bytes()):
        raise ValueError('Regional queries use a different plan')
    if run['source_manifest_sha256'] != sha((ROOT/'config/coloc-sources.json').read_bytes()):
        raise ValueError('Regional source manifest changed')
    for q in run['queries']:
        data = (ROOT/q['file']).read_bytes()
        if sha(data) != q['sha256'] or sha(gzip.decompress(data)) != q['uncompressed_sha256']:
            raise ValueError('Regional source response changed')
    return {q['id']: q for q in run['queries']}


def eligible_gwas(row, region, forward, reverse, reference):
    result = {**row, 'gene_name': region['gene_name'], 'status': '', 'snp': '', 'position': '',
              'risk_allele_grch38': '', 'canonical_effect_allele': '', 'canonical_sign': '', 'reference_risk_frequency': '',
              'control_maf': '', 'source_frequency_difference': ''}
    a0, a1 = row['allele_0'], row['allele_1']
    if len(a0) != 1 or len(a1) != 1 or a0 not in 'ACGT' or a1 not in 'ACGT' or a0 == a1:
        result['status'] = 'not_biallelic_snv'; return result
    chrom, pos = 'chr'+row['chr'], int(row['pos'])
    mapped = reciprocal_base(forward, reverse, chrom, pos-1)
    if mapped is None:
        result['status'] = 'nonunique_or_nonreciprocal_mapping'; return result
    dest, position0, strand = mapped
    if not region['start_grch38_0based'] <= position0 < region['end_grch38_exclusive']:
        result['status'] = 'outside_fixed_grch38_region'; return result
    if strand == '-': a0, a1 = complement(a0), complement(a1)
    key, _, effect = canonical_key(dest, position0+1, a0, a1)
    result.update(snp=key, position=position0+1, risk_allele_grch38=a1, canonical_effect_allele=effect,
                  canonical_sign=1 if a1 == effect else -1)
    try:
        p, frequency, odds = float(row['p']), float(row['freq_1_controls']), float(row['or'])
        if not all(math.isfinite(x) for x in (p, frequency, odds)) or not 0 < p <= 1 or not 0 < frequency < 1 or odds <= 0:
            raise ValueError('Invalid source statistics')
    except (ValueError, TypeError):
        result['status'] = 'invalid_source_statistics'; return result
    result['control_maf'] = min(frequency, 1-frequency)
    # An effect is always coded by source allele_1, explicitly called risk in the header.
    if odds < 1:
        result['status'] = 'source_risk_label_or_conflict'; return result
    if row['platform'] != 'both':
        result['status'] = 'platform_sample_size_unavailable'; return result
    if result['control_maf'] < 0.01-1e-12:
        result['status'] = 'below_fixed_maf'; return result
    ref = reference.get(key)
    if ref:
        reference_risk = ref['effect_frequency'] if a1 == effect else 1-ref['effect_frequency']
        result.update(reference_risk_frequency=reference_risk, source_frequency_difference=abs(reference_risk-frequency))
    if {a0, a1} in ({'A','T'}, {'C','G'}):
        if not ref:
            result['status'] = 'palindromic_reference_frequency_unavailable'; return result
        if result['source_frequency_difference'] > 0.15:
            result['status'] = 'palindromic_frequency_disagreement'; return result
    result['status'] = 'eligible'
    return result


def unique_alias_rows(rows):
    """Only collapse exact measurements whose sole difference is an rsID alias."""
    unique = {tuple(sorted((k,v) for k,v in row.items() if k != 'rsid')) for row in rows}
    if len(unique) != 1:
        raise ValueError('Conflicting source measurements')
    return rows[0], len(rows), ';'.join(sorted({r['rsid'] for r in rows}))


def main():
    plan = verify_plan()
    runs = {kind: load_queries(kind) for kind in ['gwas','qtl','reference']}
    if [len(runs[k]) for k in ['gwas','qtl','reference']] != [2,20,2]:
        raise ValueError('Complete regional data are required')
    forward = ChainMap.from_file(ROOT/'data/raw/coloc-hg19-to-hg38-chain.gz')
    reverse = ChainMap.from_file(ROOT/'data/raw/coloc-hg38-to-hg19-chain.gz')
    output = ROOT/'data/derived/coloc-inputs'; output.mkdir(exist_ok=True)
    local = ROOT/'data/raw/coloc-ld'; local.mkdir(exist_ok=True)
    gwas_audit, reference_audit, qtl_audit, dataset_audit, cs_rows = [], [], [], [], []
    for region in plan['regions']:
        name = region['gene_name']
        ref_file = ROOT/runs['reference']['EUR-'+name]['file']
        ref = pd.read_csv(ref_file, sep='\t', dtype={'chromosome':str})
        by_key = defaultdict(list)
        for i, row in ref.iloc[:, :5].iterrows():
            mapped = reciprocal_base(forward, reverse, 'chr'+row.chromosome, int(row.position_grch37)-1)
            item = dict(row); item['gene_name'] = name
            if mapped is None:
                item['status'] = 'nonunique_or_nonreciprocal_mapping'; reference_audit.append(item); continue
            chrom, pos0, strand = mapped
            if not region['start_grch38_0based'] <= pos0 < region['end_grch38_exclusive']:
                item['status'] = 'outside_fixed_grch38_region'; reference_audit.append(item); continue
            a, b = (complement(row.ref), complement(row.alt)) if strand == '-' else (row.ref, row.alt)
            key, _, effect = canonical_key(chrom, pos0+1, a, b)
            values = ref.iloc[i,5:].to_numpy(dtype=np.int8)
            called = values >= 0
            frequency = float(values[called].sum()/(2*called.sum())) if called.any() else float('nan')
            if b != effect: frequency = 1-frequency
            item.update(status='mapped', snp=key, position=pos0+1, effect_allele=effect,
                        effect_frequency=frequency, called_donors=int(called.sum()),
                        flip_reference_dosage=b != effect, row_index=i)
            by_key[key].append(item)
        reference = {}
        for key, entries in by_key.items():
            if len(entries) != 1:
                for entry in entries: entry['status'] = 'duplicate_reference_mapping'
            else: reference[key] = entries[0]
            reference_audit.extend(entries)
        with gzip.open(ROOT/runs['gwas']['GWAS-'+name]['file'],'rt') as stream:
            lines = stream.read().splitlines()
        fields = lines[0].lstrip('#').split()
        audited = [eligible_gwas(dict(zip(fields,line.split())),region,forward,reverse,reference) for line in lines[1:]]
        grouped = defaultdict(list)
        for row in audited:
            if row['status'] == 'eligible': grouped[row['snp']].append(row)
        for group in grouped.values():
            if len(group) > 1:
                for row in group: row['status'] = 'duplicate_gwas_mapping'
        eligible = [r for r in audited if r['status'] == 'eligible']
        gwas_audit.extend(audited)
        frame = pd.DataFrame([{'snp':r['snp'], 'position':r['position'], 'pvalue':float(r['p']),
                               'MAF':r['control_maf'], 'N':plan['gwas']['n'], 'canonical_sign':r['canonical_sign'],
                               'source_id':r['SNP'], 'source_or':r['or'], 'source_se':r['se']} for r in eligible]).sort_values(['position','snp'])
        save_frame(frame,output/('GWAS-'+name+'.csv.gz'))
        gwas_keys = set(frame.snp)
        dataset_audit.append({'gene_name':name,'dataset_id':'GWAS','source_rows':len(audited), 'eligible_unique_snvs':len(frame),
                              'statuses':dict(Counter(r['status'] for r in audited))})
        qtl_frames = {}
        for dataset in plan['datasets']:
            ds=dataset['dataset_id']
            with gzip.open(ROOT/runs['qtl'][ds+'-'+name]['file'],'rt') as stream:
                rows=list(csv.DictReader(stream,delimiter='\t'))
            grouped=defaultdict(list)
            for row in rows: grouped[row['variant']].append(row)
            kept=[]
            for variant, group in grouped.items():
                item={'gene_name':name,'dataset_id':ds,'source_variant':variant,'source_alias_rows':len(group),'status':'','snp':''}
                try:
                    row,aliases,rsids=unique_alias_rows(group)
                    if gene_id(row['gene_id'])!=region['gene_id'] or gene_id(row['molecular_trait_id'])!=region['gene_id']:
                        raise ValueError('Molecular trait is not the exact planned gene')
                    declared='chr'+row['chromosome']+'_'+row['position']+'_'+row['ref']+'_'+row['alt']
                    if variant!=declared:raise ValueError('Source variant fields disagree')
                    if row['chromosome']!=region['chromosome'][3:] or not region['start_grch38_0based']<int(row['position'])<=region['end_grch38_exclusive']:
                        raise ValueError('Source variant is outside the fixed region')
                    key,_,effect=canonical_key('chr'+row['chromosome'],int(row['position']),row['ref'],row['alt'])
                    numbers=measured_numbers(row,dataset['sample_size'])
                    maf=float(row['maf'])
                    if not math.isfinite(maf) or not 0<maf<=0.5:raise ValueError('Invalid QTL MAF')
                    beta=numbers['beta_per_alt']*(1 if row['alt']==effect else -1)
                    item.update(status='eligible',snp=key,gwas_overlap=key in gwas_keys,source_rsids=rsids,
                                source_beta=row['beta'],source_se=row['se'],source_pvalue=row['pvalue'])
                    kept.append({'snp':key,'position':int(row['position']),'beta':beta,'varbeta':numbers['standard_error']**2,
                                 'pvalue':numbers['nominal_pvalue'],'MAF':maf,'N':numbers['sample_size_from_an'],
                                 'source_variant':variant,'source_rsids':rsids,'source_beta':row['beta'],'source_se':row['se']})
                except (ValueError,TypeError) as exc:
                    item.update(status='unavailable',reason=str(exc))
                qtl_audit.append(item)
            qframe=pd.DataFrame(kept).sort_values(['position','snp'])
            if qframe.snp.duplicated().any():raise ValueError('Duplicate canonical QTL identity')
            save_frame(qframe,output/(ds+'-'+name+'.csv.gz'))
            qtl_frames[ds]=qframe
            dataset_audit.append({'gene_name':name,'dataset_id':ds,'source_rows':len(rows),'source_unique_variants':len(grouped),
                                  'eligible_unique_snvs':len(qframe),'gwas_overlap':int(qframe.snp.isin(gwas_keys).sum()),
                                  'sample_sizes':sorted(set(map(int,qframe.N)))})
        # Reference LD is an explicitly labelled sensitivity for GWAS and BLUEPRINT.
        for ds, trait in [('GWAS',frame),('QTD000031',qtl_frames['QTD000031'])]:
            selected=[]; vectors=[]
            for _,row in trait.iterrows():
                entry=reference.get(row.snp)
                if not entry or entry['called_donors']!=503 or not 0<entry['effect_frequency']<1:continue
                values=ref.iloc[entry['row_index'],5:].to_numpy(dtype=np.float64)
                if entry['flip_reference_dosage']:values=2-values
                selected.append(row.to_dict());vectors.append(values)
            selected=pd.DataFrame(selected)
            X=np.asarray(vectors,dtype=float)
            X-=X.mean(axis=1,keepdims=True)
            norm=np.linalg.norm(X,axis=1)
            if (norm==0).any():raise ValueError('Constant genotype in reference LD')
            X/=norm[:,None]
            ld=X@X.T
            if not np.isfinite(ld).all() or not np.allclose(ld,ld.T,atol=1e-12) or not np.allclose(np.diag(ld),1,atol=1e-12):
                raise ValueError('Invalid signed reference LD')
            # R reads column-major doubles; the symmetric matrix is identical in both orders.
            ld_path=local/(ds+'-'+name+'.f64')
            ld.astype('<f8').tofile(ld_path)
            save_frame(selected,output/(ds+'-'+name+'-ld.csv.gz'))
            save_json(local/(ds+'-'+name+'.json'),{'matrix_sha256':sha(ld_path.read_bytes()),'snp_file_sha256':sha((output/(ds+'-'+name+'-ld.csv.gz')).read_bytes()),
                                                       'n_variants':len(selected),'n_reference_donors':503,'coding':'signed r for lexical-maximum allele dosage',
                                                       'matrix_type':'Gram matrix of centered unit-norm donor dosage vectors; positive semidefinite by construction',
                                                       'reference_build':'GRCh37, uniquely reciprocally lifted to GRCh38', 'in_sample':False})
            dataset_audit.append({'gene_name':name,'dataset_id':ds+'-LD','eligible_unique_snvs':len(selected),'reference_donors':503})
    for dataset in plan['datasets']:
        ds=dataset['dataset_id']
        with gzip.open(ROOT/'data/raw'/('coloc-'+ds+'.credible_sets.tsv.gz'),'rt') as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                if row['molecular_trait_id'].split('.')[0] in {r['gene_id'] for r in plan['regions']}:
                    cs_rows.append({'dataset_id':ds,**row})
    save_frame(pd.DataFrame(gwas_audit),ROOT/'data/derived/coloc-gwas-alignment.csv.gz')
    save_frame(pd.DataFrame(reference_audit),ROOT/'data/derived/coloc-reference-alignment.csv.gz')
    save_frame(pd.DataFrame(qtl_audit),ROOT/'data/derived/coloc-qtl-alignment.csv.gz')
    save_frame(pd.DataFrame(cs_rows),ROOT/'data/derived/coloc-published-credible-sets.csv.gz')
    save_json(ROOT/'reports/coloc-input-audit.json',{'datasets':dataset_audit,'gwas_status_counts':dict(Counter(r['status'] for r in gwas_audit)),
              'qtl_status_counts':dict(Counter(r['status'] for r in qtl_audit)), 'reference_status_counts':dict(Counter(r['status'] for r in reference_audit)),
              'published_target_credible_set_rows':len(cs_rows),'plan_sha256':sha((ROOT/'config/coloc-plan.json').read_bytes()),
              'preparation_script_sha256':sha(__import__('pathlib').Path(__file__).read_bytes()),
              'inputs':{str(p.relative_to(ROOT)):sha(p.read_bytes()) for p in sorted(output.iterdir())}})
    print(json.dumps({'prepared_gene_dataset_pairs':20,'gwas_status_counts':dict(Counter(r['status'] for r in gwas_audit)),
                      'qtl_status_counts':dict(Counter(r['status'] for r in qtl_audit))},indent=2))


if __name__=='__main__':
    main()
