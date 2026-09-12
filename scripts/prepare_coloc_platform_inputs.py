"""Apply the documented platform-count amendment without altering the first run."""

import json
from collections import Counter

import numpy as np
import pandas as pd

from coloc_common import ROOT, ChainMap, save_json, sha, verify_plan
from prepare_coloc_inputs import eligible_gwas, save_frame


def platform_gwas(raw, region, forward, reverse, reference, counts):
    platform=raw['platform']
    # The original helper's only platform-dependent operation is its exclusion.
    # Re-run all its coordinate, allele, frequency and numeric checks unchanged.
    audited=eligible_gwas({**raw,'platform':'both'},region,forward,reverse,reference)
    audited['platform']=platform
    if platform not in counts:
        audited['status']='platform_sample_size_unavailable'
    else:
        audited.update(N=counts[platform]['N'],cases=counts[platform]['cases'],controls=counts[platform]['controls'])
    return audited


def main():
    original=verify_plan()
    plan_path=ROOT/'config/coloc-platform-plan.json'
    plan=json.loads(plan_path.read_text())
    for name,digest in plan['frozen_file_sha256'].items():
        if sha((ROOT/name).read_bytes())!=digest:raise ValueError('Platform amendment input changed')
    forward=ChainMap.from_file(ROOT/'data/raw/coloc-hg19-to-hg38-chain.gz')
    reverse=ChainMap.from_file(ROOT/'data/raw/coloc-hg38-to-hg19-chain.gz')
    previous=pd.read_csv(ROOT/'data/derived/coloc-gwas-alignment.csv.gz',dtype=str,keep_default_na=False)
    refs=pd.read_csv(ROOT/'data/derived/coloc-reference-alignment.csv.gz',dtype={'flip_reference_dosage':'string'},low_memory=False)
    output=ROOT/'data/derived/coloc-platform-inputs';output.mkdir(exist_ok=True)
    local=ROOT/'data/raw/coloc-platform-ld';local.mkdir(exist_ok=True)
    all_rows=[];metadata=[]
    for region in original['regions']:
        gene=region['gene_name']
        entries=refs[(refs.gene_name==gene)&(refs.status=='mapped')]
        if entries.snp.duplicated().any():raise ValueError('Ambiguous reference mapping')
        reference={r['snp']:r for r in entries.to_dict('records')}
        source_columns=['chr','SNP','pos','allele_0','allele_1','freq_1','freq_1_cases','freq_1_controls','mmm_var_info_nonmissing','or','se','p','platform']
        rows=[platform_gwas(row,region,forward,reverse,reference,plan['platform_counts']) for row in previous[previous.gene_name==gene][source_columns].to_dict('records')]
        keys=Counter(r['snp'] for r in rows if r['status']=='eligible')
        for r in rows:
            if r['status']=='eligible' and keys[r['snp']]>1:r['status']='duplicate_gwas_mapping'
        all_rows.extend(rows)
        frame=pd.DataFrame([{'snp':r['snp'],'position':r['position'],'pvalue':float(r['p']),'MAF':r['control_maf'],
                             'N':r['N'],'cases':r['cases'],'controls':r['controls'],'case_fraction':r['cases']/r['N'],
                             'canonical_sign':r['canonical_sign'],'source_id':r['SNP'],'platform':r['platform'],
                             'source_or':r['or'],'source_se':r['se']} for r in rows if r['status']=='eligible']).sort_values(['position','snp'])
        save_frame(frame,output/('GWAS-'+gene+'.csv.gz'))
        dosages=pd.read_csv(ROOT/'data/raw'/('coloc-EUR-'+gene+'-dosages.tsv.gz'),sep='\t')
        selected=[];vectors=[]
        for _,row in frame.iterrows():
            entry=reference.get(row.snp)
            if not entry or entry['called_donors']!=503 or not 0<entry['effect_frequency']<1:continue
            values=dosages.iloc[int(entry['row_index']),5:].to_numpy(dtype=float)
            if entry['flip_reference_dosage'] not in ('True','False'):raise ValueError('Missing reference dosage orientation')
            if entry['flip_reference_dosage']=='True':values=2-values
            selected.append(row.to_dict());vectors.append(values)
        selected=pd.DataFrame(selected); X=np.array(vectors,dtype=float);X-=X.mean(axis=1,keepdims=True)
        norm=np.linalg.norm(X,axis=1)
        if (norm==0).any():raise ValueError('Constant reference genotype')
        X/=norm[:,None];ld=X@X.T
        if not np.isfinite(ld).all() or not np.allclose(ld,ld.T,atol=1e-12) or not np.allclose(np.diag(ld),1,atol=1e-12):raise ValueError('Invalid reference correlation')
        matrix_path=local/('GWAS-'+gene+'.f64');ld.astype('<f8').tofile(matrix_path)
        save_frame(selected,output/('GWAS-'+gene+'-ld.csv.gz'))
        item={'gene_name':gene,'raw_gwas_rows':len(rows),'eligible_unique_snvs':len(frame),'reference_LD_variants':len(selected),
              'status_counts':dict(Counter(r['status'] for r in rows)), 'platform_counts_in_input':dict(Counter(frame.platform)),
              'matrix_file':str(matrix_path.relative_to(ROOT)),'matrix_sha256':sha(matrix_path.read_bytes()),'reference_donors':503}
        metadata.append(item)
    save_frame(pd.DataFrame(all_rows),ROOT/'data/derived/coloc-platform-gwas-alignment.csv.gz')
    save_json(ROOT/'reports/coloc-platform-input-audit.json',{'amendment_sha256':sha(plan_path.read_bytes()),'datasets':metadata,
              'input_sha256':{str(p.relative_to(ROOT)):sha(p.read_bytes()) for p in sorted(output.iterdir())},
              'preparation_script_sha256':sha(__import__('pathlib').Path(__file__).read_bytes())})
    print(json.dumps(metadata,indent=2))


if __name__=='__main__':main()
