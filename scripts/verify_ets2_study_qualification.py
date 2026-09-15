from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile,json,re,hashlib,math
ROOT=Path(__file__).resolve().parents[1]
raw=ROOT/'data/raw/ets2-study-qualification'
main=json.loads((raw/'matrix-structural-audit.json').read_text())
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
source=raw/'GSE255234_raw_counts_for_rna_seq_submission.xlsx'
with zipfile.ZipFile(source) as z:
    shared=[]
    if 'xl/sharedStrings.xml' in z.namelist():
        tree=ET.fromstring(z.read('xl/sharedStrings.xml'))
        shared=[''.join(n.itertext()) for n in tree.findall(NS+'si')]
    gene_ids=[]; header=[]; sums=[0]*16; numeric_cells=0; allzero=0; zero=0
    with z.open('xl/worksheets/sheet1.xml') as f:
        for _,row in ET.iterparse(f,events=['end']):
            if row.tag!=NS+'row': continue
            values=[]
            for cell in row.findall(NS+'c'):
                typ=cell.get('t'); v=cell.find(NS+'v')
                if typ=='s': value=shared[int(v.text)]
                elif typ=='inlineStr': value=''.join(cell.find(NS+'is').itertext())
                elif v is None: value=None
                else: value=float(v.text)
                col_letters=re.match(r'[A-Z]+',cell.attrib['r']).group()
                col_index=0
                for letter in col_letters: col_index=col_index*26+ord(letter)-ord('A')+1
                while len(values)<col_index: values.append(None)
                values[col_index-1]=value
            if row.get('r')=='1':
                header=values
            else:
                gene_ids.append(values[0]); numbers=values[1:]
                assert len(numbers)==16
                assert all(isinstance(n,float) and math.isfinite(n) and n>=0 and n.is_integer() for n in numbers)
                numeric_cells+=len(numbers); allzero+=all(n==0 for n in numbers); zero+=sum(n==0 for n in numbers)
                sums=[s+n for s,n in zip(sums,numbers)]
            row.clear()
result={
 'method':'Independent stdlib ZIP + ElementTree XML scan; no openpyxl',
 'gene_rows':len(gene_ids),'unique_ids':len(set(gene_ids)),
 'all_ids_ensembl':all(re.fullmatch('ENSG[0-9]{11}',g) for g in gene_ids),
 'all_zero_rows':allzero,'zero_cells':zero,'numeric_cells':numeric_cells,
 'header':header,
 'source_id_set_sorted_sha256':hashlib.sha256(('\n'.join(sorted(gene_ids))+'\n').encode()).hexdigest(),
 'library_sums_agree':sums==[r['library_total'] for r in main['sample_structural_qc']],
 'all_numeric_finite_nonnegative_integer':True,
}
assert result['gene_rows']==main['row_count']
assert result['unique_ids']==main['unique_identifiers']
assert allzero==main['all_zero_features'] and zero==main['zero_cells']
assert result['source_id_set_sorted_sha256']==main['source_id_set_sorted_sha256']
assert result['library_sums_agree'] and numeric_cells==main['numeric_cells']
assert header==main['header']
result['all_primary_comparisons_passed']=True
(raw/'matrix-independent-xml-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
