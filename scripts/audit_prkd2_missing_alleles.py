"""Check omitted RNA edits against every regional GWAS source allele record.

This is an absence/identity audit, not a rescue into any statistical input.
Both possible source REF orientations are retained when the reference permits
them, because the GWAS labels are non-risk/risk rather than REF/ALT.
"""

import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from prkd2_common import ROOT, ChainMap, save_json, sha
from prkd2_expanded_common import load_references, normalize_edit, qtl_identity
from coverage_repair_common import reciprocal_edit
from prepare_coloc_inputs import save_frame


def candidate_edits(rows, reference):
    candidates = defaultdict(set)
    counts = Counter()
    for row in rows:
        counts['all_source_rows_considered'] += 1
        a, b = row['allele_0'], row['allele_1']
        counts['SNV_rows' if len(a) == len(b) == 1 else 'non_SNV_rows'] += 1
        for ref, alt in [(a, b), (b, a)]:
            counts['orientations_considered'] += 1
            try:
                edit = normalize_edit(reference, int(row['pos']), ref, alt)
            except ValueError:
                counts['reference_invalid_orientations'] += 1
                continue
            candidates[tuple(edit)].add((int(row['source_row']), row['SNP']))
            counts['reference_valid_orientations'] += 1
    return candidates, dict(counts)


def main():
    refs = load_references()
    forward = ChainMap.from_file(ROOT / 'data/raw/coloc-hg19-to-hg38-chain.gz')
    reverse = ChainMap.from_file(ROOT / 'data/raw/coloc-hg38-to-hg19-chain.gz')
    source = ROOT / 'data/derived/prkd2-gwas-alignment.csv.gz'
    omitted = ROOT / 'data/derived/prkd2-missing-evidence.csv'
    raw = pd.read_csv(source, dtype=str, keep_default_na=False)
    candidates, counts = candidate_edits(raw.to_dict('records'), refs['GRCh37'])
    selected = pd.read_csv(omitted).query("trait == 'QTL'")
    rows = []
    for original in selected.to_dict('records'):
        item = {**original, 'genome_build': 'GRCh38', 'full_source_GWAS_possible_matches': ''}
        try:
            edit = qtl_identity(original['source_id'], refs['GRCh38'])
            mapped, bases, strand = reciprocal_edit(refs['GRCh38'], refs['GRCh37'], reverse, forward,
                    edit['position'], edit['normalized_ref'], edit['normalized_alt'], flank=100)
            normalized = normalize_edit(refs['GRCh37'], *mapped)
            matches = sorted(candidates.get(tuple(normalized), set()))
            item.update(status='possible_full_allele_source_match' if matches else 'no_full_allele_source_match',
                full_source_GWAS_possible_matches=';'.join(f'{i}:{name}' for i, name in matches),
                mapped_grch37=json.dumps(normalized), reciprocal_bases=bases, strand=strand)
        except ValueError as exc:
            item.update(status='identity_unavailable', reason=str(exc))
        rows.append(item)
    output = ROOT / 'data/derived/prkd2-missing-allele-identity-audit.csv'
    save_frame(pd.DataFrame(rows), output)
    record = {'scope': 'All top-20 omitted RNA edits in each cohort, compared with all regional GWAS source rows, including SNVs and non-SNVs',
        'statistical_inputs_changed': False, 'candidate_counts': counts, 'target_rows': len(rows),
        'status_counts': dict(Counter(r['status'] for r in rows)),
        'interpretation': 'Possible matches are candidates only; no match means no full-allele source record under these exact reference/reciprocal mapping rules. Filter causes are not inferred.',
        'development_correction': 'The preliminary uncommitted audit considered only non-SNV GWAS rows. This complete audit includes every source row; preliminary output is retained in the local work directory and was not used to change statistical inputs.',
        'inputs': {str(p.relative_to(ROOT)): sha(p.read_bytes()) for p in [source, omitted, ROOT / 'config/prkd2-expanded-sources.json']},
        'output': str(output.relative_to(ROOT)), 'output_sha256': sha(output.read_bytes()),
        'reader_sha256': sha(Path(__file__).read_bytes())}
    save_json(ROOT / 'reports/prkd2-missing-allele-verification.json', record)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
