"""Align measured molecular QTLs to the recorded pilot without rewriting it."""

import csv
import gzip
import hashlib
import io
import json
import math
import os
from collections import Counter
from pathlib import Path

from fetch_qtl_followup import source_map, verified_queries
from fetch_sources import validate

ROOT = Path(__file__).resolve().parents[1]
UBASH = 'ENSG00000160185'
PRKD2 = 'ENSG00000105287'


def read_tsv(content):
    return list(csv.DictReader(io.StringIO(content.decode()), delimiter='\t'))


def write_csv(path, rows, fields=None):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields or list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def beta_for_allele(beta, effect_allele, ref, alt, requested_allele):
    """Flip only within an already build/strand-verified allele pair."""
    if ref == alt or effect_allele not in (ref, alt) or requested_allele not in (ref, alt):
        raise ValueError('Alleles do not match the verified pair; no complement guessing')
    beta = float(beta)
    if not math.isfinite(beta):
        raise ValueError('A missing or non-finite coefficient cannot establish direction')
    return beta if requested_allele == effect_allele else -beta


def select_gene_row(rows, variant, gene):
    matches = [r for r in rows if r['variant'] == variant and r['gene_id'] == gene
               and r['molecular_trait_id'] == gene]
    if len(matches) != 1:
        raise ValueError(f'Expected one gene-level row, found {len(matches)}')
    row = matches[0]
    chrom, position, ref, alt = variant.split('_')
    if (row['chromosome'], row['position'], row['ref'], row['alt']) != (chrom.removeprefix('chr'), position, ref, alt):
        raise ValueError('Variant identifier conflicts with explicit build-38 allele fields')
    return row


def junction_coordinates(phenotype):
    """LeafCutter regtools IDs use start0 and end0+1; keep source strand separate."""
    chromosome, start, end, cluster = phenotype.split(':')
    start0, end0 = int(start), int(end) - 1
    if end0 <= start0 or cluster.rsplit('_', 1)[-1] not in ('+', '-'):
        raise ValueError('Invalid LeafCutter junction identifier')
    return chromosome, start0, end0, cluster.rsplit('_', 1)[-1]


def strand_audit(rows):
    # Count each single-gene-assigned phenotype once, not each joined gene row.
    seen = set()
    counts = Counter()
    for row in rows:
        if row['gene_count'] != '1' or row['phenotype_id'] in seen:
            continue
        seen.add(row['phenotype_id'])
        source_strand = row['phenotype_id'].rsplit('_', 1)[-1]
        gene_strand = {'1': '+', '-1': '-'}.get(row['strand'])
        if source_strand in ('+', '-') and gene_strand:
            counts['concordant' if source_strand == gene_strand else 'discordant'] += 1
    total = counts['concordant'] + counts['discordant']
    return {'unique_single_gene_phenotypes': total, 'strand_concordant': counts['concordant'],
            'strand_discordant': counts['discordant'],
            'discordant_fraction': counts['discordant'] / total if total else ''}


def main():
    manifest = verified_queries()
    sources = source_map()

    def raw(source_id):
        source = sources[source_id]
        content = (ROOT / 'data/raw' / source['file']).read_bytes()
        validate(content, source)
        return content

    pilot = {r['rsid']: r for r in csv.DictReader((ROOT / 'reports/pilot-target-summary.csv').open())}
    allele_audit = {r['rsid']: r for r in csv.DictReader((ROOT / 'data/derived/allele-audit.csv').open())}
    # These signs are taken from the original artifact; expose any schema drift.
    expression, availability, junctions, annotations, strands = [], [], [], [], []
    for query in manifest['queries']:
        meta = query['dataset']
        dataset = meta['dataset_id']
        rows = read_tsv((ROOT / query['file']).read_bytes())
        if len(rows) != query['rows']:
            raise ValueError('Subset row count mismatch')
        variant = query['query']['variant']
        target_gene = PRKD2 if variant.startswith('chr19_') else UBASH
        target_rows = [r for r in rows if r['variant'] == variant and r['gene_id'] == target_gene]
        availability.append({'dataset_id': dataset, 'study': meta['study_label'], 'sample_group': meta['sample_group'],
                             'quant_method': meta['quant_method'], 'sample_size': meta['sample_size'],
                             'queried_variant': variant, 'chromosome_in_index': query['contig_present'],
                             'returned_rows_at_position': len(rows), 'target_gene_variant_rows': len(target_rows),
                             'absence_interpretation': 'not a null result; cc files retain selected tag traits' if meta['quant_method'] == 'leafcutter' else 'all nominal gene-level rows queried'})
        if meta['quant_method'] == 'ge':
            row = select_gene_row(rows, variant, target_gene)
            n = int(row['an']) // 2
            if int(row['an']) % 2 or n != int(meta['sample_size']):
                raise ValueError('Allele number disagrees with autosomal sample count')
            beta, se = float(row['beta']), float(row['se'])
            if not math.isfinite(se) or se <= 0 or not 0 < float(row['pvalue']) <= 1:
                raise ValueError('Invalid uncertainty or p value')
            risk = allele_audit[row['rsid']]['working_risk_allele']
            recorded_model = pilot[row['rsid']]
            if recorded_model['allele_change'] != row['ref']+'>'+row['alt']:
                raise ValueError('Pilot and QTL allele changes differ')
            model_lfc = float(recorded_model['median_signed_lfc'])
            if not math.isfinite(model_lfc) or model_lfc == 0:
                raise ValueError('Pilot lacks a usable signed prediction')
            model_sign = 1 if model_lfc > 0 else -1
            expression.append({'dataset_id': dataset, 'study': meta['study_label'], 'cell_context': meta['sample_group'],
                               'gene': 'PRKD2' if target_gene == PRKD2 else 'UBASH3A', 'gene_id': target_gene,
                               'rsid': row['rsid'], 'variant_grch38': variant, 'ref': row['ref'], 'effect_allele': row['alt'],
                               'beta_per_alt': beta, 'standard_error': se, 'wald_95_lower': beta-1.96*se,
                               'wald_95_upper': beta+1.96*se, 'nominal_pvalue': row['pvalue'], 'sample_size': n,
                               'alt_allele_count': row['ac'], 'allele_number': row['an'], 'minor_allele_frequency': row['maf'],
                               'imputation_r2': row['r2'], 'working_psc_risk_allele': risk,
                               'beta_per_working_risk': beta_for_allele(beta, row['alt'], row['ref'], row['alt'], risk),
                               'model_alt_direction': 'increase' if model_sign > 0 else 'decrease',
                               'recorded_model_median_signed_lfc': model_lfc,
                               'direction_agrees': beta * model_sign > 0,
                               'comparison_scope': 'direction only; normalized QTL beta is not model log-fold change',
                               'source_file': query['file']})
        else:
            metadata_rows = read_tsv(gzip.decompress(raw(dataset+'-phenotype-metadata')))
            by_id = {r['phenotype_id']: r for r in metadata_rows if r['gene_id'] == UBASH}
            strands.append({'dataset_id': dataset, **strand_audit(metadata_rows)})
            for row in by_id.values():
                annotations.append({'dataset_id': dataset, **row})
            for row in target_rows:
                annotation = by_id[row['molecular_trait_id']]
                chromosome, start0, end0, source_strand = junction_coordinates(row['molecular_trait_id'])
                gene_strand = '+' if annotation['strand'] == '1' else '-'
                plus29 = (chromosome, start0, end0) == ('21', 42434983, 42437487)
                direction_match = float(row['beta']) > 0 if plus29 else None
                interpretation = 'different junction; not a test of the plus29 event'
                if plus29:
                    interpretation = 'same boundary; '+('direction agrees' if direction_match else 'direction does not agree')
                    if source_strand != gene_strand:
                        interpretation += '; strand unresolved'
                junctions.append({'dataset_id': dataset, 'study': meta['study_label'], 'sample_size': meta['sample_size'],
                                  'variant': row['variant'], 'rsid': row['rsid'], 'phenotype_id_original': row['molecular_trait_id'],
                                  'chromosome': chromosome, 'start0': start0, 'end0': end0,
                                  'cluster_strand_original': source_strand, 'annotated_gene_strand': gene_strand,
                                  'strand_agrees': source_strand == gene_strand, 'matches_plus29_boundaries': plus29,
                                  'plus29_model_direction_agrees': direction_match,
                                  'effect_allele': row['alt'], 'beta_per_alt': row['beta'], 'se': row['se'], 'nominal_pvalue': row['pvalue'],
                                  'interpretation': interpretation})

    reports = ROOT / 'reports'
    write_csv(reports/'qtl-expression-comparison.csv', expression)
    write_csv(reports/'qtl-query-coverage.csv', availability)
    write_csv(reports/'qtl-splice-comparison.csv', junctions)
    write_csv(reports/'qtl-strand-audit.csv', strands)
    write_csv(ROOT/'data/derived/ubash3a-qtl-junction-annotations.csv', annotations)

    import openpyxl
    workbook = openpyxl.load_workbook(io.BytesIO(raw('alphagenome-supplement-tables')), read_only=True, data_only=True)
    sheet = workbook['Suppl Table 2 Track metadata (f']
    iterator = iter(sheet.values)
    header = next(iterator)
    track_rows = []
    for values in iterator:
        row = {key: value for key, value in zip(header, values) if key is not None}
        if row['organism'] == 'human' and row['ontology_curie'] in ('CL:0000624', 'CL:0001054') and row['output_type'] in ('RNA_SEQ', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS'):
            track_rows.append(row)
    write_csv(ROOT/'data/derived/qtl-comparison-training-tracks.csv', track_rows)
    provenance = {'query_manifest_sha256': hashlib.sha256((ROOT/'config/qtl-followup-queries.json').read_bytes()).hexdigest(),
                  'source_manifest_sha256': hashlib.sha256((ROOT/'config/qtl-followup-sources.json').read_bytes()).hexdigest(),
                  'pilot_summary_sha256': hashlib.sha256((ROOT/'reports/pilot-target-summary.csv').read_bytes()).hexdigest(),
                  'allele_audit_sha256': hashlib.sha256((ROOT/'data/derived/allele-audit.csv').read_bytes()).hexdigest(),
                  'queries': len(availability), 'preserved_query_rows': sum(q['rows'] for q in manifest['queries']),
                  'expression_comparisons': len(expression), 'ubash3a_splice_rows': len(junctions),
                  'annotation_rows': len(annotations), 'training_track_rows': len(track_rows),
                  'beta_units': 'inverse-normalized molecular phenotype per ALT allele; not a percent or fold-change',
                  'confidence_intervals': 'descriptive Wald beta +/- 1.96 SE; no cross-study meta-analysis',
                  'independence': 'BLUEPRINT reuses a benchmark cohort; DICE conditions share donors; genome intervals were not held out; exact training donor overlap unresolved',
                  'junction_conversion': 'regtools BED12 to LeafCutter: A=start+firstBlockSize; B=end-lastBlockSize+1. Convert to [A,B-1).',
                  'strand_policy': 'Preserve cluster and annotated gene strands separately; no silent flip',
                  'no_new_alphagenome_requests': True}
    (reports/'qtl-followup-provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    make_figure(expression)
    print(json.dumps(provenance, indent=2))


def make_figure(rows):
    os.environ.setdefault('MPLCONFIGDIR', str(ROOT/'work/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'svg.hashsalt': 'psc-qtl-followup'})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7))
    fig.patch.set_facecolor('#faf9f5')
    for ax, gene, color, expected in zip(axes, ['PRKD2','UBASH3A'], ['#16746e','#b25439'], ['increase','decrease']):
        selected = [row for row in rows if row['gene'] == gene]
        labels = ['BLUEPRINT monocytes (191)', 'DICE monocytes (91)'] if gene == 'PRKD2' else ['DICE resting CD4 (88)', 'DICE activated CD4 (89)']
        for y, row in enumerate(selected):
            ax.errorbar(row['beta_per_alt'], y, xerr=1.96*row['standard_error'], fmt='o', color=color, capsize=5, markersize=8, linewidth=2)
        ax.axvline(0, color='#999999', linewidth=1)
        ax.set_yticks(range(len(selected)), labels, fontsize=9)
        ax.set_ylim(1.7,-.7)
        ax.set_xlim(-.15,1.55)
        ax.set_title(gene+'\nModel direction: '+expected, loc='left', fontweight='bold', fontsize=13, pad=16)
        ax.set_xlabel('Measured effect per ALT allele\n(normalized RNA; 95% Wald interval)', fontsize=10)
        ax.set_facecolor('#faf9f5')
        for side in ['top','right','left']: ax.spines[side].set_visible(False)
        ax.grid(axis='x', alpha=.16)
    fig.suptitle('Measured RNA: PRKD2 direction agrees; UBASH3A disagrees', fontsize=15, x=.05, ha='left', y=.98)
    fig.text(.05,.02,'ALT = G for PRKD2; C for UBASH3A. DICE conditions share donors. Molecular associations do not establish treatment effects.', fontsize=9, color='#555555')
    fig.tight_layout(rect=(0,.09,1,.90), w_pad=2)
    for suffix in ['png','svg']:
        path = ROOT/'reports'/('qtl-expression-comparison.'+suffix)
        fig.savefig(path, dpi=160, facecolor=fig.get_facecolor(), metadata={'Date':None} if suffix=='svg' else None)
        if suffix == 'svg': path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)


if __name__ == '__main__':
    main()
