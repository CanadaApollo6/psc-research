"""Independent Python ABF, coverage, PIP and posterior checks of the R analyses."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp, ndtri

from prkd2_common import ROOT, save_json, sha


def binary_bf(frame):
    variance = 1 / (2 * frame.N.to_numpy(float) * frame.MAF.to_numpy(float) * (1 - frame.MAF.to_numpy(float))
                    * frame.case_fraction.to_numpy(float) * (1 - frame.case_fraction.to_numpy(float)))
    shrinkage = .04 / (.04 + variance)
    z = -ndtri(frame.pvalue.to_numpy(float) / 2)
    return .5 * (np.log1p(-shrinkage) + shrinkage * z ** 2)


def expression_sd(frame):
    if len(set(frame.N)) != 1:
        raise ValueError('Variable sample size')
    precision = 1 / frame.varbeta.to_numpy(float)
    dosage_variance_n = 2 * frame.N.to_numpy(float) * frame.MAF.to_numpy(float) * (1 - frame.MAF.to_numpy(float))
    return float(np.sqrt(np.dot(precision, dosage_variance_n) / np.dot(precision, precision)))


def quantitative_bf(frame, sd=1):
    variance = frame.varbeta.to_numpy(float)
    prior = (.15 * sd) ** 2
    shrinkage = prior / (prior + variance)
    return .5 * (np.log1p(-shrinkage) + shrinkage * frame.beta.to_numpy(float) ** 2 / variance)


def weights(log_bf):
    log_bf = np.asarray(log_bf, float)
    if np.isnan(log_bf).any() or np.isposinf(log_bf).any() or not np.isfinite(logsumexp(log_bf)):
        raise ValueError('Invalid full evidence vector')
    return np.exp(log_bf - logsumexp(log_bf))


def posterior(log1, log2, p12, p1=1e-4, p2=1e-4):
    a, b = np.asarray(log1, float), np.asarray(log2, float)
    if a.shape != b.shape or a.ndim != 1:
        raise ValueError('Aligned one-dimensional vectors required')
    joint = logsumexp(a + b)
    separate = logsumexp(a) + logsumexp(b)
    if separate < joint - 1e-12:
        raise ValueError('Impossible separate-variant evidence')
    difference = separate + np.log(-np.expm1(min(0., joint - separate))) if separate > joint else -np.inf
    logs = np.array([0, np.log(p1) + logsumexp(a), np.log(p2) + logsumexp(b),
                     np.log(p1) + np.log(p2) + difference, np.log(p12) + joint])
    return np.exp(logs - logsumexp(logs))


def near(actual, expected, label, tolerance=2e-10):
    a, e = np.asarray(actual, float), np.asarray(expected, float)
    if a.shape != e.shape or not np.allclose(a, e, atol=tolerance, rtol=tolerance, equal_nan=True):
        raise ValueError('Independent calculation disagrees: ' + label)
    return float(np.nanmax(np.abs(a - e))) if a.size and not np.isnan(a).all() else 0.


def verify_baseline(prefix):
    folder = ROOT / ('data/derived/' + prefix + '-inputs')
    g = pd.read_csv(folder / 'GWAS-PRKD2.csv.gz')
    gf = pd.read_csv(folder / 'GWAS-full-evidence.csv.gz')
    table = pd.read_csv(ROOT / ('reports/' + prefix + '-baseline.csv'))
    details = pd.read_csv(ROOT / ('reports/' + prefix + '-baseline-variant-posteriors.csv.gz'))
    assert len(table) == 12
    bg, bgf = binary_bf(g), binary_bf(gf)
    max_error, variants = 0., 0
    for ds in ['QTD000021', 'QTD000504']:
        q = pd.read_csv(folder / (ds + '-PRKD2.csv.gz'))
        qf = pd.read_csv(folder / (ds + '-full-evidence.csv.gz'))
        common = sorted(set(g.snp) & set(q.snp))
        gi = g.set_index('snp').index.get_indexer(common)
        qi = q.set_index('snp').index.get_indexer(common)
        gmask = gf.source_row.isin(g.iloc[gi].source_row).to_numpy()
        qmask = qf.source_variant.isin(q.iloc[qi].source_variant).to_numpy()
        for mode in ['unit', 'estimated']:
            sd = 1 if mode == 'unit' else expression_sd(qf)
            bq, bqf = quantitative_bf(q, sd), quantitative_bf(qf, sd)
            gm, qm = float(weights(bgf)[gmask].sum()), float(weights(bqf)[qmask].sum())
            gate = len(common) >= 100 and gm >= .9 and qm >= .9
            for prior in [1e-6, 1e-7, 1e-5]:
                selected = table[(table.dataset_id == ds) & (table.sd_mode == mode) & np.isclose(table.p12, prior, atol=0, rtol=1e-12)]
                assert len(selected) == 1
                row = selected.iloc[0]
                assert row.nsnps == len(common) and bool(row.coverage_gate) == gate
                max_error = max(max_error, near([row.full_gwas_bf_mass_retained, row.full_qtl_bf_mass_retained, row.sdY_used], [gm, qm, sd], prefix + ds + mode))
                max_error = max(max_error, near(row[[f'PP.H{i}.abf' for i in range(5)]], posterior(bg[gi], bq[qi], prior), prefix + ds + mode + str(prior)))
            if mode == 'unit':
                observed = details[details.dataset_id == ds].set_index('snp').loc[common]
                near(observed.gwas_logbf, bg[gi], 'GWAS log BFs')
                near(observed.qtl_logbf, bq[qi], 'QTL log BFs')
                near(observed.conditional_H4_snp_probability, weights(bg[gi] + bq[qi]), 'Conditional variant probabilities')
                variants += len(common)
    return {'rows': len(table), 'variant_posteriors': variants, 'maximum_numeric_difference': max_error}


def published_pips(ds):
    raw = pd.read_csv(ROOT / ('data/derived/prkd2-query-rows/' + ds + '-PRKD2-published-bfs.tsv.gz'), sep='\t')
    bfcols = [c for c in raw.columns if c.startswith('lbf_variable')]
    if raw.variant.duplicated().any():
        raise ValueError('Duplicated full published variant')
    alpha = np.stack([weights(raw[c]) for c in bfcols], axis=1)
    pips = 1 - np.prod(1 - alpha, axis=1)
    cs = pd.read_csv(ROOT / 'data/derived/prkd2-published-credible-sets.csv.gz')
    cs = cs[cs.dataset_id == ds]
    indices = raw.set_index('variant').index.get_indexer(cs.variant)
    if (indices < 0).any():
        raise ValueError('Published credible-set variant absent from full vectors')
    error = near(cs.pip, pips[indices], ds + ' published marginal PIPs', tolerance=1e-8)
    return {'full_variants': len(raw), 'published_credible_set_rows': len(cs),
            'published_distinct_CS_variants': int(cs.variant.nunique()), 'components': sorted(set(cs.cs_id)),
            'maximum_PIP_difference': error}


def component_vectors(prefix, ds, n):
    gf = pd.read_csv(ROOT / ('data/derived/' + prefix + '-GWAS-component-bfs-N' + str(n) + '.csv.gz'))
    record = json.loads((ROOT / ('reports/' + prefix + '-GWAS-fit-N' + str(n) + '.json')).read_text())
    gcols = list(gf.columns[1:])
    gb = {int(cs['component']): gf[gcols[int(cs['component']) - 1]].to_numpy(float) for cs in record['credible_sets']}
    if prefix == 'prkd2-expanded':
        qf = pd.read_csv(ROOT / ('data/derived/prkd2-expanded-published-bfs/' + ds + '.csv.gz'))
        qnames = list(qf.snp)
    else:
        qf = pd.read_csv(ROOT / ('data/derived/prkd2-query-rows/' + ds + '-PRKD2-published-bfs.tsv.gz'), sep='\t')
        qnames = []
        for variant in qf.variant:
            chrom, pos, ref, alt = variant.split('_')
            qnames.append(':'.join([chrom, pos, *sorted([ref, alt])]) if len(ref) == len(alt) == 1 else 'unmatched_nonSNV:' + variant)
    cs = pd.read_csv(ROOT / 'data/derived/prkd2-published-credible-sets.csv.gz')
    ids = sorted({int(x.rsplit('_L', 1)[1]) for x in cs[cs.dataset_id == ds].cs_id})
    qb = {i: qf['lbf_variable' + str(i)].to_numpy(float) for i in ids}
    return list(gf.snp), gb, qnames, qb


def verify_components(prefix):
    table = pd.read_csv(ROOT / ('reports/' + prefix + '-signal-comparison.csv'))
    detail = pd.read_csv(ROOT / ('reports/' + prefix + '-signal-variant-posteriors.csv.gz'))
    expected_rows, max_error, checked_variants = 0, 0., 0
    for ds in ['QTD000021', 'QTD000504']:
        for n in [11386, 3504, 14890]:
            gn, gb, qn, qb = component_vectors(prefix, ds, n)
            common = sorted(set(gn) & set(qn))
            gi, qi = pd.Index(gn).get_indexer(common), pd.Index(qn).get_indexer(common)
            for i, a in gb.items():
                for j, b in qb.items():
                    gm, qm = float(weights(a)[gi].sum()), float(weights(b)[qi].sum())
                    gate = len(common) >= 100 and gm >= .9 and qm >= .9
                    for prior in [1e-6, 1e-7, 1e-5]:
                        expected_rows += 1
                        selected = table[(table.dataset_id == ds) & (table.n_approximation == n) & (table.gwas_component == i)
                                         & (table.qtl_component == j) & np.isclose(table.p12, prior, atol=0, rtol=1e-12)]
                        if len(selected) != 1:
                            raise ValueError('Missing or duplicated planned component comparison')
                        row = selected.iloc[0]
                        assert int(row.nsnps) == len(common) and bool(row.coverage_gate) == gate
                        near([row.gwas_signal_mass_retained, row.qtl_signal_mass_retained], [gm, qm], 'Full component overlap')
                        if gate:
                            max_error = max(max_error, near(row[[f'PP.H{k}.abf' for k in range(5)]], posterior(a[gi], b[qi], prior), 'Component posteriors'))
                        elif not np.isnan(row['PP.H4.abf']):
                            raise ValueError('A coverage-failed component was assigned an H4 posterior')
                    if gate and n == 11386:
                        observed = detail[(detail.dataset_id == ds) & (detail.gwas_component == i) & (detail.qtl_component == j)].set_index('snp').loc[common]
                        near(observed['SNP.PP.H4.abf'], weights(a[gi] + b[qi]), 'Component conditional variant probabilities')
                        checked_variants += len(common)
    assert len(table) == expected_rows
    return {'planned_rows_verified': expected_rows, 'variant_posteriors_verified': checked_variants, 'maximum_numeric_difference': max_error}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baselines-only', action='store_true')
    args = parser.parse_args()
    record = {'all_checks_pass': True, 'baselines': {p: verify_baseline(p) for p in ['prkd2', 'prkd2-expanded']},
              'verification_method': 'Independent NumPy/SciPy formulas; no R result used as an input coefficient or weight'}
    if not args.baselines_only:
        record['published_PIP_checks'] = {ds: published_pips(ds) for ds in ['QTD000021', 'QTD000504']}
        record['components'] = {p: verify_components(p) for p in ['prkd2', 'prkd2-expanded']}
    record['verifier_sha256'] = sha(Path(__file__).read_bytes())
    name = 'prkd2-baseline-verification.json' if args.baselines_only else 'prkd2-statistical-verification.json'
    save_json(ROOT / 'reports' / name, record)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
