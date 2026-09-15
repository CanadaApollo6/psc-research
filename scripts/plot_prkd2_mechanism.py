"""Fixed-window molecular predictions and measured peak context; no effect cropping."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd

from prkd2_mechanism_common import ROOT, sha

DATA = ROOT / 'data/derived/prkd2-mechanism-results'
RUN = ROOT / 'data/predictions/20260915T185853Z-prkd2-mechanism'
DEFAULT = ROOT / 'reports/figures/prkd2-mechanism'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.facecolor': 'white', 'axes.titleweight': 'bold', 'pdf.fonttype': 42})


def save(fig, output, name):
    for ext in ['png', 'pdf']:
        fig.savefig(output / (name + '.' + ext), dpi=160, bbox_inches='tight')
    plt.close(fig)


def overview(output):
    genes = pd.read_csv(DATA / 'primary-gene-ranks.csv')
    target = pd.read_csv(DATA / 'target-summary.csv')
    peaks = pd.read_csv(ROOT / 'data/derived/prkd2-mechanism-encode/peak-overlap-summary.csv')
    training = pd.read_csv(ROOT / 'data/derived/prkd2-mechanism-training/selected-peak-training-overlap.csv').set_index('experiment').fillna('')
    order = target.rsid.tolist()
    fig = plt.figure(figsize=(13, 7.7))
    grid = fig.add_gridspec(1, 2, width_ratios=[1, 1.65], left=.075, right=.97, bottom=.25, top=.79, wspace=.68)
    ax = fig.add_subplot(grid[0])
    for i, rsid in enumerate(order):
        group = genes[genes.rsid == rsid].sort_values('gene_id')
        ax.scatter(i + np.linspace(-.16, .16, len(group)), group.signed_rna_median, c='#9ba7ae', s=19, alpha=.65, linewidth=0)
        row = target[target.rsid == rsid].iloc[0]
        ax.scatter(i, row.signed_rna_median, c='#b33250', s=75, zorder=4, edgecolor='white')
        ax.annotate(f"PRKD2\n{row.signed_rna_median:+.5f}\nrank {int(row.model_rank)}/50",
                    (i, row.signed_rna_median), xytext=(0, 13 if i != 0 else -51), textcoords='offset points',
                    ha='center', color='#922640', fontsize=9, fontweight='bold')
    ax.axhline(0, c='#33414c', lw=.8)
    ax.set_xticks(range(3), ['rs112445263\nC>A\nnew anchor', 'rs313839\nC>G\nolder anchor', 'rs62134782\nC>G\ncomparison'])
    ax.set_ylim(-.027, .039)
    ax.set_ylabel('Predicted gene RNA effect  |  ln(ALT + ε) − ln(REF + ε)')
    ax.set_title('A  Model RNA predictions', loc='left', pad=14)
    ax.text(0, -.31, 'Gray: all 50 fixed genes. Red: PRKD2.\nOne merged CD14-monocyte RNA track.\nNo significance test or causal probability.', transform=ax.transAxes, fontsize=9, va='top', color='#45525a')
    ax = fig.add_subplot(grid[1])
    exps = sorted(peaks.experiment.unique(), key=lambda e: (training.loc[e, 'group'], e))
    columns = order + ['PRKD2_fixed_QTL_TSS']
    matrix = np.zeros((len(exps), 4))
    for i, exp in enumerate(exps):
        for j, label in enumerate(columns):
            r = peaks[(peaks.experiment == exp) & (peaks.label == label)].iloc[0]
            state = 2 if r.exact_base_peak_count > 0 else 1 if r.within_250bp_peak_count > 0 else 0
            matrix[i, j] = state
    ax.imshow(matrix, cmap=ListedColormap(['#edf0f2', '#f4dca6', '#1d827b']), vmin=0, vmax=2, aspect='auto')
    for i in range(8):
        for j in range(4):
            ax.text(j, i, ['—', 'NEAR', 'YES'][int(matrix[i,j])], ha='center', va='center',
                    color='white' if matrix[i,j] == 2 else '#38434c', fontsize=9, fontweight='bold')
    labels = []
    for exp in exps:
        r = training.loc[exp]
        mark = 'T' if r.exact_experiment_in_all_human_training_rows else 'D' if r.same_donor_training_experiments else '?'
        qc = peaks[peaks.experiment == exp].iloc[0].file_audit_levels
        flag = '*' if isinstance(qc, str) and ('WARNING' in qc or 'NOT_COMPLIANT' in qc) else ''
        labels.append(f"{r['group']} · {exp}{flag} [{mark}]")
    ax.set_yticks(range(8), labels, fontsize=9)
    ax.set_xticks(range(4), ['New\nrs112445263', 'Older\nrs313839', 'Comparison\nrs62134782', 'PRKD2\nfixed TSS'], fontsize=9)
    ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
    ax.set_yticks(np.arange(-.5, 8, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=2)
    ax.tick_params(which='minor', length=0)
    ax.set_title('B  Measured ENCODE peak overlap', loc='left', pad=14)
    ax.text(0, -.26, 'YES: edited base/TSS in a peak. NEAR: only within ±250 bp.\n[T] published training experiment; [D] training donor reuse; [?] unresolved.\n* Portal WARNING/NOT_COMPLIANT. Eight experiments represent three donors.',
            transform=ax.transAxes, fontsize=9, va='top', color='#45525a')
    fig.suptitle('PRKD2: the older candidate has the more coherent regulatory signal', x=.075, y=.96, ha='left', fontsize=17, fontweight='bold')
    fig.text(.075, .89, 'Exploratory follow-up • ALL_FOLDS • fixed common 1 Mb window • no matched comparison for the new anchor', fontsize=11, color='#45525a')
    save(fig, output, 'overview')


def fixed_tracks(output):
    inputs = json.loads((ROOT / 'config/prkd2-mechanism-inputs.json').read_text())
    plan = json.loads((ROOT / 'config/prkd2-mechanism-plan.json').read_text())
    manifest = json.loads((RUN / 'run.json').read_text())
    specifications = [('rna_seq', '-', None, 'CD14 RNA, adult minus strand'), ('dnase', '.', None, 'CD14 DNase'),
                      ('cage', '+', None, 'Generic monocyte CAGE +'), ('cage', '-', None, 'Generic monocyte CAGE −'),
                      ('chip_histone', '.', 'H3K27ac', 'CD14 H3K27ac'), ('chip_histone', '.', 'H3K4me1', 'CD14 H3K4me1'),
                      ('chip_histone', '.', 'H3K4me3', 'CD14 H3K4me3')]
    crop_rows = []
    for variant in [v for v in inputs['variants'] if v['role'] == 'anchor']:
        rsid = variant['rsid']
        entry = next(e for e in manifest['anchor_tracks'] if e['rsid'] == rsid)
        fig, axes = plt.subplots(7, 2, figsize=(12, 14), sharex='col', gridspec_kw={'hspace': .62, 'wspace': .2})
        for i, (mod, strand, mark, label) in enumerate(specifications):
            traces = []
            for allele in ['REF', 'ALT']:
                record = next(r for r in entry['outputs'] if r['modality'] == mod and r['allele'] == allele)
                meta = pd.read_csv(RUN / record['metadata_file']).fillna('')
                selector = meta.strand == strand
                if mod == 'rna_seq': selector &= meta.biosample_life_stage == 'adult'
                if mark: selector &= meta.histone_mark == mark
                indices = np.flatnonzero(selector)
                if len(indices) != 1: raise ValueError('Raw track identity is not unique')
                values = np.load(RUN / record['values_file'])['values'][:, indices[0]]
                traces.append((record, values))
            ymax = 0
            for j, (window, point1) in enumerate([('Variant', variant['position_grch38_1based']), ('Fixed PRKD2 TSS', plan['region']['tss_grch38_1based'])]):
                ax = axes[i,j]
                center0 = point1 - 1
                cropped = []
                for allele, (record, values), color, style in zip(['REF','ALT'], traces, ['#535f69','#068a89'], ['-','--']):
                    resolution = record['resolution']
                    starts = record['start0'] + np.arange(len(values)) * resolution
                    keep = (starts < center0 + 2501) & (starts + resolution > center0 - 2500)
                    x = starts[keep] + resolution/2 - center0
                    y = values[keep]
                    ax.plot(x, y, color=color, ls=style, lw=1.1, label=allele)
                    cropped.append(y)
                    ymax = max(ymax, float(y.max()))
                delta = float(np.max(np.abs(cropped[1] - cropped[0])))
                crop_rows.append({'rsid': rsid, 'modality': mod, 'strand': strand, 'histone_mark': mark or '',
                                  'window': window, 'center_position_GRCh38_1based': point1, 'requested_start0': center0-2500,
                                  'requested_end0': center0+2501, 'resolution_bp': resolution, 'bins': len(cropped[0]),
                                  'max_absolute_ALT_minus_REF_in_crop': delta, 'REF_sum_of_intersecting_bins': float(cropped[0].sum(dtype=np.float64)),
                                  'ALT_sum_of_intersecting_bins': float(cropped[1].sum(dtype=np.float64))})
                ax.axvline(0, c='#c0925e', lw=.8)
                ax.text(.02,.98, f'max |ALT − REF| = {delta:.4g}', transform=ax.transAxes, va='top', fontsize=8, color='#45525a')
                ax.set_xlim(-2500,2500)
                ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,3))
                if i == 0: ax.set_title(f'{window}\nchr19:{point1:,} ±2,500 bp', fontsize=11, pad=12)
            for ax in axes[i,:]: ax.set_ylim(0, max(ymax*1.25, 1e-5))
            axes[i,0].set_ylabel(label + '\nmodel units', fontsize=9)
        for ax in axes[-1,:]: ax.set_xlabel('Genomic offset from stated center (bp; increasing coordinate)')
        axes[0,1].legend(loc='upper right', fontsize=8, frameon=False)
        fig.suptitle(f'{rsid} {variant["reference"]}>{variant["alternate"]}: fixed prediction windows', fontsize=16, fontweight='bold', y=.985)
        fig.subplots_adjust(top=.91, bottom=.145, left=.12, right=.97)
        fig.text(.12,.035,'Reference and alternate predictions use the same 1 Mb sequence context and ALL_FOLDS.\nCAGE has mixed monocyte/progenitor-cell source labels. Histone bins are 128 bp; other tracks are 1 bp.\nCrops include bins intersecting the fixed window; these summaries are not the variant-score masks.',fontsize=9,color='#45525a')
        save(fig, output, rsid + '-fixed-tracks')
    pd.DataFrame(crop_rows).to_csv(output / 'fixed-crop-summary.csv',index=False,lineterminator='\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-directory', type=Path, default=DEFAULT)
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    overview(args.output_directory)
    fixed_tracks(args.output_directory)
    (args.output_directory / 'figure-provenance.json').write_text(json.dumps({
        'script_sha256': sha(Path(__file__)), 'model_run_sha256': sha(RUN/'run.json'),
        'summary_sha256': sha(DATA/'summary.json'), 'crop_rule': 'Prespecified ±2500 bp, inclusive integer bases; all intersecting model bins retained.',
        'outputs': {p.name: sha(p) for p in sorted(args.output_directory.iterdir()) if p.name != 'figure-provenance.json'}}, indent=2)+'\n')
    print('Saved overview and both anchor fixed-window figures (PNG/PDF).')
