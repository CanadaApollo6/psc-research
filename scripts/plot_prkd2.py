"""Static research figures drawn directly from the complete PRKD2 audit tables."""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

from prkd2_common import ROOT, save_json, sha


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False,
                         'axes.spines.right': False, 'svg.hashsalt': 'prkd2-v1'})
    directory = ROOT / 'reports/figures'
    directory.mkdir(exist_ok=True)
    inputs = []
    def read(name):
        path = ROOT/name
        inputs.append(path)
        return pd.read_csv(path)
    baseline = read('reports/prkd2-baseline.csv').query("sd_mode == 'unit' and p12 == 0.000001")
    expanded = read('reports/prkd2-expanded-baseline.csv').query("sd_mode == 'unit' and p12 == 0.000001")
    signals = read('reports/prkd2-expanded-signal-comparison.csv').query('n_approximation == 11386 and p12 == 0.000001')
    fig, axes = plt.subplots(1, 2, figsize=(13.3, 5.8), gridspec_kw={'width_ratios': [1, 1.1]})
    fig.suptitle('PRKD2: incomplete shared evidence prevents a qualified signal comparison', x=.06, ha='left', fontsize=15, weight='bold')
    labels, original, repaired = [], [], []
    for ds, label in [('QTD000021', 'BLUEPRINT'), ('QTD000504', 'DICE')]:
        for trait in ['gwas', 'qtl']:
            labels.append(label+' / '+('PSC' if trait == 'gwas' else 'RNA'))
            original.append(baseline.set_index('dataset_id').loc[ds, 'full_'+trait+'_bf_mass_retained'])
            repaired.append(expanded.set_index('dataset_id').loc[ds, 'full_'+trait+'_bf_mass_retained'])
    y = np.arange(len(labels))
    axes[0].barh(y-.16, original, height=.3, color='#9ba8b1', label='Original SNP-only scope')
    axes[0].barh(y+.16, repaired, height=.3, color='#167d9a', label='Full-allele extension')
    for position, value in zip(y+.16, repaired):
        axes[0].text(value+.012, position, f'{value:.1%}', va='center', fontsize=9)
    axes[0].set_yticks(y, labels); axes[0].invert_yaxis()
    axes[0].set_title('A. Full regional marginal evidence retained', loc='left', pad=15)
    axes[0].legend(loc='lower left', bbox_to_anchor=(0, -.32), frameon=False, fontsize=9)
    signal_labels, gvals, qvals = [], [], []
    for row in signals.itertuples():
        signal_labels.append(('BLUEPRINT' if row.dataset_id == 'QTD000021' else 'DICE')+f' / RNA signal {row.qtl_component}')
        gvals.append(row.gwas_signal_mass_retained); qvals.append(row.qtl_signal_mass_retained)
    axes[1].barh(y-.16, gvals, height=.3, color='#4b6475', label='PSC component')
    axes[1].barh(y+.16, qvals, height=.3, color='#c76731', label='Published RNA component')
    for position, value in zip(y+.16, qvals):
        axes[1].text(value+.012, position, f'{value:.1%}', va='center', fontsize=9)
    axes[1].set_yticks(y, signal_labels); axes[1].invert_yaxis()
    axes[1].set_title('B. Full component weight retained after extension', loc='left', pad=15)
    axes[1].legend(loc='lower left', bbox_to_anchor=(0, -.32), frameon=False, fontsize=9)
    for ax in axes:
        ax.set_xlim(0, 1.04); ax.xaxis.set_major_formatter(PercentFormatter(1))
        ax.axvline(.9, color='#333333', linestyle='--', linewidth=1)
        ax.set_xlabel('Weight represented by shared variants')
        ax.grid(axis='x', alpha=.15); ax.set_axisbelow(True)
    fig.text(.06, .015, 'Dashed line: prespecified 90% minimum in BOTH traits. B uses the primary scalar-N approximation (11,386).\nAll 72 component/prior/N comparisons fail coverage; missing shared-signal probabilities are unavailable, not zero.', fontsize=9)
    fig.subplots_adjust(left=.12, right=.97, bottom=.29, top=.82, wspace=.61)
    outputs=[]
    for extension in ['png', 'svg']:
        path=directory/('prkd2-evidence-coverage.'+extension)
        fig.savefig(path, dpi=180, metadata={'Date': None} if extension == 'svg' else None)
        if extension == 'svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
        outputs.append(path)
    plt.close(fig)

    g=read('data/derived/prkd2-expanded-inputs/GWAS-PRKD2.csv.gz')
    q1=read('data/derived/prkd2-expanded-inputs/QTD000021-PRKD2.csv.gz')
    q2=read('data/derived/prkd2-expanded-inputs/QTD000504-PRKD2.csv.gz')
    fig, axes=plt.subplots(3, 1, figsize=(11.5, 8), sharex=True)
    fig.suptitle('PRKD2 regional association patterns and missing comparisons', x=.10, ha='left', fontsize=15, weight='bold')
    pairs=[(g, set(q1.snp)&set(q2.snp), 'PSC discovery GWAS: 4,570 identity-qualified variants'),
           (q1, set(g.snp), 'BLUEPRINT monocyte RNA: 7,505 distinct variants'),
           (q2, set(g.snp), 'DICE classical-monocyte RNA: 6,261 distinct variants')]
    for ax, (frame, shared, label) in zip(axes, pairs):
        mask=frame.snp.isin(shared)
        ax.scatter(frame.loc[mask,'position']/1e6, -np.log10(frame.loc[mask,'pvalue']), s=8, alpha=.45, color='#167d9a', rasterized=True)
        ax.scatter(frame.loc[~mask,'position']/1e6, -np.log10(frame.loc[~mask,'pvalue']), s=13, alpha=.8, marker='x', color='#c76731', rasterized=True)
        ax.axvline(46.717127, color='#555555', linestyle=':', linewidth=1)
        ax.set_title(label, loc='left', fontsize=11); ax.set_ylabel('−log10(P)')
        ax.grid(alpha=.12); ax.set_axisbelow(True)
    missing=g[g.source_id=='rs112445263'].iloc[0]
    axes[0].annotate('rs112445263: absent from both RNA files', xy=(missing.position/1e6, -np.log10(missing.pvalue)),
                     xytext=(46.05, 7.9), arrowprops={'arrowstyle':'->','color':'#333333'}, fontsize=9)
    axes[-1].set_xlabel('Chromosome 19 position (GRCh38, Mb)')
    axes[-1].set_xlim(45.717127,47.717127)
    fig.text(.10,.035,'Blue dots: shared variants; orange crosses: unavailable in a comparison. PSC sharing requires both RNA datasets.\nThe full frozen ±1-Mb window is shown. Dotted line: PRKD2 TSS. Plotted PSC records passed full-allele identity checks.',fontsize=9)
    fig.subplots_adjust(left=.10,right=.97,bottom=.14,top=.90,hspace=.38)
    for extension in ['png','svg']:
        path=directory/('prkd2-regional-associations.'+extension)
        fig.savefig(path,dpi=180,metadata={'Date':None} if extension=='svg' else None)
        if extension == 'svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
        outputs.append(path)
    plt.close(fig)
    save_json(ROOT/'reports/prkd2-figure-provenance.json',{'inputs':{str(p.relative_to(ROOT)):sha(p.read_bytes()) for p in inputs},
        'outputs':{str(p.relative_to(ROOT)):sha(p.read_bytes()) for p in outputs},'plot_script_sha256':sha(Path(__file__).read_bytes())})


if __name__=='__main__':
    main()
