"""Check nominal-P/beta/SE consistency without assuming a hidden sample size."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, t

from analyze_prkd2_recovery import save_csv
from coloc_common import ROOT, save_json, sha


def main():
    plan_path = ROOT / 'config/prkd2-recovery-ibdverse-format-plan.json'
    plan = json.loads(plan_path.read_text())
    rows = []
    for filename, digest in plan['frozen_source_files'].items():
        path = ROOT / filename
        if sha(path.read_bytes()) != digest:
            raise ValueError('Changed source for distribution-convention audit')
        frame = pd.read_csv(path, sep='\t')
        if not (np.isfinite(frame.slope) & np.isfinite(frame.slope_se) & frame.slope_se.gt(0) & frame.pval_nominal.gt(0) & frame.pval_nominal.le(1)).all():
            raise ValueError('Invalid exported association values')
        statistic = np.abs(frame.slope.to_numpy() / frame.slope_se.to_numpy())
        observed = np.log(frame.pval_nominal.to_numpy())
        errors = [np.abs((np.log(2) + t.logsf(statistic, df) - observed) / np.log(10)) for df in range(1, 501)]
        inferred = int(np.argmin([np.median(x) for x in errors])) + 1
        error = errors[inferred - 1]
        normal = np.abs((np.log(2) + norm.logsf(statistic) - observed) / np.log(10))
        rows.append({'dataset_id': path.name.removesuffix('-PRKD2.tsv.gz'), 'source_rows_checked': len(frame),
                     'best_matching_integer_residual_df': inferred, 'median_abs_log10_P_error': float(np.median(error)),
                     'p95_abs_log10_P_error': float(np.quantile(error, .95)), 'max_abs_log10_P_error': float(np.max(error)),
                     'normal_median_abs_log10_P_error': float(np.median(normal)), 'normal_max_abs_log10_P_error': float(np.max(normal)),
                     'very_low_residual_df_flag': inferred < 30, 'actual_sample_N_recovered': False,
                     'qualification': 'distribution-consistency inference only; not observed sample N or covariate count'})
    output = save_csv(pd.DataFrame(rows), ROOT / 'data/derived/prkd2-recovery-ibdverse-format-audit.csv')
    record = {'plan_sha256': sha(plan_path.read_bytes()), 'script_sha256': sha(Path(__file__).read_bytes()),
              'contexts': rows, 'output': output, 'new_H4_posteriors': 0,
              'interpretation': 'Retain all source values and numerical coverage diagnostics. A very-low-df context does not qualify a Gaussian-BF interpretation merely because its numerical overlap exceeds 90%.'}
    save_json(ROOT / 'reports/prkd2-recovery-ibdverse-format-audit.json', record)
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    main()
