"""Scientific reporting gates must not turn unavailable posteriors into findings."""

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from summarize_coloc import curate_signals


class SignalCoverageTests(unittest.TestCase):
    def inputs(self,coverage=.99):
        row={'gene_name':'TEST','dataset_id':'QTD_TEST','n_approximation':1000,'gwas_component':1,'qtl_component':1}
        overlap=pd.DataFrame([{**row,'common_variants':200,'gwas_signal_mass_retained':1.,'qtl_signal_mass_retained':coverage}])
        raw=pd.DataFrame([{**row,'p12':p,'status':'raw','reason':'',**{f'PP.H{i}.abf':v for i,v in enumerate([0.,0.,0.,.01,.99])}} for p in [1e-6,1e-7,1e-5]])
        return raw,overlap

    def test_high_H4_cannot_override_missing_signal_evidence(self):
        raw,overlap=self.inputs(.51)
        result=curate_signals(raw,overlap)
        self.assertEqual(set(result.status),{'insufficient_signal_overlap'})
        self.assertNotIn('PP.H4.abf',result)

    def test_incomplete_posterior_with_passing_coverage_is_an_error(self):
        raw,overlap=self.inputs()
        raw.loc[0,'PP.H4.abf']=np.nan
        with self.assertRaises(ValueError):curate_signals(raw,overlap)

    def test_filtered_signal_is_restored_as_unavailable(self):
        raw,overlap=self.inputs(.51)
        result=curate_signals(raw.iloc[:0],overlap)
        self.assertEqual(len(result),3)
        self.assertEqual(set(result.status),{'insufficient_signal_overlap'})


if __name__=='__main__':unittest.main()
