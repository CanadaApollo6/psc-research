"""Independent expected mathematics; primary helpers are exercised, not reused as oracles."""
from __future__ import annotations
import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts import verify_macromap_program_context as v
from scripts import analyze_macromap_program_context as m


def pair_fixture():
    rows=[]
    for protocol in v.PROTOCOLS:
        for run,n_lines in [('01',1),('02',2),('10',3),('20',1),('30',2)]:
            for number in range(n_lines):
                line=f'{protocol}_{run}_L{number}'
                for time in v.TIMES:
                    for s,stimulus in enumerate(v.STIMULI):
                        if number==1 and s in {1,3,5}:continue
                        base=2+.03*s+.1*number+.02*int(run)+.1*(time==24)
                        gain=.2*number-.04*int(run)+.1*(time==24)+(.3 if protocol==v.PROTOCOLS[0] else -.2)
                        rows.append({'protocol':protocol,'run':run,'line':line,'time':time,'stimulus':stimulus,
                            'status':'available','baseline_loss':base,'extended_loss':base-gain,'gain':gain})
    return rows


def fitted_fixture():
    rng=np.random.default_rng(626217)
    train=[]; test=[]
    for line in range(8):
        for label,stimulus in enumerate(v.STIMULI):
            if line==0 and label%2:continue
            record={'protocol':v.PROTOCOLS[0],'run':f'r{line}','line':f'line{line}',
                'time':6,'stimulus':stimulus,'status':'available'}
            (train if line<6 else test).append(record)
    y=np.asarray([v.STIMULI.index(r['stimulus']) for r in train]); yt=np.asarray([v.STIMULI.index(r['stimulus']) for r in test])
    x=rng.normal(size=(len(train),5)); xt=rng.normal(size=(len(test),5))
    x[:,4]+=y*.7;xt[:,4]+=yt*.7
    result=m.heldout_increment(x,xt,train,test)
    if result['status']!='available':raise AssertionError(result['status'])
    arrays={'train_x':x,'test_x':xt,'train_y':y,'test_y':yt,
        'train_lines':np.asarray([r['line'] for r in train]),'test_lines':np.asarray([r['line'] for r in test]),
        'train_runs':np.asarray([r['run'] for r in train]),'test_runs':np.asarray([r['run'] for r in test]),
        'scaler_mean':result['scaler'].mean,'scaler_scale':result['scaler'].scale,
        'baseline_coefficients':result['baseline'].coefficients,'extended_coefficients':result['extended'].coefficients,
        'baseline_logp':np.asarray([r['baseline_log_probabilities'] for r in result['predictions']]),
        'extended_logp':np.asarray([r['extended_log_probabilities'] for r in result['predictions']])}
    return arrays,result


class TestIndependentIdentityAndReferences(unittest.TestCase):
    def test_whole_run_hash_assignment_including_nomatch_inventory(self):
        samples=[{'line':f'L{r}{j}','run':str(r),'protocol':p} for p in v.PROTOCOLS for r in range(11) for j in range(2)]
        # Source lines must be distinct across protocols in this synthetic panel.
        for row in samples:row['line']=row['protocol']+row['line']
        expected=v.reconstruct_folds(samples)
        self.assertEqual(m.build_folds(samples),expected)
        self.assertEqual(m.build_folds(samples[::-1]),expected)
        bad=samples+[dict(samples[0],run='newrun')]
        with self.assertRaises(v.VerificationError):v.reconstruct_folds(bad)

    def test_reference_uses_training_controls_once_per_line_not_stimuli_or_test(self):
        values=np.arange(12*10,dtype=float).reshape(12,10)
        pairs=[]; folds={}
        for i in range(5):
            folds['P',str(i)]=1 if i==4 else 0
            for s in range(i+1):pairs.append({'protocol':'P','run':str(i),'line':f'L{i}','time':6,'stimulus':f's{s}',
                'status':'available','control_index':i,'treatment_index':i+5})
        # Poison every treatment and the held-out control; matching must not look.
        values[:,4:]=np.nan
        expected,indices=v.training_reference(values,pairs,folds,'P',6,1)
        observed,qc=m.training_control_reference(values,pairs,folds,'P',6,1)
        np.testing.assert_allclose(expected,np.arange(12)*10+1.5)
        np.testing.assert_allclose(observed,expected)
        self.assertEqual(indices,[0,1,2,3]);self.assertEqual(qc['reference_lines'],4)
        altered=copy.deepcopy(pairs);altered[1]['control_index']=8
        with self.assertRaises(v.VerificationError):v.training_reference(values,altered,folds,'P',6,1)

    def test_matching_seed_bins_exclusions_no_replacement_and_linear_score(self):
        ids=[f'ENSG{i:011d}.1' for i in range(160)]
        reference=np.asarray([i//3 for i in range(160)],float)
        target=[8,23,46,70,109,135];excluded=target+[10,12,15,34,80]
        expected=v.independent_matches(reference,ids,target,excluded,program='toy',protocol='P',time=6,fold=2,bins=8,draws=61)
        observed=m.matched_weights(reference,ids,{'toy':target},excluded,'P',6,2,bins=8,replicates=61)['toy']
        self.assertEqual(observed['selection_sha256'],expected['selection_sha256'])
        np.testing.assert_allclose(observed['background_weights'],expected['background_weights'],atol=1e-15)
        for draw in expected['selections']:
            self.assertEqual(len(set(draw)),len(target));self.assertFalse(set(draw)&set(excluded))
        x=np.random.default_rng(33).normal(size=(160,7))
        explicit=np.asarray([x[draw].mean(axis=0) for draw in expected['selections']]).mean(axis=0)
        actual=m.score_with_weights(x,observed,list(range(7)))
        np.testing.assert_allclose(actual['background'],explicit,atol=1e-14)
        np.testing.assert_allclose(actual['adjusted'],x[target].mean(axis=0)-explicit,atol=1e-14)
        impossible=v.independent_matches(reference,ids,target,list(range(160)),program='toy',protocol='P',time=6,fold=2,bins=8,draws=61)
        self.assertFalse(impossible['available'])

    def test_training_scaler_has_equal_line_not_equal_row_weights(self):
        x=np.asarray([[0,1],[10,3],[10,5],[10,7]],float);lines=['A','B','B','B']
        mean,scale,weights=v.weighted_scaler(x,lines)
        expected=m.fit_scaler(x,lines)
        np.testing.assert_allclose(mean,[5,3]);np.testing.assert_allclose(weights,[.5,1/6,1/6,1/6])
        np.testing.assert_allclose(expected.mean,mean);np.testing.assert_allclose(expected.scale,scale)
        self.assertNotEqual(mean[0],x[:,0].mean())
        constant=m.fit_scaler(np.ones((4,1)),lines)
        np.testing.assert_equal(constant.scale,[1]);np.testing.assert_allclose(constant.mean,[1],rtol=0,atol=1e-15)


class TestIndependentSoftmax(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.arrays,cls.primary=fitted_fixture()

    def test_objective_gradient_and_Hessian_against_finite_differences(self):
        rng=np.random.default_rng(31);x=rng.normal(size=(11,2));y=np.arange(11)%3
        lines=['A']*3+['B']*4+['C']*4;weights=v.weighted_scaler(x,lines)[2];coef=rng.normal(size=(3,3))*.2
        obj,grad,hess=v.softmax_calculus(coef,x,y,weights)
        actual_obj,actual_grad=m.ridge_objective(coef.ravel(),x,y,weights,3)
        self.assertAlmostEqual(obj,actual_obj,places=14);np.testing.assert_allclose(grad.ravel(),actual_grad,atol=1e-14)
        epsilon=1e-5; numerical=np.empty(9);numerical_hessian=np.empty((9,9))
        for j in range(9):
            high=coef.copy().ravel();low=high.copy();high[j]+=epsilon;low[j]-=epsilon
            plus=v.softmax_calculus(high.reshape(3,3),x,y,weights);minus=v.softmax_calculus(low.reshape(3,3),x,y,weights)
            numerical[j]=(plus[0]-minus[0])/(2*epsilon)
            numerical_hessian[:,j]=(plus[1].ravel()-minus[1].ravel())/(2*epsilon)
        np.testing.assert_allclose(grad.ravel(),numerical,rtol=1e-7,atol=1e-9)
        np.testing.assert_allclose(hess,numerical_hessian,rtol=1e-7,atol=1e-9)

    def test_intercept_unpenalized_symmetric_slopes_and_mean_not_sum_loss(self):
        rng=np.random.default_rng(9);x=rng.normal(size=(18,2));y=np.arange(18)%9;lines=[str(i//9) for i in range(18)]
        coef=rng.normal(size=(3,9))*.1;w=v.weighted_scaler(x,lines)[2]
        obj,gradient,_=v.softmax_calculus(coef,x,y,w)
        shifted=coef.copy();shifted[0]+=13
        self.assertAlmostEqual(v.softmax_calculus(shifted,x,y,w)[0],obj,places=13)
        repeat=np.repeat(np.arange(18),3)
        repeated=v.softmax_calculus(coef,x[repeat],y[repeat],v.weighted_scaler(x[repeat],np.asarray(lines)[repeat].tolist())[2])[0]
        self.assertAlmostEqual(obj,repeated,places=13)
        perm=np.arange(9)[::-1];inverse=np.argsort(perm)
        self.assertAlmostEqual(v.softmax_calculus(coef[:,perm],x,inverse[y],w)[0],obj,places=13)
        self.assertAlmostEqual(gradient[0].sum(),0,places=14)

    def test_independent_trust_exact_refits_both_models_and_logprobabilities(self):
        checked=v.verify_fit_arrays(self.arrays,refit=True)
        self.assertEqual(checked['status'],'fit_state_verified')
        for state in checked['models'].values():self.assertLess(state['independent_coefficient_max_error'],2e-5)
        mean,scale,_=v.weighted_scaler(self.arrays['train_x'],self.arrays['train_lines'].tolist())
        tr=(self.arrays['train_x']-mean)/scale
        baseline=v.independent_softmax_fit(tr[:,:4],self.arrays['train_y'],self.arrays['train_lines'].tolist())
        np.testing.assert_allclose(baseline['coefficients'],self.primary['baseline'].coefficients,atol=2e-5)

    def test_failed_state_leakage_scaling_and_coefficient_tampering_are_rejected(self):
        altered={key:val.copy() for key,val in self.arrays.items()};altered['test_lines'][0]=altered['train_lines'][0]
        with self.assertRaises(v.VerificationError):v.verify_fit_arrays(altered,refit=False)
        altered={key:val.copy() for key,val in self.arrays.items()};altered['scaler_mean'][0]+=.1
        with self.assertRaises(AssertionError):v.verify_fit_arrays(altered,refit=False)
        altered={key:val.copy() for key,val in self.arrays.items()};altered['extended_coefficients'][1,0]+=.1
        with self.assertRaises(AssertionError):v.verify_fit_arrays(altered,refit=False)
        failed=m.heldout_increment(self.arrays['train_x'],self.arrays['test_x'],
            [{'line':l,'run':r,'time':6,'protocol':'P','stimulus':v.STIMULI[int(y)]} for l,r,y in zip(self.arrays['train_lines'],self.arrays['train_runs'],self.arrays['train_y'])],
            [{'line':l,'run':r,'time':6,'protocol':'P','stimulus':v.STIMULI[int(y)]} for l,r,y in zip(self.arrays['test_lines'],self.arrays['test_runs'],self.arrays['test_y'])],max_iter=0)
        self.assertEqual(failed['status'],'failed_paired_fit');self.assertEqual(failed['predictions'],[])


class TestIndependentRiskAndUncertainty(unittest.TestCase):
    def test_line_time_and_protocol_weights_not_treatment_rows_or_folds(self):
        records=pair_fixture();expected=v.aggregate_losses(records,records)
        actual=m.aggregate_prediction(records,expected_pairs=records)
        self.assertAlmostEqual(actual['primary_gain'],expected['primary_gain'],places=14)
        self.assertEqual(actual['protocol_weights'],expected['protocol_weights'])
        for key,cell in expected['cells'].items():
            for field,value in cell.items():self.assertAlmostEqual(actual['cells'][key][field],value,places=14)
        rawrow=np.mean([r['gain'] for r in records])
        self.assertNotAlmostEqual(rawrow,expected['primary_gain'],places=4)
        # Duplicate stimulus rows cannot silently turn into extra independent units.
        with self.assertRaises(v.VerificationError):v.aggregate_losses(records+[records[0]],records+[records[0]])

    def test_unbalanced_protocol_line_shares_and_once_only_time_average(self):
        records=[r for r in pair_fixture() if r['protocol']==v.PROTOCOLS[0] or r['run'] in {'01','02','10'}]
        actual=m.aggregate_prediction(records,expected_pairs=records)
        expected=v.aggregate_losses(records,records)
        self.assertAlmostEqual(actual['primary_gain'],expected['primary_gain'],places=14)
        self.assertEqual(expected['protocol_weights'],{v.PROTOCOLS[0]:9/15,v.PROTOCOLS[1]:6/15})
        equal_protocol=np.mean([cell['gain'] for cell in expected['cells'].values()])
        self.assertNotAlmostEqual(expected['primary_gain'],equal_protocol,places=4)

    def test_run_cluster_bootstrap_matches_expanded_runs_fixed_protocol_shares(self):
        records=pair_fixture();expected=v.conditional_run_bootstrap(records,records,draws=151)
        actual=m.cluster_bootstrap_gain(records,expected_pairs=records,replicates=151)
        np.testing.assert_allclose(actual['interval'],expected['interval'],rtol=1e-12,atol=1e-12)
        for key,interval in expected['cell_intervals'].items():np.testing.assert_allclose(actual['cell_intervals'][key],interval,rtol=1e-12,atol=1e-12)
        shuffled=list(reversed(records));again=v.conditional_run_bootstrap(shuffled,shuffled,draws=151)
        np.testing.assert_allclose(expected['interval'],again['interval'],atol=1e-12)

    def test_missing_expected_predictions_classes_times_and_nonfinite_never_drop(self):
        records=pair_fixture()
        for bad in [records[1:],records+records[:1]]:
            with self.assertRaises(v.VerificationError):v.aggregate_losses(bad,records)
        bad=[r for r in records if r['time']==6]
        with self.assertRaises(v.VerificationError):v.aggregate_losses(bad,bad)
        bad=[r for r in records if r['stimulus']!='CIL']
        with self.assertRaises(v.VerificationError):v.aggregate_losses(bad,bad)
        bad=copy.deepcopy(records);bad[0]['baseline_loss']=float('nan')
        with self.assertRaises(v.VerificationError):v.aggregate_losses(bad,records)
        self.assertEqual(m.aggregate_prediction(records[1:],expected_pairs=records)['status'],'incomplete_predictions')

    def test_response_bootstrap_keeps_clustered_lines_and_no_Pvalue_screen(self):
        values=np.asarray([-2.,-1,3,5,8,10]);lines=list('abcdef');runs=['r1','r1','r2','r3','r3','r3']
        actual=m.clustered_response_summary(values,lines,runs,replicates=101)
        rng=np.random.Generator(np.random.PCG64(2026091601));medians=[];means=[]
        for _ in range(101):
            sample=np.concatenate([values[np.asarray(runs)==f'r{i+1}'] for i in rng.integers(3,size=3)])
            medians.append(np.median(sample));means.append(np.mean(sample))
        np.testing.assert_allclose(actual['median_interval'],np.quantile(medians,[.025,.975],method='linear'))
        np.testing.assert_allclose(actual['mean_interval'],np.quantile(means,[.025,.975],method='linear'))
        self.assertNotIn('pvalue',actual)
        self.assertEqual(actual['leave_one_run_out_median_range'],[-1.,6.5])


def synthetic_source_design():
    samples=[]
    for p in v.PROTOCOLS:
        for run in range(5):
            for donor in range(2):
                line=f"{p}-r{run}-line{donor}"
                mapped=not(p==v.PROTOCOLS[1] and run==0 and donor==0)
                for time in v.TIMES:
                    for stimulus in v.STIMULI+("Ctrl","PIC"):
                        condition=f"{stimulus}_{time}";index=len(samples)
                        samples.append({"line":line,"run":str(run),"protocol":p,"condition":condition,
                            "sample_id":f"{line}-{condition}","sample_index":index,"identity_mapped":mapped,
                            "HipsciID":line if mapped else "NOMATCH"})
    folds=v.reconstruct_folds(samples)
    # A same-time reference line need not belong to the both-time classifier cohort.
    for p in v.PROTOCOLS:
        run=next(r for q,r in folds if q==p and folds[q,r]!=0)
        line=f"{p}-reference-only"
        for stimulus in ["Ctrl","CIL"]:
            samples.append({"line":line,"run":run,"protocol":p,"condition":f"{stimulus}_6",
                "sample_id":f"{line}-{stimulus}_6","sample_index":len(samples),"identity_mapped":True,"HipsciID":line})
    programs={p:list(range(j*10,(j+1)*10)) for j,p in enumerate(v.ORIGINAL_PROGRAMS)}
    programs[v.NONOVERLAP]=programs['ets2_g1_dn']
    return {"samples":samples,"sample_ids":[r['sample_id'] for r in samples],
        "feature_ids":[f'ENSG{i:011d}.1' for i in range(1000)],"folds":folds,
        "program_indices":programs,"mapping_gates":{p:True for p in programs},"pool_exclusion_indices":list(range(80))}


def write_synthetic_counts(path,counts,design):
    import gzip
    with gzip.open(path,'wt') as stream:
        stream.write('\t'.join(['Chr','Geneid','Start','End','Strand','Length']+design['sample_ids'])+'\n')
        for i,row in enumerate(counts):
            stream.write('\t'.join(['1',design['feature_ids'][i],'1','10','+','10']+[str(int(value)) for value in row])+'\n')


class TestBoundedSourceToScore(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from unittest.mock import patch
        from contextlib import ExitStack
        cls.temp=tempfile.TemporaryDirectory();cls.directory=Path(cls.temp.name)
        cls.design=synthetic_source_design()
        cls.counts=np.random.default_rng(1762).integers(0,9,size=(1000,len(cls.design['samples'])))
        cls.source=cls.directory/'synthetic.gz';write_synthetic_counts(cls.source,cls.counts,cls.design)
        cls.cache=cls.directory/'normalized.f64'
        m.normalize_count_cache(cls.source,cls.cache,cls.design['feature_ids'],cls.design['sample_ids'])
        cls.normalized=np.memmap(cls.cache,dtype='<f8',mode='r',shape=cls.counts.shape)
        cls.output=cls.directory/'engine';cls.output.mkdir()
        functions={'matched_weights':(m.matched_weights,7),'clustered_response_summary':(m.clustered_response_summary,11),
            'cluster_bootstrap_gain':(m.cluster_bootstrap_gain,11)}
        with ExitStack() as stack:
            for name,(subject,draws) in functions.items():
                stack.enter_context(patch.object(m,name,side_effect=lambda *a,_subject=subject,_draws=draws,**k:_subject(*a,**{**k,'replicates':_draws})))
            cls.result=m.execute_fixed_context(cls.normalized,cls.design,cls.output)

    @classmethod
    def tearDownClass(cls):
        cls.normalized._mmap.close();cls.temp.cleanup()

    def test_raw_all_gene_totals_and_normalized_cells_for_hash_selected_samples(self):
        result=v.verify_normalized_subset(self.source,self.normalized,self.design['feature_ids'],self.design['sample_ids'])
        self.assertEqual(result['sample_count'],8);self.assertEqual(result['gene_count'],1000)
        self.assertLess(result['max_absolute_error'],1e-12)
        self.assertEqual(v.deterministic_subset(self.design['sample_ids'],8),v.deterministic_subset(self.design['sample_ids'][::-1],8))
        corrupted=np.asarray(self.normalized).copy();chosen=v.deterministic_subset(self.design['sample_ids'],8)[0]
        corrupted[3,self.design['sample_ids'].index(chosen)]+=.01
        with self.assertRaises(AssertionError):v.verify_normalized_subset(self.source,corrupted,self.design['feature_ids'],self.design['sample_ids'])

    def test_independent_outer_reference_draws_sample_scores_and_training_test_X(self):
        for scope in ['mapped_primary','all_source_lines_sensitivity']:
            pairs=v.reconstruct_pairs(self.design['samples'],include_nomatch=scope!='mapped_primary')
            predictive=v.reconstruct_prediction_cohort(pairs)
            self.assertFalse(any('reference-only' in r['line'] for r in predictive))
            for protocol in v.PROTOCOLS:
                for time in v.TIMES:
                    directory=self.output/scope/f'{protocol}-t{time}-fold0'
                    with np.load(directory/'score-reference-arrays.npz',allow_pickle=False) as state:arrays=dict(state)
                    audit=json.loads((directory/'score-reference-audit.json').read_text())
                    checked,paired=v.verify_context_arrays(arrays,audit,self.normalized,self.design,scope=scope,
                        protocol=protocol,time=time,fold=0,draws=7)
                    self.assertEqual(checked['programs'],9);self.assertEqual(checked['failed_programs'],[])
                    if time==6:self.assertTrue(any('reference-only' in value for value in arrays['training_control_ids']))
                    for variant in ['ets2_g1_dn',v.NONOVERLAP]:
                        with np.load(directory/variant/'fit-arrays.npz',allow_pickle=False) as state:fits=dict(state)
                        fit_audit=json.loads((directory/variant/'fit-audit.json').read_text())
                        selected=[r for r in predictive if r['protocol']==protocol and r['time']==time]
                        train=[r for r in selected if self.design['folds'][protocol,r['run']]!=0]
                        test=[r for r in selected if self.design['folds'][protocol,r['run']]==0]
                        result,records=v.verify_fit_bundle(fits,fit_audit,train,test,variant=variant,paired_subset=paired)
                        self.assertEqual(result['status'],'fit_state_verified');self.assertEqual(len(records),len(test))
                    # Same mathematical audit must reject a changed training control or weight.
                    bad={key:value.copy() for key,value in arrays.items()};bad['reference'][0]+=.01
                    with self.assertRaises(AssertionError):v.verify_context_arrays(bad,audit,self.normalized,self.design,
                        scope=scope,protocol=protocol,time=time,fold=0,draws=7)

    def test_complete_fit_manifest_and_all_expected_prediction_keys(self):
        manifest=json.loads((self.output/'fit-manifest.json').read_text())
        self.assertEqual(len(manifest),80)
        keys=[(r['identity_scope'],r['protocol'],r['time'],r['fold'],r['variant']) for r in manifest]
        self.assertEqual(len(set(keys)),80)
        for scope in ['mapped_primary','all_source_lines_sensitivity']:
            expected=v.reconstruct_prediction_cohort(v.reconstruct_pairs(self.design['samples'],include_nomatch=scope!='mapped_primary'))
            for variant in ['ets2_g1_dn',v.NONOVERLAP]:
                rows=v.read_rows(self.output/scope/variant/'heldout-loss-records.csv')
                for row in rows:
                    for key in ['baseline_loss','extended_loss','gain']:row[key]=float(row[key])
                own=v.conditional_run_bootstrap(rows,expected,draws=11)
                saved=json.loads((self.output/scope/variant/'incremental-endpoint.json').read_text())
                self.assertAlmostEqual(saved['primary_gain'],own['point']['primary_gain'],places=12)
                np.testing.assert_allclose(saved['interval'],own['interval'],atol=1e-12)

    def test_final_pooled_response_means_and_medians_are_descriptive_only(self):
        rows=v.read_rows(self.output/'mapped_primary/line-responses.csv')
        chosen=[r for r in rows if r['stimulus']=='CIL' and r['time']=='6' and r['program_id']=='ets2_g1_dn']
        values=[float(r['adjusted_change']) for r in chosen];lines=[r['line'] for r in chosen]
        runs=[r['run'] for r in chosen];protocols=[r['protocol'] for r in chosen]
        actual=m.pooled_response_summary(values,lines,runs,protocols)
        self.assertAlmostEqual(actual['mean'],np.mean(values));self.assertEqual(actual['median'],np.median(values))
        self.assertIsNone(actual['median_interval']);self.assertIsNone(actual['mean_interval'])
        self.assertEqual(actual['interval_status'],'not_computed_prespecified')
        shares={p:protocols.count(p)/len(protocols) for p in v.PROTOCOLS}
        self.assertEqual(actual['original_contrast_protocol_line_shares'],shares)
        self.assertNotEqual(shares,self.result['mapped_primary']['ets2_g1_dn']['protocol_weights'])

    def test_readonly_complete_engine_verifier_and_missing_expected_cell(self):
        result=v.verify_engine_bundle(self.output,self.normalized,self.design,matching_draws=7,bootstrap_draws=11)
        self.assertEqual(result['fit_cells'],80);self.assertEqual(len(result['source_score_contexts']),8)
        self.assertEqual(len(result['endpoints']),4)
        manifest=self.output/'fit-manifest.json';original=manifest.read_text()
        try:
            manifest.write_text(json.dumps(json.loads(original)[1:]))
            with self.assertRaises(v.VerificationError):v.verify_engine_bundle(self.output,self.normalized,self.design,matching_draws=7,bootstrap_draws=11,refit=False)
        finally:manifest.write_text(original)



class TestReadOnlyVerifierGuards(unittest.TestCase):
    def test_unauthorized_numerical_CLI_rejects_before_metadata_or_arrays(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory);plan=directory/'draft.json'
            plan.write_text(json.dumps({'status':'draft','execute_real_expression':False,'frozen_at_utc':'2026-01-01T00:00:00Z'}))
            argv=['verify','--preparation',str(directory/'prep'),'--preparation-sha256','a'*64,
                '--execution-directory',str(directory/'uncomputed'),'--plan',str(plan),'--plan-sha256',v.digest(plan),
                '--output',str(directory/'no-output.json')]
            with patch.object(sys,'argv',argv),patch.object(v,'verify_preparation') as metadata,patch.object(v.np,'load') as arrays:
                with self.assertRaises(v.VerificationError):v.main()
                metadata.assert_not_called();arrays.assert_not_called()
            self.assertFalse((directory/'no-output.json').exists())

    def test_verification_report_cannot_modify_immutable_input_directory(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory);prep=directory/'prep';prep.mkdir()
            argv=['verify','--preparation',str(prep),'--preparation-sha256','a'*64,'--output',str(prep/'new.json')]
            with patch.object(sys,'argv',argv),patch.object(v,'verify_preparation') as metadata:
                with self.assertRaises(v.VerificationError):v.main()
                metadata.assert_not_called()
            self.assertFalse((prep/'new.json').exists())

    def test_readonly_bundle_paths_never_follow_alias_or_parent_escape(self):
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory);inside=directory/'inside';inside.mkdir()
            target=directory/'target';target.write_text('unchanged')
            (inside/'alias').symlink_to(target)
            for relative in ['../target','alias','/absolute']:
                with self.subTest(relative=relative),self.assertRaises(v.VerificationError):v.contained_artifact(inside,relative)
            self.assertEqual(target.read_text(),'unchanged')

    def test_pins_and_literal_CSV_NA_duplicate_headers_and_ragged_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'file.csv';path.write_text('id,value\nNA,1\n')
            self.assertEqual(v.read_rows(path)[0]['id'],'NA')
            pin=v.digest(path);self.assertEqual(v.pinned(path,pin),path)
            path.write_text('id,id\na,b\n')
            with self.assertRaises(v.VerificationError):v.read_rows(path)
            with self.assertRaises(v.VerificationError):v.pinned(path,pin)
            path.write_text('id,value\na,b,c\n')
            with self.assertRaises(v.VerificationError):v.read_rows(path)

    def test_real_output_gate_requires_frozen_external_plan_pin(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'plan.json'
            for plan in [{'status':'draft','frozen_at_utc':'2026-01-01T00:00:00Z'},
                         {'status':'frozen','frozen_at_utc':'2999-01-01T00:00:00Z'}]:
                path.write_text(json.dumps(plan))
                with self.assertRaises(v.VerificationError):v._frozen_plan(path,v.digest(path))
            path.write_text(json.dumps({'status':'frozen','frozen_at_utc':'2026-01-01T00:00:00Z'}))
            self.assertEqual(v._frozen_plan(path,v.digest(path))['status'],'frozen')
            with self.assertRaises(v.VerificationError):v._frozen_plan(path,'a'*64)


if __name__=='__main__':unittest.main()
