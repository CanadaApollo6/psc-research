"""Synthetic aggregate-only tests. No scientific model/real expression imports."""
import copy
import csv
import gzip
import hashlib
import io
import itertools
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import report_macromap_program_context as r


def fixture():
    class MemoryReader:
        def aggregate(self, name): return copy.deepcopy(self.payload[name])
        def aggregate_csv(self, name): return copy.deepcopy(self.csvs[name])
    reader = MemoryReader(); reader.payload = {}; reader.csvs = {}; reader.complete = {"endpoints": {s: {} for s in r.SCOPES}}
    manifest = []
    for scope in r.SCOPES:
        panel = []
        counts = {r.PROTOCOLS[0]: 117, r.PROTOCOLS[1]: 77 if scope == r.SCOPES[0] else 81}
        runs = {r.PROTOCOLS[0]: 34, r.PROTOCOLS[1]: 23}
        for p, t, s, g in itertools.product((*r.PROTOCOLS, "pooled_descriptive"), r.TIMES, r.STIMULI, r.PROGRAMS):
            n = sum(counts.values()) if p == "pooled_descriptive" else counts[p]
            nr = sum(runs.values()) if p == "pooled_descriptive" else runs[p]
            identity = dict(identity_scope=scope, protocol=p, time=t, stimulus=s, program_id=g,
                            expected_paired_lines=0 if s == "PIC" else n, expected_paired_runs=0 if s == "PIC" else nr)
            if s == "PIC": panel.append({**identity, "status": "unavailable_mock_identity"}); continue
            for q in r.QUANTITIES:
                value = (r.PROGRAMS.index(g) % 3 - 1) / 10
                row = {**identity, "quantity": q, "input_statuses": {"available": n}, "status": "available_descriptive",
                       "lines": n, "runs": nr, "median": value, "mean": value, "minimum": value, "maximum": value,
                       "positive": n if value > 0 else 0, "negative": n if value < 0 else 0, "zero": n if value == 0 else 0,
                       "leave_one_line_out_median_range": [value, value], "leave_one_run_out_median_range": [value, value]}
                if p == "pooled_descriptive":
                    row.update(interval_status=r.NO_INTERVAL, median_interval=None, mean_interval=None,
                               median_interval_status=r.NO_INTERVAL, mean_interval_status=r.NO_INTERVAL,
                               original_contrast_protocol_line_shares={k: v / sum(counts.values()) for k, v in counts.items()})
                else:
                    row.update(interval_status=r.RESPONSE_INTERVAL, median_interval=[value-.02, value+.02], mean_interval=[value-.02,value+.02], replicates=5000, seed=2026091601)
                panel.append(row)
        reader.payload[f"{scope}/response-panel.json"] = panel
        folds = []
        for variant in r.VARIANTS:
            cells = {}
            for p, t in itertools.product(r.PROTOCOLS, r.TIMES):
                gain = .1 + .02 * r.VARIANTS.index(variant) + .01 * r.PROTOCOLS.index(p) + .01 * r.TIMES.index(t)
                cells[str((p,t))] = {"baseline_loss": 2., "extended_loss": 2.-gain, "gain": gain}
            ns = r.COUNTS[scope]; weights = {p:ns[p]/sum(ns.values()) for p in r.PROTOCOLS}
            pooled = {k:sum(weights[p]*sum(cells[str((p,t))][k] for t in r.TIMES)/2 for p in r.PROTOCOLS) for k in ("baseline_loss","extended_loss","gain")}
            point = {"status":"available", "primary_gain":pooled["gain"], "pooled":pooled, "cells":cells, "line_counts":ns,
                     "protocol_weights":weights, "interval_status":r.PREDICTION_INTERVAL, "interval":[pooled["gain"]-.02,pooled["gain"]+.02],
                     "cell_intervals":{k:[v["gain"]-.02,v["gain"]+.02] for k,v in cells.items()}, "replicates":5000,"seed":2026091601}
            reader.complete["endpoints"][scope][variant] = copy.deepcopy(point)
            reader.payload[f"{scope}/{variant}/incremental-endpoint.json"] = point
            influence = [{"diagnostic":"leave_one_run_out_fixed_predictions", "protocol":p, "omitted_run":f"SECRET_RUN_{p}_{i}", "status":"available", "gain":pooled["gain"] + (i%3-1)/1000} for p,n in zip(r.PROTOCOLS,(34,22)) for i in range(n)]
            reader.payload[f"{scope}/{variant}/run-influence.json"] = influence
            for p,t,f in itertools.product(r.PROTOCOLS,r.TIMES,range(5)):
                n = ns[p]; testn = n//5 + int(f<n%5)
                allruns=34 if p==r.PROTOCOLS[0] else 22; testr=allruns//5+int(f<allruns%5)
                gate = {"status":"available", "train_pairs":(n-testn)*9,"test_pairs":testn*9,"train_lines":n-testn,"test_lines":testn,
                        "train_runs":allruns-testr,"test_runs":testr,"train_class_lines":{s:n-testn for s in r.STIMULI[:-1]},
                        "train_class_runs":{s:allruns-testr for s in r.STIMULI[:-1]},"missing_test_stimuli":[]}
                folds.append({"identity_scope":scope,"protocol":p,"time":t,"fold":f,"predictor_program":variant,
                              **{k:v for k,v in gate.items() if k!='status'},"fold_gate_status":"available","status":"available","line_weighted_gain":cells[str((p,t))]["gain"]})
                context=f"{scope}/{p}-t{t}-fold{f}"; audit_name=f"{context}/{variant}/fit-audit.json"
                manifest.append({"identity_scope":scope,"variant":variant,"protocol":p,"time":t,"fold":f,"status":"available",
                    "fit_npz":f"{context}/{variant}/fit-arrays.npz","fit_npz_sha256":"a"*64,"fit_audit":audit_name,
                    "stratum_npz":f"{context}/score-reference-arrays.npz","stratum_audit":f"{context}/score-reference-audit.json"})
                reader.payload[audit_name]={"schema":"macromap-private-audit-v1","status":"available","fold_gate":gate,
                    "identity_scope":scope,"variant":variant,"protocol":p,"time":t,"fold":f,
                    "feature_programs":["inflammation","interferon_gamma","oxidative_stress","apoptotic_signaling",variant],
                    "class_order":list(r.STIMULI[:-1]),"penalty":.1,"loss_normalization":"sum_weights_one_equal_lines","intercepts":"unpenalized_symmetric_softmax",
                    "fits":{m:{"status":"available","iterations":5,"gradient_inf":1e-8,"objective":1.8,"optimizer_message":"SECRET_SAMPLE_MUST_NOT_EXPORT"} for m in ('baseline','extended')},
                    "preprocessing_scope":"both_train_and_test_X_use_this_outer_fold_training_control_reference"}
                reader.payload[f"{context}/score-reference-audit.json"]={"scope":scope,"protocol":p,"time":t,"heldout_fold":f,
                    "program_status":{g:"available" for g in r.PROGRAMS},"reference_error":None,"reference":{"reference_lines":n-testn,"reference_runs":allruns-testr,"control_indices":[987654321]},
                    "matching":{g:{"available":True,"reason":"available","seed":42,"target_bin_counts":{"0":4},"pool_bin_counts":{str(i):20 for i in range(20)},"underflow_bins":[],"selection_sha256":"SECRET_SELECTION","baseline_target":.9,"baseline_background":.8} for g in r.PROGRAMS}}
        reader.payload[f"{scope}/fold-diagnostics.json"]=folds
    reader.payload['fit-manifest.json']=manifest
    reader.summary={"protocol_design":{r.PROTOCOLS[0]:{"source_lines":124,"mapped_lines":124,"runs":34,"predictive_lines_both_times":114},r.PROTOCOLS[1]:{"source_lines":85,"mapped_lines":81,"runs":23,"predictive_lines_both_times":71}}}
    reader.csvs['paired-design-aggregate.csv']=[{"protocol":p,"time":str(t),"stimulus":s,"paired_lines":str(0 if s=='PIC' else 117 if p==r.PROTOCOLS[0] else 77),"status":"unavailable_mock_identity" if s=='PIC' else 'source_qualified'} for p,t,s in itertools.product(r.PROTOCOLS,r.TIMES,r.STIMULI)]
    reader.csvs['mapping-inventory.csv']=[]
    for g in r.PROGRAMS:
        derived=g==r.VARIANTS[1]
        row={"dataset_key":"macromap","program_id":g,"source_unique_before_exclusion":"10","excluded_target_gene_count":"0","intended_unique_source_genes":"10","mapped_unique_genes":"8" if derived else "10","missing_source_genes":"0","ambiguous_source_genes":"0","mapping_fraction":"0.8" if derived else "1.0","mapping_gate_pass":"True","role":"SECRET_ROLE_NOT_PUBLIC",
             "derived_mapping_fraction":"0.8" if derived else "","parent_mapping_fraction":"1.0" if derived else "","parent_mapping_gate_pass":"True" if derived else "","parent_program_id":r.VARIANTS[0] if derived else "","parent_mapped_unique_genes":"10.0" if derived else "","comparator_union_mapped_genes":"4.0" if derived else "","parent_comparator_overlap_genes":"2.0" if derived else "","derived_nonoverlap_mapped_genes":"8.0" if derived else ""}
        reader.csvs['mapping-inventory.csv'].append(row)
    return reader


def qc_fixture(memory, plan, plan_sha, complete_sha, payload_pins):
    rows=[]
    for i,entry in enumerate(memory.payload['fit-manifest.json']):
        audit=memory.payload[entry['fit_audit']]
        for model in ('baseline','extended'):
            flag=(len(rows)%4 != 0)
            rows.append({**{k:entry[k] for k in ('identity_scope','variant','protocol','time','fold')},'model':model,
                'coefficient_agreement':True,'coefficient_max_absolute_error':1e-6,'constant_training_feature_count':0,'constant_training_feature_names':[],
                'frozen_numeric_criteria_pass':True,'gradient_inf':1e-9,'iterations':6,'objective':1.8,'optimizer_success':flag,
                'saved_constant_flags_verified':True,'solver_flag_status':'reported_success' if flag else 'reported_unsuccessful',
                'training_feature_count':4 if model=='baseline' else 5,'training_rows':audit['fold_gate']['train_pairs']})
    qc={'created_at_utc':'synthetic','existing_numerical_verification_status':'unchanged','limitations':['SECRET_LIMITATION_NOT_EXPORTED'],
        'metric_definitions':{'unused':'SECRET_DEFINITION_NOT_EXPORTED'},
        'provenance':{'plan_sha256':plan_sha,'primary_completion_sha256':complete_sha,'frozen_verifier_sha256':'c'*64,
            'fit_manifest_sha256':payload_pins['fit-manifest.json'],'diagnostic_wrapper_sha256':'d'*64,'runtime':plan['runtime'],
            'fit_npz_pins':[{**{k:x[k] for k in ('identity_scope','variant','protocol','time','fold')},'sha256':x['fit_npz_sha256']} for x in memory.payload['fit-manifest.json']]},
        'purpose':'post_fit_flag_capture_not_replacement_verification','rows':rows,'schema':'macromap-independent-optimizer-flag-capture-v1',
        'status':'captured_frozen_numeric_criteria_pass','unchanged_criteria':{'coefficient_atol':2e-5,'coefficient_rtol':2e-5,'constant_training_sd_threshold':1e-12,
            'independent_gradient_limit':1e-7,'method':'frozen independent_softmax_fit defaults','raw_solver_success_is_separate_from_numeric_acceptance':True,'scaler_atol':1e-12,'scaler_rtol':1e-10}}
    qc['summary']=qc_summary(rows)
    return qc


def qc_summary(rows):
    names={}
    for row in rows:
        for name in row['constant_training_feature_names']:names[name]=names.get(name,0)+1
    return {'fit_cells':80,'models':160,'optimizer_success_true':sum(x['optimizer_success'] for x in rows),'optimizer_success_false':sum(not x['optimizer_success'] for x in rows),
        'frozen_numeric_criteria_pass':sum(x['frozen_numeric_criteria_pass'] for x in rows),'max_gradient_inf':max(x['gradient_inf'] for x in rows),
        'max_coefficient_absolute_error':max(x['coefficient_max_absolute_error'] for x in rows),'models_with_constant_training_features':sum(x['constant_training_feature_count']>0 for x in rows),
        'constant_training_feature_name_counts':names,'all_saved_constant_flags_verified':True}


def disk_fixture(root):
    memory=fixture(); directory=root/'execution'; prepared=root/'preparation'; directory.mkdir();prepared.mkdir()
    def save(path,obj):
        path.parent.mkdir(parents=True,exist_ok=True);data=(json.dumps(obj,sort_keys=True,allow_nan=False)+'\n').encode();path.write_bytes(data);return r.digest(data)
    payload_pins={name:save(directory/name,obj) for name,obj in memory.payload.items()}
    # Declared NPZ commitments only: there are no model-array files to open.
    payload_pins.update({entry['fit_npz']:entry['fit_npz_sha256'] for entry in memory.payload['fit-manifest.json']})
    prepared_pins={}
    for name,rows in memory.csvs.items():
        text=io.StringIO(newline='');writer=csv.DictWriter(text,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
        data=text.getvalue().encode();(prepared/name).write_bytes(data);prepared_pins[name]=r.digest(data)
    summary={**memory.summary,"artifact_sha256":prepared_pins,"source_pins":{"synthetic_source":{"sha256":"b"*64,"bytes":1}}}
    summary_pin=save(prepared/'source-only-readiness.json',summary)
    contract={"identity_scopes":list(r.SCOPES),"protocols":list(r.PROTOCOLS),"times":list(r.TIMES),"programs":list(r.PROGRAMS),"predictive_variants":list(r.VARIANTS),"folds":5,"pooled_response_mean_interval":False,"pooled_response_median_interval":False,"bootstrap_replicates":5000,"bootstrap_seed":2026091601,"time_weights":[.5,.5]}
    plan={"status":"frozen","execute_real_expression":True,"analysis_id":"macromap-fixed-program-context-v1","method_contract":contract,"output_directory":"execution","prepared_directory":"preparation","prepared_summary_sha256":summary_pin,"preparation_artifacts":prepared_pins,"script_sha256":"a"*64,"runner_script_sha256":"b"*64,"runtime":{"synthetic":True},"code_and_protocol_artifacts":[{"path":"scripts/verify_macromap_program_context.py","sha256":"c"*64,"bytes":1}]}
    plan_path=root/'plan.json';plan_sha=save(plan_path,plan)
    complete={"stage":"execution_complete_all_prespecified_statuses_retained","plan_sha256":plan_sha,**{k:plan[k] for k in ('script_sha256','runner_script_sha256','runtime','prepared_summary_sha256','method_contract')},"artifact_sha256":payload_pins,"source_pins":summary['source_pins'],"endpoints":memory.complete['endpoints']}
    complete_sha=save(directory/'execution-complete.json',complete)
    receipt={"status":r.VERIFIED,"plan_sha256":plan_sha,"execution_receipt_sha256":complete_sha,"verification_is_read_only":True,"effects_read":True,"fit_cells":80,"fit_checks":[{"synthetic":True}]*80,"endpoints":[{}]*4,"responses":[{"identity_scope":s,"complete_cells":1512,"all_points_statuses_checked":True,"pooled_response_intervals":False,"protocol_CIL6h_interval_cells_checked":54} for s in r.SCOPES],"source_score_contexts":[{}]*8,"scope_limits":r.SCOPE_LIMITS,"normalization":{"status":"normalized_source_subset_verified","gene_count":58243,"sample_count":8,"max_absolute_error":1e-15}}
    receipt_path=root/'private-receipt.json';receipt_sha=save(receipt_path,receipt)
    public_path=root/'public-approved-receipt.json';public_sha=save(public_path,receipt)
    qc=qc_fixture(memory,plan,plan_sha,complete_sha,payload_pins);qc_path=root/'optimizer-qc.json';qc_sha=save(qc_path,qc)
    kwargs=dict(root=root,plan_path=plan_path,plan_sha=plan_sha,complete_sha=complete_sha,receipt_path=receipt_path,receipt_sha=receipt_sha,public_path=public_path,public_sha=public_sha,optimizer_path=qc_path,optimizer_sha=qc_sha)
    return kwargs,memory


class AggregateProjectionTests(unittest.TestCase):
    def setUp(self): self.reader=fixture()

    def test_all_keys_counts_labels_and_no_private_fields(self):
        tables=r.assemble(self.reader)
        expected={'responses':3024,'predictions':20,'leave-run-summary':12,'model-folds':80,'model-fits':160,'matching-status':360,'scale-status':80,'mapping':9,'paired-design':80,'cohort-design':4}
        self.assertEqual({k:len(v) for k,v in tables.items()},expected)
        for name,rows in tables.items():
            output=gzip.decompress(r.csv_gzip(rows,r.COLUMNS[name])).decode()
            for secret in ('SECRET_RUN','SECRET_SAMPLE','SECRET_SELECTION','987654321','SECRET_ROLE'):
                self.assertNotIn(secret,output)
        primary=[x for x in tables['predictions'] if x['is_primary_endpoint']]
        self.assertEqual(len(primary),1);self.assertEqual(primary[0]['original_total_lines'],185)
        self.assertEqual({x['original_total_lines'] for x in tables['predictions'] if x['identity_scope']==r.SCOPES[1]},{189})

    def test_missing_duplicate_extra_response_keys_fail(self):
        base=self.reader.payload[f'{r.SCOPES[0]}/response-panel.json']
        for rows in (base[:-1],base+[base[0]],base+[dict(base[0],program_id='foreign')]):
            with self.subTest(size=len(rows)),self.assertRaises(ValueError):r.response_table(rows,r.SCOPES[0])

    def test_pic_never_zero_but_real_zero_stays_zero_in_heatmap(self):
        import numpy as np
        table,_=r.response_table(self.reader.payload[f'{r.SCOPES[0]}/response-panel.json'],r.SCOPES[0])
        matrix=r.response_matrix(table,r.PROTOCOLS[0])
        self.assertEqual(matrix.shape,(9,20));self.assertTrue(np.isnan(matrix[:,[9,19]]).all());self.assertEqual(matrix[1,0],0.)
        bad=copy.deepcopy(self.reader.payload[f'{r.SCOPES[0]}/response-panel.json'])
        next(x for x in bad if x['stimulus']=='PIC')['median']=0.
        with self.assertRaises(ValueError):r.response_table(bad,r.SCOPES[0])

    def test_no_pooled_response_ci_but_predictive_ci_retained(self):
        rows=self.reader.payload[f'{r.SCOPES[0]}/response-panel.json'];row=next(x for x in rows if x['protocol']=='pooled_descriptive' and x['stimulus']!='PIC');row['mean_interval']=[0.,1.]
        with self.assertRaises(ValueError):r.response_table(rows,r.SCOPES[0])
        ep=self.reader.payload[f'{r.SCOPES[0]}/{r.VARIANTS[0]}/incremental-endpoint.json'];data=r.prediction_table(ep,r.SCOPES[0],r.VARIANTS[0])
        self.assertEqual(data[-1]['gain_interval_low'],ep['interval'][0])
        ep.pop('interval')
        with self.assertRaises(ValueError):r.prediction_table(ep,r.SCOPES[0],r.VARIANTS[0])

    def test_response_and_prediction_weights_not_exchangeable(self):
        panel=self.reader.payload[f'{r.SCOPES[0]}/response-panel.json'];pooled=next(x for x in panel if x['protocol']=='pooled_descriptive' and x['stimulus']!='PIC')
        pooled['original_contrast_protocol_line_shares']={p:n/185 for p,n in r.COUNTS[r.SCOPES[0]].items()}
        with self.assertRaises(ValueError):r.response_table(panel,r.SCOPES[0])
        ep=self.reader.payload[f'{r.SCOPES[1]}/{r.VARIANTS[0]}/incremental-endpoint.json'];ep['line_counts']=r.COUNTS[r.SCOPES[0]]
        with self.assertRaises(ValueError):r.prediction_table(ep,r.SCOPES[1],r.VARIANTS[0])

    def test_response_failure_shapes_retained_and_bootstrap_metadata_private(self):
        panel=self.reader.payload[f'{r.SCOPES[0]}/response-panel.json'];original=next(x for x in panel if x.get('quantity')=='adjusted_change' and x['protocol']==r.PROTOCOLS[0])
        minimal={k:original[k] for k in ('identity_scope','protocol','time','stimulus','program_id','expected_paired_lines','expected_paired_runs','quantity','input_statuses')}
        original.clear();original.update(minimal,status='missing_response',median=None)
        table,_=r.response_table(panel,r.SCOPES[0]);failed=[x for x in table if x['status']=='missing_response'];self.assertEqual(len(failed),1);self.assertIsNone(failed[0]['mean'])
        for key in ('seed','replicates','median_interval_status','mean_interval_status'):
            original[key]='SECRET_ID'
            with self.assertRaises(ValueError):r.response_table(panel,r.SCOPES[0])
            original.pop(key)

    def test_failed_and_insufficient_interval_endpoints_retain_five_cells(self):
        for st in ('incomplete_predictions','missing_protocol','unpaired_times'):
            rows=r.prediction_table({'status':st,'primary_gain':None},r.SCOPES[0],r.VARIANTS[0])
            self.assertEqual(len(rows),5);self.assertTrue(all(x['gain'] is None for x in rows))
            with self.assertRaises(ValueError):r.prediction_table({'status':st,'primary_gain':None,'cell_intervals':{"private":[0,1]}},r.SCOPES[0],r.VARIANTS[0])
        ep=self.reader.payload[f'{r.SCOPES[0]}/{r.VARIANTS[0]}/incremental-endpoint.json']
        ep.update(interval_status='insufficient_runs',interval=None)
        for key in ('cell_intervals','replicates','seed'):ep.pop(key)
        rows=r.prediction_table(ep,r.SCOPES[0],r.VARIANTS[0]);self.assertTrue(all(x['gain'] is not None and x['gain_interval_low'] is None for x in rows))

    def test_all_80_fit_keys_required_and_failed_component_not_dropped(self):
        self.reader.payload['fit-manifest.json'].pop()
        with self.assertRaises(ValueError):r.assemble(self.reader)
        self.reader=fixture();name=self.reader.payload['fit-manifest.json'][0]['fit_audit'];self.reader.payload[name]['fits']['extended']['status']='optimizer_failed'
        with self.assertRaises(ValueError):r.assemble(self.reader)

    def test_failed_fit_and_matching_status_preserved_no_scale_invention(self):
        entry=self.reader.payload['fit-manifest.json'][0];s,v,p,t,f=(entry[k] for k in ('identity_scope','variant','protocol','time','fold'))
        entry['status']='missing_predictor';audit=self.reader.payload[entry['fit_audit']];audit['status']='missing_predictor';audit['fits']={m:{'status':'missing_predictor'} for m in ('baseline','extended')}
        fold=next(x for x in self.reader.payload[f'{s}/fold-diagnostics.json'] if (x['predictor_program'],x['protocol'],x['time'],x['fold'])==(v,p,t,f));fold.update(status='missing_predictor',line_weighted_gain=None)
        ep={'status':'incomplete_predictions','primary_gain':None};self.reader.payload[f'{s}/{v}/incremental-endpoint.json']=ep;self.reader.complete['endpoints'][s][v]=ep
        self.reader.payload[f'{s}/{v}/run-influence.json']=[{'status':'incomplete_predictions','diagnostic':'unavailable_incomplete_predictions'}]
        match=self.reader.payload[entry['stratum_audit']];match['program_status'][r.PROGRAMS[0]]='insufficient_control_pool';match['matching'][r.PROGRAMS[0]]={'available':False,'reason':'insufficient_control_pool','seed':42,'target_bin_counts':{'0':9},'pool_bin_counts':{'0':1},'underflow_bins':[0]}
        tables=r.assemble(self.reader)
        self.assertEqual(len(tables['model-fits']),160);self.assertEqual(len(tables['scale-status']),80)
        failed=[x for x in tables['scale-status'] if x['paired_status']=='missing_predictor'];self.assertEqual(failed[0]['scale_status'],'not_recorded_in_aggregate_audit')
        self.assertIn('insufficient_control_pool',{x['matching_reason'] for x in tables['matching-status']})

    def test_leave_run_failed_status_counts_and_identity_removal(self):
        s,v=r.SCOPES[0],r.VARIANTS[0];ep=self.reader.complete['endpoints'][s][v];data=self.reader.payload[f'{s}/{v}/run-influence.json']
        data[0].update(status='missing_protocol',gain=None)
        result=r.influence_table(data,ep,s,v,{r.PROTOCOLS[0]:34,r.PROTOCOLS[1]:22})
        self.assertEqual(result[-1]['unavailable_omissions'],1);self.assertIn('missing_protocol',result[-1]['diagnostic_status_counts'])
        self.assertNotIn('SECRET_RUN',json.dumps(result))
        data[0]['gain']=0.
        with self.assertRaises(ValueError):r.influence_table(data,ep,s,v,{r.PROTOCOLS[0]:34,r.PROTOCOLS[1]:22})

    def test_exact_mapping_counts_blank_null_and_fraction_rejected(self):
        result=r.assemble(self.reader);ordinary=result['mapping'][0];self.assertIsNone(ordinary['parent_mapped_unique_genes'])
        self.assertEqual(result['mapping'][-1]['parent_mapped_unique_genes'],10)
        self.reader.csvs['mapping-inventory.csv'][0]['mapped_unique_genes']='9.9999999999999999'
        with self.assertRaises(ValueError):r.assemble(self.reader)

    def test_unknown_status_bool_count_nan_and_field_privacy(self):
        panel=self.reader.payload[f'{r.SCOPES[0]}/response-panel.json']
        for key,value in (('expected_paired_lines',True),('mean',float('nan')),('status','SECRET_PERSON'),('source_line','SECRET_PERSON')):
            mutated=copy.deepcopy(panel);mutated[0][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):r.response_table(mutated,r.SCOPES[0])
        with self.assertRaises(ValueError):r.csv_gzip([{'subject_id':'SECRET'}],('status',))


class AuthenticatedReportingTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name);self.kwargs,self.memory=disk_fixture(self.root)
    def tearDown(self):self.tmp.cleanup()
    def reader(self):return r.Reader(**self.kwargs)


    def test_optimizer_qc_flags_preserved_separate_from_numeric_pass(self):
        reader=self.reader();tables=r.assemble(reader);rows=tables['independent-optimizer-qc']
        self.assertEqual(len(rows),160);self.assertEqual(sum(not x['optimizer_success'] for x in rows),40)
        self.assertTrue(all(x['frozen_numeric_criteria_pass'] for x in rows))
        self.assertTrue(all(x['constant_training_feature_count']==0 and x['saved_constant_flags_verified'] for x in tables['scale-status']))
        public=json.dumps(reader.public)
        self.assertNotIn('SECRET',public);self.assertNotIn('fit_npz',public)
        self.assertEqual(reader.public['optimizer_diagnostic']['summary']['optimizer_success_false'],40)
        self.assertIn('rerun',rows[0]['diagnostic_scope'])

    def test_optimizer_qc_whitelisted_names_correct_four_five_masks(self):
        reader=self.reader();qc=reader.optimizer_qc
        pair=tuple(qc['rows'][0][k] for k in ('identity_scope','variant','protocol','time','fold'))
        chosen=[x for x in qc['rows'] if tuple(x[k] for k in ('identity_scope','variant','protocol','time','fold'))==pair]
        for row in chosen:
            names=['inflammation']+([pair[1]] if row['model']=='extended' else [])
            row.update(constant_training_feature_names=names,constant_training_feature_count=len(names))
        qc['summary']=qc_summary(qc['rows']);tables=r.assemble(reader)
        scale=next(x for x in tables['scale-status'] if tuple(x[k] for k in ('identity_scope','variant','protocol','time','fold'))==pair)
        self.assertEqual(scale['constant_training_feature_count'],2);self.assertEqual(scale['baseline_constant_training_feature_count'],1)
        baseline=next(x for x in chosen if x['model']=='baseline');baseline['constant_training_feature_names']=[];baseline['constant_training_feature_count']=0
        qc['summary']=qc_summary(qc['rows'])
        with self.assertRaises(ValueError):r.assemble(reader)

    def test_optimizer_qc_missing_tampered_foreign_or_private_shape_rejected(self):
        for field,value in [('constant_training_feature_names',['SECRET_SAMPLE']),('optimizer_success','False'),('solver_flag_status','precision_loss_guess'),('gradient_inf',1e-5),('training_feature_count',5),('training_rows',999),('saved_constant_flags_verified',False)]:
            reader=self.reader();reader.optimizer_qc['rows'][0][field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):r.assemble(reader)
        reader=self.reader();reader.optimizer_qc['rows'].pop()
        with self.assertRaises(ValueError):r.assemble(reader)
        reader=self.reader();reader.optimizer_qc['provenance']['fit_npz_pins'][0]['sha256']='0'*64
        with self.assertRaises(ValueError):r.assemble(reader)
        reader=self.reader();reader.optimizer_qc['provenance']['primary_completion_sha256']='0'*64
        with self.assertRaises(ValueError):r.assemble(reader)
        reader=self.reader();reader.optimizer_qc['rows'][0]['patient_id']='SECRET'
        with self.assertRaises(ValueError):r.assemble(reader)

    def test_optimizer_external_hash_rechecked_and_publication_requires_it(self):
        with self.assertRaises(ValueError):r.Reader(**{**self.kwargs,'optimizer_sha':'0'*64})
        reader=self.reader();path=self.kwargs['optimizer_path'];path.write_bytes(path.read_bytes()+b' ')
        with self.assertRaises(ValueError):reader.recheck()
        # Restore byte pin for a separate deliberately omitted optional API argument.
        no_qc={k:v for k,v in self.kwargs.items() if not k.startswith('optimizer_')}
        with self.assertRaises(ValueError):r.publish(r.Reader(**no_qc),self.root/'no-qc-data',self.root/'no-qc-fig')

    def test_wrong_status_layer_or_uncertainty_metadata_rejected(self):
        with self.assertRaises(ValueError):r.prediction_table({'status':'optimizer_failed','primary_gain':None},r.SCOPES[0],r.VARIANTS[0])
        panel=self.memory.payload[f'{r.SCOPES[0]}/response-panel.json']
        first=panel[0];first.update(interval_status='insufficient_runs',median_interval=None,mean_interval=None)
        with self.assertRaises(ValueError):r.response_table(panel,r.SCOPES[0])
        with self.assertRaises(ValueError):r.Reader(**{**self.kwargs,'public_sha':'x'*64})

    def test_passed_authentication_projection_no_forbidden_inputs(self):
        reader=self.reader();tables=r.assemble(reader);self.assertEqual(len(tables['responses']),3024)
        self.assertTrue(all(not path.suffix=='.npz' and 'log2cpm' not in path.name for path,_ in reader.inputs))
        for name in ('released-log2cpm.f64','mapped_primary/line-responses.csv','mapped_primary/ets2_g1_dn/heldout-loss-records.csv'):
            with self.assertRaises(ValueError):reader.aggregate(name)
        self.assertEqual(reader.public['normalization_max_absolute_error'],1e-15)

    def test_plan_completion_receipt_and_payload_tamper(self):
        for which in ('plan_sha','complete_sha','receipt_sha','public_sha'):
            kwargs={**self.kwargs,which:'0'*64}
            with self.subTest(which=which),self.assertRaises(ValueError):r.Reader(**kwargs)
        reader=self.reader();target=reader.directory/'mapped_primary/response-panel.json';target.write_bytes(target.read_bytes()+b' ')
        with self.assertRaises(ValueError):r.assemble(reader)

    def test_failed_or_foreign_receipt_with_valid_external_hash_rejected(self):
        path=self.kwargs['receipt_path'];obj=json.loads(path.read_text());obj['status']='failed';data=json.dumps(obj).encode();path.write_bytes(data)
        with self.assertRaises(ValueError):r.Reader(**{**self.kwargs,'receipt_sha':r.digest(data)})

    def test_strict_json_duplicate_nonfinite(self):
        for data in ('{"a":1,"a":2}','{"a":NaN}','{"a":Infinity}'):
            with self.assertRaises(ValueError):r.strict_json(data)

    def test_symlink_dangling_and_literal_parent_alias_rejected(self):
        target=self.root/'alias';target.symlink_to(self.root/'plan.json')
        with self.assertRaises(ValueError):r.read_pin(target,self.kwargs['plan_sha'])
        dangling=self.root/'dangling';dangling.symlink_to(self.root/'missing')
        with self.assertRaises(ValueError):r.absolute(dangling)
        with self.assertRaises(ValueError):r.absolute(self.root/'x'/'..'/'new')

    def test_replay_gzip_figures_manifest_are_byte_deterministic(self):
        first=r.publish(self.reader(),self.root/'public-a',self.root/'fig-a')
        second=r.publish(self.reader(),self.root/'public-b',self.root/'fig-b')
        self.assertEqual(first,second)
        for p in (self.root/'public-a').iterdir():self.assertEqual(p.read_bytes(),(self.root/'public-b'/p.name).read_bytes())
        for name in ('prediction','responses'):
            for ext in ('png','svg'):self.assertEqual((self.root/f'fig-a-{name}.{ext}').read_bytes(),(self.root/f'fig-b-{name}.{ext}').read_bytes())
        svg=(self.root/'fig-a-responses.svg').read_text();self.assertIn('N/A',svg);self.assertIn('PIC',svg);self.assertNotIn('SECRET',svg)
        self.assertEqual(json.loads((self.root/'public-a/reporting-manifest.json').read_text())['table_rows']['responses'],3024)

    def test_fresh_outputs_no_overwrite_and_immutable_overlap(self):
        (self.root/'taken').mkdir()
        with self.assertRaises(ValueError):r.publish(self.reader(),self.root/'taken',self.root/'newfig')
        existing=self.root/'reserved-prediction.png';existing.write_bytes(b'keep')
        with self.assertRaises(ValueError):r.publish(self.reader(),self.root/'new',self.root/'reserved')
        self.assertEqual(existing.read_bytes(),b'keep')
        with self.assertRaises(ValueError):r.publish(self.reader(),self.root/'execution/new',self.root/'unused')

    def test_input_changed_during_presentation_blocks_publication(self):
        reader=self.reader();target=reader.directory/'fit-manifest.json'
        def changed(_):target.write_bytes(target.read_bytes()+b' ');return {'prediction_png':b'fake'}
        with patch.object(r,'render_figures',side_effect=changed),self.assertRaises(ValueError):r.publish(reader,self.root/'no-public',self.root/'no-fig')
        self.assertFalse((self.root/'no-public').exists())


if __name__=='__main__':unittest.main()
