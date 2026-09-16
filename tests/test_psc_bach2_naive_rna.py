"""Synthetic method tests. Real-source tests are annotation-only and opt-in CLI.

These tests never execute a real count extraction. Fake matrices use BytesIO or
TemporaryDirectory. No closed count qualification replay is called here.
"""
import copy
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

MODULE = Path(__file__).resolve().parents[1] / 'scripts/analyze_psc_bach2_naive_rna.py'
spec = importlib.util.spec_from_file_location('bach2_naive_method', MODULE)
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


def synthetic_matrix(lines='1 1 2\n2 1 5\n3 1 999\n1 2 3\n2 2 7\n3 2 888\n1 3 1000\n2 3 2000\n', rows=3, cols=3, nnz=None):
    if nnz is None:
        nnz = len(lines.splitlines())
    body = (a.MEX_HEADER + '\n%metadata_json: ' + json.dumps(a.EXPORTER) + f'\n{rows} {cols} {nnz}\n' + lines).encode()
    return body, {'rows': rows, 'columns': cols, 'bytes': len(body), 'sha256': hashlib.sha256(body).hexdigest()}


def donor_selection():
    return {'donors': [{'sample': sample, 'donor': str(i + 1), 'group': 'SNP' if i < 4 else 'noSNP',
        'state_counts': {'C0': 1, 'C1': 1},
        'selected': [{'column_1based': 1, 'barcode': 'AAAA-1', 'state': 'C0'}, {'column_1based': 2, 'barcode': 'CCCC-1', 'state': 'C1'}]}
        for i, sample in enumerate(a.SAMPLES)]}


def fake_plan():
    return {'selection_sha256': '1' * 64, 'source_manifest_sha256': '2' * 64,
            'implementation': [{'sha256': '3' * 64}], 'runtime': {'python': 'synthetic'},
            'public_schema_sha256': a.digest(a.PUBLIC_SCHEMA)}


def auth_for(plan, output):
    return {'schema_version': 1, 'status': a.AUTH_STATUS, 'expression_authorized': True,
        'authorizer': 'root', 'one_fixed_extraction_only': True, 'output_directory': output,
        'plan_sha256': a.digest(plan), 'selection_sha256': plan['selection_sha256'],
        'source_manifest_sha256': plan['source_manifest_sha256'],
        'implementation_sha256': plan['implementation'][0]['sha256'], 'runtime_sha256': a.digest(plan['runtime']),
        'public_schema_sha256': plan['public_schema_sha256'], 'candidate_sha256': a.CANDIDATE_SHA}


def meta_row(sample='sample01', barcode='AAAA-1', state='C0', group='SNP'):
    label, seurat = ('C0: CD4+ TN RTE', '0') if state == 'C0' else ('C1: CD4+ TN mature', '1')
    return {'sample_name': sample, 'barcode': sample + '_' + barcode, 'condition': group, 'sex': 'male',
            'paper_clusters': label, 'paper_clusters_short': state, 'seurat_clusters': seurat}


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.donors = {'sample01': {'donor': '1', 'group': 'SNP'}, 'sample02': {'donor': '2', 'group': 'noSNP'}}
        self.axes = {'sample01': ['AAAA-1', 'CCCC-1', 'GGGG-1'], 'sample02': ['AAAA-1', 'CCCC-1']}
        self.rows = [meta_row(), meta_row('sample01', 'CCCC-1', 'C1'), meta_row('sample02', group='noSNP')]

    def test_sample_key_preserves_cross_sample_barcode_collision(self):
        result = a.select_metadata(self.rows, self.donors, self.axes)
        self.assertEqual(len(result['retained_keys']), 3)
        self.assertEqual([d['selected'][0]['column_1based'] for d in result['donors']], [1, 1])

    def test_no_source_extra_cells_added(self):
        result = a.select_metadata(self.rows, self.donors, self.axes)
        self.assertNotIn(['sample01', 'GGGG-1'], result['retained_keys'])

    def test_wrong_sample_prefix(self):
        self.rows[0]['barcode'] = 'sample02_AAAA-1'
        with self.assertRaises(a.GateError): a.select_metadata(self.rows, self.donors, self.axes)

    def test_prefix_not_arbitrary_split(self):
        self.rows[0]['barcode'] = 'other_sample01_AAAA-1'
        with self.assertRaises(a.GateError): a.select_metadata(self.rows, self.donors, self.axes)

    def test_missing_retained_cell_invalidates_selection(self):
        self.axes['sample01'] = ['CCCC-1']
        with self.assertRaises(a.GateError): a.select_metadata(self.rows, self.donors, self.axes)

    def test_duplicate_retained_key(self):
        self.rows.append(self.rows[0])
        with self.assertRaises(a.GateError): a.select_metadata(self.rows, self.donors, self.axes)

    def test_crosswalk_change(self):
        self.rows[0]['seurat_clusters'] = '1'
        with self.assertRaises(a.GateError): a.select_metadata(self.rows, self.donors, self.axes)

    def test_disguised_C0_label_rejected(self):
        self.rows[0]['paper_clusters'] = 'unapproved naive state'
        with self.assertRaises(a.GateError): a.select_metadata(self.rows, self.donors, self.axes)

    def test_wrong_group_rejected(self):
        self.rows[0]['condition'] = 'noSNP'
        with self.assertRaises(a.GateError): a.select_metadata(self.rows, self.donors, self.axes)

    def test_no_other_state_selected(self):
        other = meta_row('sample01', 'GGGG-1')
        other.update(paper_clusters='C2: TREG', paper_clusters_short='C2', seurat_clusters='2')
        result = a.select_metadata(self.rows + [other], self.donors, self.axes)
        self.assertEqual(result['donors'][0]['retained_count'], 3)
        self.assertEqual(len(result['donors'][0]['selected']), 2)

    def test_whole_donor_selection_required(self):
        with self.assertRaises(a.GateError): a.select_metadata(self.rows[:2], self.donors, self.axes)

    def test_source_barcode_duplicates(self):
        with self.assertRaises(a.GateError): a.parse_barcodes('AAAA-1\nAAAA-1\n')

    def test_full_axis_feature_type_and_one_based_position(self):
        text = '\t'.join(a.TARGET) + '\nENSG_OTHER\tOTHER\tGene Expression\nADT_X\tX\tAntibody Capture\n'
        parsed = a.parse_features(text, target_row=1)
        self.assertEqual(parsed[0], a.TARGET)
        with self.assertRaises(a.GateError): a.parse_features(text, target_row=0)

    def test_duplicate_ID_and_false_modality_rejected(self):
        for text in ('\t'.join(a.TARGET) + '\n' + '\t'.join(a.TARGET), '\t'.join(a.TARGET[:2]) + '\tAntibody Capture'):
            with self.subTest(text=text), self.assertRaises(a.GateError): a.parse_features(text, target_row=1)

    def test_reciprocal_reference_identity(self):
        header = 'Gene stable ID\tNCBI gene (formerly Entrezgene) ID\tHGNC symbol\n'
        self.assertTrue(a.reference_identity(header + 'ENSG00000112182\t60468\tBACH2\n')['annotation_only'])
        with self.assertRaises(a.GateError): a.reference_identity(header + 'ENSG00000112182\t60468\tOTHER\n')


class ExactMatrixTests(unittest.TestCase):
    def parse(self, lines=None, **kwargs):
        body, pin = synthetic_matrix(**({} if lines is None else {'lines': lines}), **kwargs)
        return a.stream_selected_matrix(io.BytesIO(body), pin, [1, 2], [1, 2], target_row=1)

    def test_RNA_only_full_denominator_same_selected_cells(self):
        # ADT999+888 and outside-population RNA1000+2000 must never enter R.
        self.assertEqual(self.parse(), {'B': 5, 'R': 17, 'B_positive_cells': 2, 'selected_cells': 2})

    def test_zero_BACH2_retained_no_detection_gate(self):
        self.assertEqual(self.parse('2 1 5\n3 1 99\n2 2 7\n')['B'], 0)

    def test_exact_integer_above_float_precision(self):
        n = 2 ** 53 + 1
        result = self.parse(f'1 1 {n}\n2 2 2\n')
        self.assertEqual(result['B'], n)
        self.assertEqual(result['R'], n + 2)

    def test_exact_exponent_and_decimal_integer(self):
        self.assertEqual(a.exact_positive_integer('2.0e3'), 2000)
        self.assertEqual(a.exact_positive_integer('+0002'), 2)

    def test_nonintegral_nonfinite_zero_negative_oversized(self):
        for token in ('1.2', 'NaN', 'Infinity', '0', '-1', '1e100', '1e999999999999999999999', '1' * 101, '１'):
            with self.subTest(token=token), self.assertRaises(a.GateError): a.exact_positive_integer(token)

    def test_duplicate_unsorted_and_bounds(self):
        for lines in ('1 1 1\n1 1 2\n', '2 1 2\n1 1 1\n', '4 1 2\n', '1 0 1\n'):
            with self.subTest(lines=lines), self.assertRaises(a.GateError): self.parse(lines)

    def test_nnz_mismatch(self):
        with self.assertRaises(a.GateError): self.parse('1 1 1\n', nnz=2)

    def test_hash_tamper_after_full_scan_fails(self):
        body, pin = synthetic_matrix()
        pin['sha256'] = '0' * 64
        with self.assertRaises(a.GateError): a.stream_selected_matrix(io.BytesIO(body), pin, [1], [1, 2], target_row=1)

    def test_missing_RNA_target(self):
        body, pin = synthetic_matrix()
        with self.assertRaises(a.GateError): a.stream_selected_matrix(io.BytesIO(body), pin, [1], [2], target_row=1)

    def test_duplicate_selected_columns(self):
        body, pin = synthetic_matrix()
        with self.assertRaises(a.GateError): a.stream_selected_matrix(io.BytesIO(body), pin, [1, 1], [1, 2], target_row=1)

    def test_header_metadata_dimensions_and_line_limit(self):
        body, pin = synthetic_matrix()
        for replacement in (body.replace(b'integer general', b'real general'), body.replace(b'cellranger-6.1.1', b'cellranger-3.0.2'), body.replace(b'3 3 8', b'4 3 8'), body + b'1' * 513):
            p = {**pin, 'bytes': len(replacement), 'sha256': hashlib.sha256(replacement).hexdigest()}
            with self.subTest(body=replacement[:50]), self.assertRaises(a.GateError): a.stream_selected_matrix(io.BytesIO(replacement), p, [1], [1, 2], target_row=1)


class NumericalTests(unittest.TestCase):
    def test_zero_and_equal_exact_ratio(self):
        self.assertEqual(a.score_units(0, 3), 0)
        self.assertEqual(a.score_units(1, 10), a.score_units(2, 20))

    def test_integer_units_gate(self):
        for B, R in ((1, 0), (0, 0), (-1, 3), (4, 3), (True, 3), (1.0, 3), (1, 10 ** 109)):
            with self.subTest(B=B, R=R), self.assertRaises(a.GateError): a.score_units(B, R)

    def test_no_spurious_variance_for_constant_point_one(self):
        m, v = a.moment([Decimal('0.1')] * 4)
        self.assertEqual(m, Fraction(1, 10)); self.assertEqual(v, 0)

    def test_both_constant_zero_or_distinct_leave_inference_unavailable(self):
        for other in ('0', '0.1', '3'):
            result = a.welch({'SNP': [Decimal('0.1')] * 4, 'noSNP': [Decimal(other)] * 4})
            self.assertIsNone(result['CI95']); self.assertIsNone(result['P'])
            self.assertEqual(result['SE'], '0')

    def test_one_constant_arm_valid_actual_df_three(self):
        result = a.welch({'SNP': [Decimal('0.1')] * 4, 'noSNP': list(map(Decimal, ('1', '2', '3', '4')))})
        self.assertEqual(Decimal(result['df']), 3)
        self.assertIsNotNone(result['CI95']); self.assertIsNotNone(result['P'])

    def test_equal_variances_actual_df_six(self):
        result = a.welch({'SNP': list(map(Decimal, ('1', '2', '3', '4'))), 'noSNP': list(map(Decimal, ('0', '1', '2', '3')))})
        self.assertEqual(Decimal(result['df']), 6)
        self.assertEqual(Decimal(result['delta']), 1)
        self.assertAlmostEqual(float(result['SE']), (5 / 6) ** .5, places=14)
        self.assertAlmostEqual(float(result['P']), .3153335962012296, places=12)

    def test_unequal_variance_actual_df_not_fixed(self):
        result = a.welch({'SNP': list(map(Decimal, ('1', '2', '3', '4'))), 'noSNP': list(map(Decimal, ('0', '2', '4', '6')))})
        self.assertLess(abs(Fraction(Decimal(result['df'])) - Fraction(75, 17)), Fraction(1, 10 ** 300))

    def test_tiny_variance_below_float_square_retained(self):
        values = [Decimal(i) * Decimal('1e-200') for i in range(4)]
        result = a.welch({'SNP': values, 'noSNP': [Decimal(0)] * 4})
        self.assertGreater(Decimal(result['SE']), 0)
        self.assertEqual(Decimal(result['df']), 3)
        self.assertIsNotNone(result['P'])

    def test_genuine_tiny_variation_around_nonzero_anchor(self):
        with localcontext() as ctx:
            ctx.prec = 320
            values = [Decimal(1) + Decimal(i) * Decimal('1e-250') for i in range(4)]
        result = a.welch({'SNP': values, 'noSNP': [Decimal(1)] * 4})
        self.assertGreater(Decimal(result['SE']), 0)
        self.assertGreater(Decimal(result['delta']), 0)
        self.assertEqual(Decimal(result['df']), 3)

    def test_tiny_realizable_integer_ratio_difference_not_collapsed(self):
        n = 10 ** 99
        self.assertNotEqual(a.score_units(n - 1, n), a.score_units(n, n + 1))

    def test_tail_underflow_not_forced_P_zero(self):
        result = a.welch({'SNP': list(map(Decimal, ('1e-200', '2e-200', '3e-200', '4e-200'))), 'noSNP': [Decimal(1)] * 4})
        self.assertIsNone(result['P']); self.assertIsNotNone(result['CI95'])
        self.assertEqual(result['inference_status'], 'CI_available_P_unavailable_numeric_tail')

    def test_valid_delta_zero_P_one_is_not_degenerate_rescue(self):
        values = list(map(Decimal, ('0', '1', '2', '3')))
        result = a.welch({'SNP': values, 'noSNP': values})
        self.assertEqual(Decimal(result['P']), 1)
        self.assertGreater(Decimal(result['SE']), 0)

    def test_nonfinite_quantile_keeps_point_only(self):
        values = list(map(Decimal, ('0', '1', '2', '3')))
        with patch('scipy.stats.t.ppf', return_value=float('nan')):
            result = a.welch({'SNP': values, 'noSNP': values})
        self.assertEqual(result['delta'], '0'); self.assertIsNone(result['CI95']); self.assertIsNone(result['P'])

    def test_seven_donors_and_nonfinite_rejected(self):
        with self.assertRaises(a.GateError): a.welch({'SNP': [Decimal(0)] * 3, 'noSNP': [Decimal(0)] * 4})
        with self.assertRaises(a.GateError): a.welch({'SNP': [Decimal('NaN')] * 4, 'noSNP': [Decimal(0)] * 4})

    def test_all_eight_leaveout_points_and_sign_counts(self):
        records = [{'sample': str(i), 'group': 'SNP' if i < 4 else 'noSNP', 'Y': str(x)} for i, x in enumerate((0, 0, 0, 8, 1, 1, 1, 1))]
        private, public = a.leaveouts(records)
        self.assertEqual(len(private), 8); self.assertEqual(public['count'], 8)
        self.assertEqual(Decimal(public['minimum']), -1)
        self.assertEqual(public['strict_sign_reversal_count'], 1)
        self.assertEqual(public['positive_count'] + public['negative_count'] + public['zero_count'], 8)
        self.assertFalse(any('P' in row or 'CI' in row for row in private))

    def test_zero_full_point_sign_is_explicit(self):
        records = [{'sample': str(i), 'group': 'SNP' if i < 4 else 'noSNP', 'Y': '0.1'} for i in range(8)]
        _, public = a.leaveouts(records)
        self.assertEqual(public['zero_count'], 8); self.assertEqual(public['full_point_sign'], 0)
        self.assertEqual(public['strict_sign_reversal_count'], 0)

    def test_analyze_whole_eight_and_zero_gene_retained(self):
        measurements = {s: {'B': 0, 'R': 10, 'selected_cells': 2, 'B_positive_cells': 0} for s in a.SAMPLES}
        private = a.analyze_donors(measurements, donor_selection())
        self.assertEqual(private['zero_B_donors'], 8); self.assertEqual(len(private['donors']), 8)
        del measurements[a.SAMPLES[-1]]
        with self.assertRaises(a.GateError): a.analyze_donors(measurements, donor_selection())


    def test_product_equal_group_logs_give_exact_zero_not_roundoff_sign(self):
        Qs = (2, 3, 5, 7, 1, 1, 1, 210)
        measurements = {s: {'B': q - 1, 'R': 1_000_000, 'selected_cells': 2, 'B_positive_cells': int(q > 1)} for s, q in zip(a.SAMPLES, Qs)}
        private = a.analyze_donors(measurements, donor_selection())
        self.assertEqual(private['inference']['delta'], '0')
        self.assertEqual(private['inference']['groups']['SNP']['mean'], private['inference']['groups']['noSNP']['mean'])
        self.assertEqual(private['LOO_summary']['full_point_sign'], 0)
        self.assertEqual(private['LOO_summary']['strict_sign_reversal_count'], 0)

    def test_stable_exact_log_of_fraction_very_near_one(self):
        epsilon = Fraction(1, 10 ** 1000)
        positive = a.log_positive_fraction(1 + epsilon)
        negative = a.log_positive_fraction(1 - epsilon)
        self.assertGreater(positive, 0); self.assertLess(negative, 0)
        self.assertEqual(a.log_positive_fraction(Fraction(1)), 0)
        self.assertLess(abs(Fraction(positive) / epsilon - 1), Fraction(1, 10 ** 300))

    def test_product_contrast_tiny_sign_and_all_eight_source_LOOs(self):
        epsilon = Fraction(1, 10 ** 99)
        args = {'SNP': [Fraction(2), Fraction(3), Fraction(5), Fraction(7) + epsilon],
                'noSNP': [Fraction(1), Fraction(1), Fraction(1), Fraction(210)]}
        _, delta = a.logarithmic_means_and_contrast(args)
        self.assertGreater(delta, 0)
        swapped = {'SNP': args['noSNP'], 'noSNP': args['SNP']}
        _, opposite = a.logarithmic_means_and_contrast(swapped)
        self.assertLess(opposite, 0)
        self.assertLess(abs(Fraction(delta) + Fraction(opposite)), Fraction(1, 10 ** 300))


    def test_fourth_order_product_cancellation_preserves_tiny_true_contrast(self):
        offsets = (1, -1, 8, -8, 4, -4, 7, -7)
        m = {s: {'B': 10 ** 93 + offset, 'R': 10 ** 99, 'selected_cells': 2, 'B_positive_cells': 1}
             for s, offset in zip(a.SAMPLES, offsets)}
        private = a.analyze_donors(m, donor_selection())
        delta = Decimal(private['inference']['delta'])
        self.assertLess(delta, 0)
        self.assertGreater(delta, Decimal('-2e-371'))
        self.assertLess(delta, Decimal('-1e-371'))
        self.assertGreater(Decimal(private['inference']['SE']), 0)
        self.assertEqual(private['LOO_summary']['full_point_sign'], -1)

    def test_decimal_global_precision_and_rounding_cannot_change_method(self):
        from decimal import ROUND_UP
        expected = a.score_units(1, 3)
        with localcontext() as ctx:
            ctx.prec = 9; ctx.rounding = ROUND_UP
            observed = a.score_units(1, 3)
        self.assertEqual(observed, expected)

    def test_any_invalid_RNA_denominator_invalidates_whole_endpoint(self):
        measurements = {s: {'B': 0, 'R': 10, 'selected_cells': 2, 'B_positive_cells': 0} for s in a.SAMPLES}
        measurements[a.SAMPLES[-1]]['R'] = 0
        with self.assertRaises(a.GateError): a.analyze_donors(measurements, donor_selection())


class GuardTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.root_patch = patch.object(a, 'ROOT', self.root)
        self.root_patch.start()
        self.addCleanup(self.root_patch.stop); self.addCleanup(self.tmp.cleanup)
        self.output = a.BASE + '/runs/frozen_synthetic'
        self.plan = fake_plan()
        self.auth = auth_for(self.plan, self.output)

    def check(self, auth=None, plan=None):
        auth = self.auth if auth is None else auth
        return a.check_authorization(auth, a.digest(auth), a.digest(auth), self.plan if plan is None else plan, self.output)

    def test_missing_execution_authorization_does_not_read_any_file(self):
        with patch.object(Path, 'open', side_effect=AssertionError('must not open')):
            with self.assertRaises(a.GateError): a.execute(a.DEFAULT_PLAN, None, None, None)


    def test_unapproved_authorization_file_role_rejected_before_any_open(self):
        bad_paths = ('work/psc-bach2-cohort-qualification/sample01/members/matrix.txt',
            'data/raw/psc-bach2-qualification/cohort/sample01_filtered_feature_bc_matrix.tar.gz',
            'data/derived/psc-bach2-cohort-qualification.json', a.DEFAULT_PLAN + '/plan.json',
            'work/psc-bach2-followup-review/root-method-decision.json', '../authorization.json')
        with patch.object(Path, 'open', side_effect=AssertionError('must not open any body')):
            for path in bad_paths:
                with self.subTest(path=path), self.assertRaises(a.GateError): a.execute(a.DEFAULT_PLAN, path, '0' * 64, self.output)

    def test_root_auth_namespace_size_and_symlink_gate(self):
        relative = a.AUTHORIZATION_ROOT + 'synthetic.json'
        path = self.root / relative
        path.parent.mkdir(parents=True)
        path.write_bytes(b'{}')
        self.assertEqual(a.read_root_authorization(relative), b'{}')
        path.write_bytes(b'x' * (a.MAX_AUTHORIZATION_BYTES + 1))
        with patch.object(Path, 'open', side_effect=AssertionError('no oversized body open')):
            with self.assertRaises(a.GateError): a.read_root_authorization(relative)
        path.unlink()
        (self.root / 'other.json').write_bytes(b'{}')
        path.symlink_to(self.root / 'other.json')
        with self.assertRaises(a.GateError): a.read_root_authorization(relative)
        for bad in (a.AUTHORIZATION_ROOT + 'x/m.json', a.AUTHORIZATION_ROOT + 'matrix.txt', a.AUTHORIZATION_ROOT + '../x.json'):
            with self.subTest(bad=bad), self.assertRaises(a.GateError): a.read_root_authorization(bad)

    def test_preparation_only_decision_cannot_authorize(self):
        bad = {**self.auth, 'expression_authorized': False, 'status': 'APPROVED_METHOD_PREPARATION_ONLY_NOT_EXPRESSION_EXECUTION'}
        with self.assertRaises(a.GateError): self.check(bad)

    def test_auth_SHA_missing_or_tampered(self):
        for expected in (None, '', '0' * 64):
            with self.subTest(expected=expected), self.assertRaises(a.GateError): a.check_authorization(self.auth, expected, a.digest(self.auth), self.plan, self.output)


    def test_auth_boolean_not_replaced_by_integer(self):
        for key in ('expression_authorized', 'one_fixed_extraction_only'):
            bad = {**self.auth, key: 1}
            with self.subTest(key=key), self.assertRaises(a.GateError): self.check(bad)
        with self.assertRaises(a.GateError): self.check({**self.auth, 'schema_version': True})

    def test_every_authorization_field_frozen(self):
        for key in self.auth:
            bad = {**self.auth, key: 'tampered'}
            with self.subTest(key=key), self.assertRaises(a.GateError): self.check(bad)

    def test_code_runtime_plan_source_selection_schema_tamper(self):
        for key in ('implementation', 'runtime', 'source_manifest_sha256', 'selection_sha256', 'public_schema_sha256'):
            changed = copy.deepcopy(self.plan)
            if key == 'implementation': changed[key][0]['sha256'] = '0' * 64
            elif key == 'runtime': changed[key]['python'] = 'changed'
            else: changed[key] = '0' * 64
            with self.subTest(key=key), self.assertRaises(a.GateError): self.check(plan=changed)

    def test_public_schema_allowlist_ignores_private_injected_values(self):
        m = {s: {'B': 0, 'R': 10, 'selected_cells': 2, 'B_positive_cells': 0} for s in a.SAMPLES}
        private = a.analyze_donors(m, donor_selection())
        private['unapproved_gene_endpoint'] = 'SECRET_INDIVIDUAL_VALUE'
        private['inference']['groups']['SNP']['private_join'] = 'SECRET_INDIVIDUAL_VALUE'
        private['LOO_summary']['individual_points'] = 'SECRET_INDIVIDUAL_VALUE'
        public = a.public_projection(private, self.plan, 'a' * 64)
        self.assertNotIn('SECRET_INDIVIDUAL_VALUE', json.dumps(public))
        self.assertEqual(set(public), set(a.PUBLIC_SCHEMA))
        self.assertEqual(set(public['groups']['SNP']), {'n', 'mean', 'sample_sd'})
        self.assertNotIn('donors', public); self.assertNotIn('LOO_private', public)


    def test_nested_canary_under_allowlisted_keys_rejected(self):
        m = {s: {'B': 0, 'R': 10, 'selected_cells': 2, 'B_positive_cells': 0} for s in a.SAMPLES}
        for field in ('delta', 'CI95', 'P', 'SE', 'df'):
            private = a.analyze_donors(m, donor_selection())
            private['inference'][field] = {'sample01': 'CANARY_DONOR_SECRET'}
            with self.subTest(field=field), self.assertRaises(a.GateError): a.public_projection(private, self.plan, 'a' * 64)

    def test_nonfinite_or_array_scalar_and_extra_public_keys_rejected(self):
        m = {s: {'B': 0, 'R': 10, 'selected_cells': 2, 'B_positive_cells': 0} for s in a.SAMPLES}
        original = a.public_projection(a.analyze_donors(m, donor_selection()), self.plan, 'a' * 64)
        for bad in ('NaN', 'Infinity', 'sample01', ['1', '2'], {'sample': '1'}, True):
            changed = copy.deepcopy(original)
            changed['contrast']['SNP_minus_noSNP'] = bad
            with self.subTest(bad=bad), self.assertRaises(a.GateError): a.validate_public_result(changed)
        original['donor_values'] = ['secret']
        with self.assertRaises(a.GateError): a.validate_public_result(original)

    def test_ignored_private_work_policy_is_explicit(self):
        a.validate_private_ignore('work/\n!.env.example\ndata/raw/\n')
        for text in ('data/raw/\n', 'work/\n!work/\n', 'work/\n!*\n', 'work/\n!**/measurements.private.json\n'):
            with self.subTest(text=text), self.assertRaises(a.GateError): a.validate_private_ignore(text)


    def test_installed_runtime_file_tamper_detected_even_if_RECORD_unchanged(self):
        import base64
        file = self.root / 'component.py'
        file.write_bytes(b'fixed runtime')
        encoded = base64.urlsafe_b64encode(hashlib.sha256(file.read_bytes()).digest()).decode().rstrip('=')
        class Distribution:
            version = 'synthetic'
            def read_text(inner, name):
                return 'component.py,sha256=' + encoded + ',13\npackage.dist-info/RECORD,,\n'
            def locate_file(inner, name):
                return self.root / name
        self.assertEqual(a.verify_distribution(Distribution())['verified_installed_files'], 1)
        file.write_bytes(b'changed bytes')
        with self.assertRaises(a.GateError): a.verify_distribution(Distribution())

    def test_public_nonzero_synthetic_inference_passes_recursive_schema(self):
        m = {s: {'B': i + 1, 'R': 1000, 'selected_cells': 2, 'B_positive_cells': 1} for i, s in enumerate(a.SAMPLES)}
        public = a.public_projection(a.analyze_donors(m, donor_selection()), self.plan, 'a' * 64)
        self.assertIsNotNone(public['contrast']['CI95'])
        self.assertGreater(Decimal(public['contrast']['P_two_sided']), 0)

    def test_unsafe_output_paths(self):
        for path in ('/tmp/run', 'data/derived/run', a.BASE + '/runs/../bad', a.BASE + '/runs/x/y', a.BASE + '/runs/.', a.BASE + '/runs/x.csv', a.BASE + '/runs//x'):
            with self.subTest(path=path), self.assertRaises(a.GateError): a.fresh_directory(path)

    def test_fresh_output_exclusive_and_reuse_refused(self):
        self.check()
        path = a.fresh_directory(self.output, create=True)
        self.assertTrue(path.is_dir())
        with self.assertRaises(a.GateError): self.check()

    def test_symlink_parent_rejected(self):
        base = self.root / a.BASE
        base.mkdir(parents=True)
        (base / 'runs').symlink_to(self.root)
        with self.assertRaises(a.GateError): self.check()

    def test_source_path_lexical_escape_or_symlink_rejected(self):
        (self.root / 'link').symlink_to(self.root)
        for path in ('../x', '/tmp/x', 'a/../x', 'a//x', 'a\\x', 'link/x'):
            with self.subTest(path=path), self.assertRaises(a.GateError): a.source_path(path)

    def test_default_validate_only_cannot_execute_or_write_results(self):
        with patch.object(a, 'validate_plan', return_value=(self.plan, {})) as validation, patch.object(a, 'execute', side_effect=AssertionError('no execution')), patch.object(Path, 'write_bytes', side_effect=AssertionError('no writes')), patch.object(Path, 'open', side_effect=AssertionError('no counts')), patch('sys.stdout', new_callable=io.StringIO):
            self.assertEqual(a.main([]), 0)
            validation.assert_called_once()

    def test_validate_failure_still_cannot_execute_or_write(self):
        with patch.object(a, 'validate_plan', side_effect=a.GateError('tampered')), patch.object(a, 'execute', side_effect=AssertionError('no execution')), patch.object(Path, 'write_bytes', side_effect=AssertionError('no writes')), patch('sys.stderr', new_callable=io.StringIO):
            self.assertEqual(a.main([]), 2)

    def test_execute_flags_not_silently_accepted_in_default_mode(self):
        with patch.object(a, 'validate_plan', side_effect=AssertionError('not reached')), patch('sys.stderr', new_callable=io.StringIO):
            self.assertEqual(a.main(['--authorization', 'fake.json']), 2)


    def test_wrong_plan_directory_role_rejected_before_any_open(self):
        paths = ('data/raw/psc-bach2-qualification', 'work/psc-bach2-cohort-qualification/sample01/members',
                 a.BASE + '/runs', a.BASE + '/runs/old-run', a.BASE + '/independent-review', a.BASE + '/../other')
        with patch.object(Path, 'open', side_effect=AssertionError('no wrong-role file body open')):
            for path in paths:
                with self.subTest(path=path), self.assertRaises(a.GateError): a.validate_plan(path)

    def test_source_only_compilation_and_exact_replay(self):
        selection = {'synthetic_annotation': True}
        with patch.object(a, 'build_source_plan', return_value=(self.plan, selection)), patch.object(a, 'stream_selected_matrix', side_effect=AssertionError('no matrices')):
            a.compile_plan(a.DEFAULT_PLAN)
            before = {p.name: p.read_bytes() for p in (self.root / a.DEFAULT_PLAN).iterdir()}
            self.assertEqual(a.validate_plan(a.DEFAULT_PLAN), (self.plan, selection))
            after = {p.name: p.read_bytes() for p in (self.root / a.DEFAULT_PLAN).iterdir()}
            self.assertEqual(before, after)
            with self.assertRaises(a.GateError): a.compile_plan(a.DEFAULT_PLAN)

    def test_selection_or_plan_file_tamper_rejected(self):
        with patch.object(a, 'build_source_plan', return_value=(self.plan, {'annotation': True})):
            a.compile_plan(a.DEFAULT_PLAN)
            selection = self.root / a.DEFAULT_PLAN / 'selection.private.json'
            selection.write_text('{}')
            with self.assertRaises(a.GateError): a.validate_plan(a.DEFAULT_PLAN)



class CompletionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.addCleanup(self.temp.cleanup)
        self.patch = patch.object(a, 'ROOT', self.root)
        self.patch.start(); self.addCleanup(self.patch.stop)
        self.plan = fake_plan()
        self.output = a.BASE + '/runs/completion_fixture'
        self.directory = self.root / self.output
        self.directory.mkdir(parents=True)
        for name in a.RUN_PAYLOADS:
            (self.directory / name).write_bytes(a.canonical({'synthetic': True, 'role': name}))
        self.auth_sha = 'a' * 64

    def complete(self):
        return a.write_completion(self.directory, self.plan, self.auth_sha, self.output)

    def validate(self):
        return a.validate_completed_run(self.directory, self.plan, self.auth_sha, self.output)

    def test_completion_exact_payload_bytes_hashes_inventory_and_provenance(self):
        record = self.complete()
        self.assertEqual(record, self.validate())
        self.assertEqual(record['allowed_inventory'], list(a.RUN_INVENTORY))
        self.assertEqual(record['provenance'], a.result_provenance(self.plan, self.auth_sha))
        for name in a.RUN_PAYLOADS:
            body = (self.directory / name).read_bytes()
            self.assertEqual(record['payloads'][name], {'bytes': len(body), 'sha256': hashlib.sha256(body).hexdigest()})

    def test_missing_completion_is_not_a_completed_result(self):
        with self.assertRaises(a.GateError): self.validate()

    def test_missing_or_modified_payload_invalidates_completion(self):
        self.complete()
        path = self.directory / a.RUN_PAYLOADS[0]
        body = path.read_bytes()
        path.write_bytes(body + b' ')
        with self.assertRaises(a.GateError): self.validate()
        path.unlink()
        with self.assertRaises(a.GateError): self.validate()

    def test_extra_file_or_directory_invalidates_completion(self):
        self.complete()
        extra = self.directory / 'unapproved.json'
        extra.write_text('{}')
        with self.assertRaises(a.GateError): self.validate()
        extra.unlink(); extra.mkdir()
        with self.assertRaises(a.GateError): self.validate()

    def test_modified_completion_status_or_provenance_rejected(self):
        record = self.complete()
        for mutation in ('status', 'binding'):
            changed = copy.deepcopy(record)
            if mutation == 'status': changed['status'] = 'partial'
            else: changed['provenance']['authorization_sha256'] = 'b' * 64
            (self.directory / 'completion.json').write_bytes(a.canonical(changed))
            with self.subTest(mutation=mutation), self.assertRaises(a.GateError): self.validate()

    def test_completion_publish_failure_leaves_no_completion(self):
        with patch.object(Path, 'replace', side_effect=OSError('synthetic publication failure')):
            with self.assertRaises(OSError): self.complete()
        self.assertFalse((self.directory / 'completion.json').exists())
        self.assertFalse((self.directory / 'completion.pending').exists())
        with self.assertRaises(a.GateError): self.validate()

    def test_execute_failure_before_completion_preserves_failure_marker(self):
        # Fully synthetic outputs and extraction adapter. No real source is read.
        output = a.BASE + '/runs/failed_execution'
        plan = fake_plan()
        plan['source_manifest'] = {'matrix_expected_pins_NOT_read': []}
        auth = auth_for(plan, output)
        relative = a.AUTHORIZATION_ROOT + 'synthetic.json'
        path = self.root / relative
        path.parent.mkdir(parents=True)
        path.write_bytes(a.canonical(auth))
        measurements = {s: {'B': 0, 'R': 10, 'selected_cells': 2, 'B_positive_cells': 0} for s in a.SAMPLES}
        private = a.analyze_donors(measurements, donor_selection())
        with patch.object(a, 'validate_plan', return_value=(plan, donor_selection())), patch.object(a, 'analyze_donors', return_value=private), patch.object(a, 'write_completion', side_effect=OSError('synthetic pre-completion failure')):
            with self.assertRaises(OSError): a.execute(a.DEFAULT_PLAN, relative, a.digest(auth), output)
        destination = self.root / output
        self.assertTrue((destination / 'failure.json').exists())
        self.assertFalse((destination / 'completion.json').exists())
        with self.assertRaises(a.GateError): a.validate_completed_run(destination, plan, a.digest(auth), output)


if __name__ == '__main__':
    unittest.main()
