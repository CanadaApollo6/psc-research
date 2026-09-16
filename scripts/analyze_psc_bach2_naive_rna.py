#!/usr/bin/env python3
"""One fixed BACH2 RNA endpoint. Default is annotation-only validation.

No real execution is authorized by this file, a compiled plan, or the method-
preparation decision. --execute needs a separately root-frozen authorization
and its externally supplied SHA256. All run artifacts stay under ignored work/.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
import sys
from collections import Counter
from decimal import Context, Decimal, InvalidOperation, localcontext
from fractions import Fraction

ROOT = Path(__file__).resolve().parents[1]
BASE = 'work/psc-bach2-naive-rna-plan'
DEFAULT_PLAN = BASE + '/candidate-v4'
AUTHORIZATION_ROOT = 'work/psc-bach2-naive-rna-freeze/'
MAX_AUTHORIZATION_BYTES = 8192
RUN_PAYLOADS = ('measurements.private.json', 'public-aggregate.json', 'authorization.used.json')
RUN_INVENTORY = tuple(sorted((*RUN_PAYLOADS, 'completion.json')))
COMPLETION_STATUS = 'completed_fixed_endpoint_bound_payload_inventory'
CANDIDATE = 'work/psc-bach2-followup-review/candidate-specification.json'
CANDIDATE_SHA = 'b54c4fc47a93afc48b6463410902e467e459049cef759dc3ba50d8fb6149d627'
REVIEW = 'reports/psc-bach2-followup-design-review.md'
REVIEW_SHA = 'eef554dbbef39aab08feff1a41b9238eff317c10d5cd62a46b46e77b14e5c368'
DECISION = 'work/psc-bach2-followup-review/root-method-decision.json'
METADATA = 'data/raw/psc-bach2-qualification/05-scrna_meta.tsv.body'
META_SHA = '0ba598b3edd30dc7f2ba76f0e7779725e60cce50daa01a0e2e04faebcec7c26c'
SDRF = 'data/raw/psc-bach2-qualification/01-E-MTAB-14013.sdrf.txt.body'
SDRF_SHA = '58273fa29690bdcaff0060161519ff6bbbee9e2af0df4bef03ef91a94184b10d'
REFERENCE = 'work/psc-liver-program-inputs/source-audit/ensembl116-human-gene-crosswalk.reconstructed.tsv'
REFERENCE_SHA = 'a60cf0e0ac5782e7c39269a577d1ee24be08d26c100324fe2c38b77a7f81f167'
FEATURE_SHA = 'fb11f5aa4ca7f366ab3ff0ceb22579a6d9bfd386ce66babc7519c280853b3e1e'
TARGET = ('ENSG00000112182', 'BACH2', 'Gene Expression')
TARGET_ROW = 12314
SAMPLES = tuple(f'sample{i:02d}' for i in range(1, 9))
GROUPS = ('SNP', 'noSNP')
STATES = {'C0: CD4+ TN RTE': ('C0', '0'), 'C1: CD4+ TN mature': ('C1', '1')}
PRECISION = 320
MAX_COORDINATES = 50_000_000
MAX_VALUE_DIGITS = 100
MEX_HEADER = '%%MatrixMarket matrix coordinate integer general'
EXPORTER = {'software_version': 'cellranger-6.1.1', 'format_version': 2}
AUTH_STATUS = 'ROOT_AUTHORIZED_ONE_FIXED_BACH2_NAIVE_RNA_EXTRACTION'
# Structural allowlist, including nested keys. No passthrough dictionaries.
PUBLIC_SCHEMA = {
    'schema_version': 'integer:1', 'status': 'literal:fixed_endpoint_complete',
    'endpoint': 'literal:ENSG00000112182/BACH2; fixed pooled C0+C1 naive CD4; eight source donors',
    'quantity': 'literal:mean donor log2(1+1000000*released_BACH2_RNA/all36601_released_RNA_units)',
    'groups': {'SNP': {'n': 'integer:4', 'mean': 'decimal-string', 'sample_sd': 'decimal-string'},
               'noSNP': {'n': 'integer:4', 'mean': 'decimal-string', 'sample_sd': 'decimal-string'}},
    'contrast': {'SNP_minus_noSNP': 'decimal-string', 'SE': 'decimal-string',
                 'df': 'decimal-string-or-null', 'CI95': 'two-decimal-strings-or-null',
                 'P_two_sided': 'decimal-string-or-null', 'inference_status': 'enum:unavailable_zero_SE|unavailable_nonfinite_inference|CI_available_P_unavailable_numeric_tail|available_nominal_Welch'},
    'leave_one_donor_out': {'count': 'integer:8', 'minimum': 'decimal-string', 'maximum': 'decimal-string',
        'positive_count': 'integer', 'negative_count': 'integer', 'zero_count': 'integer',
        'strict_sign_reversal_count': 'integer', 'full_point_sign': 'integer:-1/0/1'},
    'QC': {'donors': 'integer:8', 'selected_cells': 'integer:16659', 'retained_cells': 'integer:55460',
        'RNA_rows': 'integer:36601', 'zero_BACH2_donors_retained': 'integer',
        'BACH2_positive_selected_cells': 'integer', 'source_and_membership_gates': 'literal:passed'},
    'provenance': {'plan_sha256': 'sha256', 'selection_sha256': 'sha256', 'source_manifest_sha256': 'sha256',
        'implementation_sha256': 'sha256', 'runtime_sha256': 'sha256', 'public_schema_sha256': 'sha256',
        'authorization_sha256': 'sha256', 'candidate_sha256': 'sha256'},
}


class GateError(ValueError):
    """A fixed gate failed. No donor or endpoint repair is allowed."""


def require(condition, code):
    if not condition:
        raise GateError(code)


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def sha_bytes(value):
    return hashlib.sha256(value).hexdigest()


def source_path(relative, root=None):
    """No lexical traversal, absolute paths, or symlink components."""
    root = ROOT if root is None else Path(root)
    require(isinstance(relative, str) and '\\' not in relative, 'unsafe_source_path')
    require(relative == PurePosixPath(relative).as_posix(), 'noncanonical_source_path')
    parts = PurePosixPath(relative).parts
    require(parts and not PurePosixPath(relative).is_absolute() and all(p not in ('.', '..') for p in parts), 'unsafe_source_path')
    path = root
    for part in parts:
        path = path / part
        require(not path.is_symlink(), 'symlink_source_path')
    require(path.resolve().is_relative_to(root.resolve()), 'source_path_escape')
    return path


def read_pin(relative, expected_sha=None):
    """Annotation, implementation or receipt ONLY. Never called on a matrix."""
    data = source_path(relative).read_bytes()
    pin = {'path': relative, 'bytes': len(data), 'sha256': sha_bytes(data)}
    require(expected_sha is None or pin['sha256'] == expected_sha, 'annotation_pin_mismatch')
    return data, pin


def verify_distribution(distribution):
    """Verify installed wheel hashes, not just a possibly stale RECORD digest."""
    import base64
    record = distribution.read_text('RECORD')
    require(record is not None, 'runtime_record_missing')
    verified = 0
    for relative, expected_hash, expected_size in csv.reader(io.StringIO(record)):
        if not expected_hash:
            require(relative.endswith('/RECORD') or relative.endswith('.pyc'), 'runtime_unhashed_file')
            continue
        require(expected_hash.startswith('sha256='), 'runtime_hash_algorithm')
        body = Path(distribution.locate_file(relative)).read_bytes()
        encoded = base64.urlsafe_b64encode(hashlib.sha256(body).digest()).decode().rstrip('=')
        require(encoded == expected_hash.removeprefix('sha256=') and str(len(body)) == expected_size, 'runtime_installed_file_tampered')
        verified += 1
    return {'version': distribution.version, 'RECORD_sha256': sha_bytes(record.encode()), 'verified_installed_files': verified}


def runtime_snapshot():
    """Native versions, executable bytes and verified installed distributions."""
    import numpy
    import scipy
    packages = {}
    for name, module in (('numpy', numpy), ('scipy', scipy)):
        distribution = importlib.metadata.distribution(name)
        require(Path(module.__file__).resolve() == Path(distribution.locate_file(name + '/__init__.py')).resolve(), 'runtime_import_shadowed')
        packages[name] = verify_distribution(distribution)
    executable = Path(sys.executable)
    return {'python': sys.version, 'cache_tag': sys.implementation.cache_tag,
            'executable': str(executable.absolute()), 'executable_sha256': sha_bytes(executable.read_bytes()),
            'packages': packages}


def parse_donors(text):
    rows = list(csv.reader(io.StringIO(text), delimiter='\t'))
    require(len(rows) == 9, 'eight_SDRF_rows_required')
    header, rows = rows[0], rows[1:]
    needed = ('Characteristics[individual]', 'Characteristics[sex]', 'Characteristics[genotype]', 'Factor Value[genotype]')
    require(all(header.count(k) == 1 for k in needed), 'SDRF_columns')
    file_columns = [i for i, name in enumerate(header) if name == 'Derived Array Data File']
    donors = {}
    for row in rows:
        require(len(row) == len(header), 'SDRF_width')
        fields = {key: row[header.index(key)] for key in needed}
        matches = [re.fullmatch(r'(sample0[1-8])_filtered_feature_bc_matrix\.tar\.gz', row[i]) for i in file_columns]
        samples = [m.group(1) for m in matches if m]
        require(len(samples) == 1 and samples[0] not in donors, 'SDRF_sample_mapping')
        group = fields['Characteristics[genotype]']
        require(group in GROUPS and group == fields['Factor Value[genotype]'], 'SDRF_genotype')
        require(fields['Characteristics[sex]'] == 'male', 'SDRF_sex')
        donors[samples[0]] = {'donor': fields['Characteristics[individual]'], 'group': group}
    require(tuple(sorted(donors)) == SAMPLES, 'eight_source_samples_required')
    require(len({d['donor'] for d in donors.values()}) == 8, 'eight_unique_source_donors_required')
    require(Counter(d['group'] for d in donors.values()) == {'SNP': 4, 'noSNP': 4}, 'four_four_required')
    return donors


def parse_features(text, target_row=TARGET_ROW):
    features = [tuple(line.split('\t')) for line in text.splitlines()]
    require(features and all(len(f) == 3 and all(f) for f in features), 'feature_axis_width')
    require(len({f[0] for f in features}) == len(features), 'duplicate_feature_ID')
    require(all(f[2] in ('Gene Expression', 'Antibody Capture') for f in features), 'unknown_feature_type')
    require([i for i, f in enumerate(features, 1) if f == TARGET] == [target_row], 'exact_BACH2_row_required')
    require([f for f in features if f[0] == TARGET[0] or f[1] == TARGET[1]] == [TARGET], 'ambiguous_BACH2_identity')
    return features


def parse_barcodes(text):
    barcodes = text.splitlines()
    require(barcodes and all(re.fullmatch(r'[ACGT]+-[0-9]+', b) for b in barcodes), 'barcode_syntax')
    require(len(set(barcodes)) == len(barcodes), 'duplicate_source_barcode')
    return barcodes


def select_metadata(rows, donors, barcode_axes):
    """Pure annotation selection; sample identity is never stripped from keys."""
    require(set(donors) == set(barcode_axes), 'sample_axis_membership')
    positions = {s: {b: i for i, b in enumerate(axis, 1)} for s, axis in barcode_axes.items()}
    require(all(len(positions[s]) == len(barcode_axes[s]) for s in donors), 'duplicate_source_barcode')
    seen = set()
    selection = {s: {'sample': s, **donors[s], 'retained_count': 0, 'selected': [], 'state_counts': {'C0': 0, 'C1': 0}} for s in sorted(donors)}
    retained = []
    for row in rows:
        sample = row['sample_name']
        require(sample in donors and row['condition'] == donors[sample]['group'] and row['sex'] == 'male', 'metadata_donor_group')
        prefix = sample + '_'
        require(row['barcode'].startswith(prefix), 'exact_sample_prefix_required')
        barcode = row['barcode'][len(prefix):]
        key = (sample, barcode)
        require(key not in seen and barcode in positions[sample], 'fixed_retained_join_failed')
        seen.add(key)
        retained.append([sample, barcode])
        selection[sample]['retained_count'] += 1
        label = row['paper_clusters']
        is_selected = label in STATES
        if is_selected or row['paper_clusters_short'] in ('C0', 'C1') or row['seurat_clusters'] in ('0', '1'):
            require(is_selected and (row['paper_clusters_short'], row['seurat_clusters']) == STATES[label], 'state_crosswalk_changed')
        if is_selected:
            short = STATES[label][0]
            selection[sample]['state_counts'][short] += 1
            selection[sample]['selected'].append({'barcode': barcode, 'column_1based': positions[sample][barcode], 'state': short})
    require(all(d['retained_count'] > 0 and d['selected'] for d in selection.values()), 'whole_donor_membership_required')
    return {'retained_keys': retained, 'donors': list(selection.values())}


def reference_identity(text):
    rows = list(csv.DictReader(io.StringIO(text), delimiter='\t'))
    matches = [r for r in rows if r['Gene stable ID'] == TARGET[0] or r['HGNC symbol'] == TARGET[1]]
    require(matches and all(r['Gene stable ID'] == TARGET[0] and r['HGNC symbol'] == TARGET[1] for r in matches), 'reference_identity_conflict')
    return {'Ensembl_release': 116, 'gene_ID': TARGET[0], 'HGNC_symbol': TARGET[1],
            'reciprocal_relationships': len(matches), 'annotation_only': True,
            'source_executed_build_or_variant_alleles_reconciled': False,
            'original_retrieval_UTC': '2026-09-12T17:40:41Z',
            'source_URL': 'https://jun2026.archive.ensembl.org/biomart/martservice',
            'cache_provenance': 'data/derived/psc-liver-program-inputs.json:source_provenance.reference_raw_response_reconstruction; exact historical bytes reconstructed, no new retrieval'}



def validate_private_ignore(text):
    patterns = [line.strip() for line in text.splitlines() if line.strip() and not line.lstrip().startswith('#')]
    require('work/' in patterns and all(not p.startswith('!') or p == '!.env.example' for p in patterns), 'ignored_private_work_policy_changed')


def validate_public_result(result):
    """Recursive key, scalar type, finite decimal and fixed-length checks."""
    def decimal(value):
        require(type(value) is str and len(value) <= PRECISION + 32 and re.fullmatch(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', value) is not None, 'public_decimal_type')
        try:
            require(Decimal(value).is_finite(), 'public_decimal_nonfinite')
        except InvalidOperation as exc:
            raise GateError('public_decimal_invalid') from exc
    def check(value, specification):
        if isinstance(specification, dict):
            require(type(value) is dict and set(value) == set(specification), 'public_object_keys')
            for key, child in specification.items():
                check(value[key], child)
        elif specification.startswith('literal:'):
            require(type(value) is str and value == specification[len('literal:'):], 'public_literal')
        elif specification.startswith('enum:'):
            require(type(value) is str and value in specification[len('enum:'):].split('|'), 'public_enum')
        elif specification.startswith('integer'):
            require(type(value) is int, 'public_integer_type')
            if ':' in specification:
                require(value in [int(v) for v in specification.split(':')[1].split('/')], 'public_fixed_integer')
            else:
                require(value >= 0, 'public_nonnegative_integer')
        elif specification == 'sha256':
            require(type(value) is str and re.fullmatch('[a-f0-9]{64}', value) is not None, 'public_sha256')
        elif specification in ('decimal-string', 'decimal-string-or-null'):
            if value is not None or specification == 'decimal-string':
                decimal(value)
        elif specification == 'two-decimal-strings-or-null':
            if value is not None:
                require(type(value) is list and len(value) == 2, 'public_CI_shape')
                for bound in value:
                    decimal(bound)
                require(Decimal(value[0]) <= Decimal(value[1]), 'public_CI_order')
        else:
            raise GateError('unknown_public_schema_type')
    check(result, PUBLIC_SCHEMA)
    loo, qc, contrast = result['leave_one_donor_out'], result['QC'], result['contrast']
    require(sum(loo[k] for k in ('positive_count', 'negative_count', 'zero_count')) == 8, 'public_LOO_count')
    require(loo['strict_sign_reversal_count'] == (loo['negative_count'] if loo['full_point_sign'] > 0 else loo['positive_count'] if loo['full_point_sign'] < 0 else 0), 'public_LOO_sign_summary')
    require(Decimal(loo['minimum']) <= Decimal(loo['maximum']), 'public_LOO_range')
    require(qc['zero_BACH2_donors_retained'] <= 8 and qc['BACH2_positive_selected_cells'] <= qc['selected_cells'], 'public_QC_bounds')
    require(Decimal(contrast['SE']) >= 0, 'public_negative_SE')
    for summary in result['groups'].values():
        require(0 <= Decimal(summary['mean']) < 20 and Decimal(summary['sample_sd']) >= 0, 'public_group_bounds')
    if contrast['df'] is not None:
        require(3 <= Decimal(contrast['df']) <= 6, 'public_df_bounds')
    if contrast['P_two_sided'] is not None:
        require(0 < Decimal(contrast['P_two_sided']) <= 1, 'public_P_bounds')
    expected_missing = {
        'unavailable_zero_SE': (True, True), 'unavailable_nonfinite_inference': (True, True),
        'CI_available_P_unavailable_numeric_tail': (False, True), 'available_nominal_Welch': (False, False),
    }
    require((contrast['CI95'] is None, contrast['P_two_sided'] is None) == expected_missing[contrast['inference_status']], 'public_inference_availability')
    if contrast['inference_status'] == 'unavailable_zero_SE':
        require(Decimal(contrast['SE']) == 0 and contrast['df'] is None, 'public_zero_SE_status')
    else:
        require(Decimal(contrast['SE']) > 0 and contrast['df'] is not None, 'public_positive_SE_status')
    return result


def build_source_plan():
    """Read only metadata, features, barcodes, receipts, and local code/runtime."""
    pins = []
    def read(relative, expected=None):
        data, pin = read_pin(relative, expected)
        pins.append(pin)
        return data
    validate_private_ignore(read('.gitignore').decode())
    candidate = json.loads(read(CANDIDATE, CANDIDATE_SHA))
    read(REVIEW, REVIEW_SHA)
    decision = json.loads(read(DECISION))
    require(decision['status'] == 'APPROVED_METHOD_PREPARATION_ONLY_NOT_EXPRESSION_EXECUTION'
            and decision['expression_authorized'] is False and decision['candidate_sha256'] == CANDIDATE_SHA,
            'preparation_decision_changed')
    donors = parse_donors(read(SDRF, SDRF_SHA).decode())
    rows = list(csv.DictReader(io.StringIO(read(METADATA, META_SHA).decode()), delimiter='\t'))
    require(len(rows) == 55460, 'fixed_55460_required')
    reference = reference_identity(read(REFERENCE, REFERENCE_SHA).decode())
    axes, barcode_axes, matrices, archives = [], {}, [], []
    for sample in SAMPLES:
        member_base = f'work/psc-bach2-cohort-qualification/{sample}/members'
        for kind in ('features', 'barcodes', 'matrix'):
            receipt = json.loads(read(f'{member_base}/{kind}-integrity.json'))
            require(receipt['nested_gzip_CRC_ISIZE_and_EOF_passed'] is True and receipt['source_payload_matches_tar'] is True, 'closed_member_receipt_failed')
            expected_path = f'{member_base}/{kind}.txt'
            require(receipt['decoded_path'] == expected_path, 'member_receipt_path')
            require(isinstance(receipt['decoded_bytes'], int) and receipt['decoded_bytes'] > 0 and re.fullmatch('[a-f0-9]{64}', receipt['decoded_sha256']), 'member_receipt_pin')
            if kind == 'matrix':
                # NO matrix open, stat, or hashing in preparation/validation.
                matrices.append({'sample': sample, 'path': expected_path, 'bytes': receipt['decoded_bytes'], 'sha256': receipt['decoded_sha256']})
            else:
                data = read(expected_path, receipt['decoded_sha256'])
                require(len(data) == receipt['decoded_bytes'], 'annotation_byte_count')
                if kind == 'features':
                    require(sha_bytes(data) == FEATURE_SHA, 'complete_feature_axis_pin')
                    features = parse_features(data.decode())
                    require(Counter(f[2] for f in features) == {'Gene Expression': 36601, 'Antibody Capture': 38}, 'full_source_feature_universe')
                    axes.append(features)
                else:
                    barcode_axes[sample] = parse_barcodes(data.decode())
        filename = sample + '_filtered_feature_bc_matrix.tar.gz'
        archive_path = 'data/raw/psc-bach2-qualification/' + ('counts/' if sample == 'sample05' else 'cohort/') + filename
        receipt_path = ('data/raw/psc-bach2-qualification/counts/acquisition-receipt.json'
                        if sample == 'sample05' else archive_path + '.receipt.json')
        acq = json.loads(read(receipt_path))
        expected_bytes = acq['authorized_archive_bytes'] if sample == 'sample05' else acq['expected_bytes']
        require(acq['complete'] is True and acq['body_bytes'] == expected_bytes, 'closed_acquisition_receipt')
        if sample == 'sample05':
            require(acq['archive_path'] == archive_path, 'reused_archive_path')
        archives.append({'sample': sample, 'path': archive_path, 'bytes': acq['body_bytes'], 'sha256': acq['sha256'],
                         'source_URL': acq['url'], 'retrieved_at_UTC': acq['completed_at_utc']})
    require(all(axis == axes[0] for axis in axes), 'ordered_axes_changed')
    selection = select_metadata(rows, donors, barcode_axes)
    states = Counter(entry['state'] for donor in selection['donors'] for entry in donor['selected'])
    require(states == {'C0': 9320, 'C1': 7339}, 'fixed_16659_union_required')
    require(all(all(d['state_counts'][s] > 0 for s in ('C0', 'C1')) for d in selection['donors']), 'fixed_state_membership')
    for matrix in matrices:
        matrix['rows'] = len(axes[0])
        matrix['columns'] = len(barcode_axes[matrix['sample']])
    implementation = []
    for relative in ('scripts/analyze_psc_bach2_naive_rna.py', 'tests/test_psc_bach2_naive_rna.py'):
        _, pin = read_pin(relative)
        implementation.append(pin)
    manifest = {'annotation_and_receipt_pins': pins, 'matrix_expected_pins_NOT_read': matrices, 'archive_expected_pins_NOT_read': archives}
    plan = {'schema_version': 1, 'status': 'METHOD_PREPARATION_ONLY_NOT_EXECUTION_AUTHORIZATION',
            'candidate_sha256': CANDIDATE_SHA, 'fixed_candidate': candidate, 'reference_identity': reference,
            'implementation': implementation, 'runtime': runtime_snapshot(),
            'source_manifest': manifest, 'source_manifest_sha256': digest(manifest),
            'selection_sha256': digest(selection), 'public_schema': PUBLIC_SCHEMA, 'public_schema_sha256': digest(PUBLIC_SCHEMA),
            'matrix_header': MEX_HEADER, 'exporter_comment': EXPORTER,
            'nnz_gate': 'Declared nnz is bound by frozen matrix SHA; require exact streamed count and strict column-major uniqueness',
            'RNA_rows_1based': [i for i, f in enumerate(axes[0], 1) if f[2] == 'Gene Expression'],
            'numerics': {'Decimal_precision': PRECISION, 'maximum_value_digits': MAX_VALUE_DIGITS, 'maximum_coordinates_per_matrix': MAX_COORDINATES,
                'variance': 'exact Fraction moments of fixed-precision Decimal Y; no constancy tolerance',
                'point_and_LOO': 'algebraically identical exact rational products of 1+1e6B/R; stable logarithm preserves true zero/sign',
                'Student_t': 'native SciPy actual Welch df; P tail underflow is unavailable, not zero', 'numeric_JSON': 'finite Decimal strings; working precision is not biological precision'},
            'closed_limits': ['released RNA source-count-format units, not certified untouched UMIs', 'no source build/alleles or protein/causal interpretation', 'same-cohort pooled mixture; no clinical adjustment or extra endpoint'],
            'permitted_future_run_root': BASE + '/runs/',
            'future_run_inventory': list(RUN_INVENTORY), 'completion_status': COMPLETION_STATUS}
    return plan, selection


def fresh_directory(relative, kind='run', create=False):
    base = BASE + ('/runs/' if kind == 'run' else '/')
    require(isinstance(relative, str) and relative.startswith(base), 'unapproved_output_root')
    name = relative[len(base):]
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,79}', name) is not None, 'unsafe_output_name')
    if kind != 'run':
        require(name not in ('runs', 'independent-review'), 'reserved_output_name')
    path = source_path(relative)
    require(not path.exists(), 'output_must_be_fresh')
    if create:
        path.parent.mkdir(parents=True, exist_ok=True)
        # Recheck after parent creation; exclusive mkdir makes authorization one-use.
        source_path(relative)
        path.mkdir(mode=0o700)
    return path


def compile_plan(relative=DEFAULT_PLAN):
    fresh_directory(relative, kind='plan')
    plan, selection = build_source_plan()
    destination = fresh_directory(relative, kind='plan', create=True)
    for name, value in (('plan.json', plan), ('selection.private.json', selection), ('public-schema.json', PUBLIC_SCHEMA)):
        (destination / name).write_bytes(canonical(value))
    return {'status': 'source_only_plan_compiled', 'plan_sha256': digest(plan), 'selection_sha256': digest(selection), 'matrix_bodies_opened': 0, 'measurements_executed': False}



def plan_directory_path(relative):
    """Only named preparation plans may be read, never source or run roles."""
    require(type(relative) is str and relative.startswith(BASE + '/'), 'plan_requires_preparation_namespace')
    name = relative[len(BASE) + 1:]
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,79}', name) is not None
            and name not in ('runs', 'independent-review'), 'invalid_plan_directory_role')
    return source_path(relative)


def validate_plan(relative=DEFAULT_PLAN):
    directory = plan_directory_path(relative)
    plan_bytes = source_path(relative + '/plan.json').read_bytes()
    selection_bytes = source_path(relative + '/selection.private.json').read_bytes()
    plan, selection = build_source_plan()
    require(plan_bytes == canonical(plan), 'compiled_plan_or_current_inputs_changed')
    require(selection_bytes == canonical(selection), 'compiled_selection_changed')
    require(source_path(relative + '/public-schema.json').read_bytes() == canonical(PUBLIC_SCHEMA), 'public_schema_changed')
    return plan, selection


def check_authorization(auth, supplied_sha, raw_sha, plan, output):
    require(re.fullmatch('[a-f0-9]{64}', supplied_sha or '') is not None and supplied_sha == raw_sha, 'authorization_SHA_required_or_mismatch')
    expected = {'schema_version': 1, 'status': AUTH_STATUS, 'expression_authorized': True,
        'authorizer': 'root', 'one_fixed_extraction_only': True, 'output_directory': output,
        'plan_sha256': digest(plan), 'selection_sha256': plan['selection_sha256'],
        'source_manifest_sha256': plan['source_manifest_sha256'],
        'implementation_sha256': plan['implementation'][0]['sha256'], 'runtime_sha256': digest(plan['runtime']),
        'public_schema_sha256': plan['public_schema_sha256'], 'candidate_sha256': CANDIDATE_SHA}
    require(canonical(auth) == canonical(expected), 'root_frozen_authorization_binding_failed')
    fresh_directory(output)
    return expected


def exact_positive_integer(token):
    require(isinstance(token, str) and len(token) <= 100 and re.fullmatch(r'[+]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?', token, flags=re.ASCII) is not None, 'invalid_exact_integer_token')
    try:
        value = Decimal(token)
        require(value.is_finite() and 0 <= value.adjusted() < MAX_VALUE_DIGITS and value == value.to_integral_value() and value > 0, 'nonpositive_nonintegral_or_oversized_value')
        return int(value)
    except (InvalidOperation, OverflowError, ValueError) as exc:
        raise GateError('invalid_exact_integer_token') from exc


def stream_selected_matrix(stream, pin, selected_columns, RNA_rows, target_row=TARGET_ROW):
    """Future authorized single pass only; tests use synthetic bytes only.

    Parse all coordinates for exact integrity, but accumulate no ADT, unselected
    cell, or other gene-level measurement. Matrix SHA must pass before any Y.
    """
    require(target_row in RNA_rows and len(set(RNA_rows)) == len(RNA_rows), 'RNA_target_or_axis')
    require(selected_columns and len(set(selected_columns)) == len(selected_columns) and all(type(c) is int and 1 <= c <= pin['columns'] for c in selected_columns), 'selected_column_axis')
    RNA_rows, selected_columns = set(RNA_rows), set(selected_columns)
    hasher, size = hashlib.sha256(), 0
    def line():
        nonlocal size
        raw = stream.readline(513)
        require(len(raw) <= 512, 'matrix_line_bound')
        size += len(raw)
        require(size <= pin['bytes'], 'matrix_byte_bound')
        hasher.update(raw)
        try:
            return raw.decode('ascii').rstrip('\r\n') if raw else None
        except UnicodeDecodeError as exc:
            raise GateError('matrix_not_ASCII') from exc
    require(line() == MEX_HEADER, 'matrix_header_changed')
    comment = line()
    require(comment is not None and comment.startswith('%metadata_json: '), 'exporter_comment_missing')
    require(json.loads(comment[len('%metadata_json: '):]) == EXPORTER, 'exporter_comment_changed')
    dimensions = line()
    require(dimensions is not None and re.fullmatch(r'[0-9]+ [0-9]+ [0-9]+', dimensions) is not None, 'matrix_dimensions_syntax')
    nr, nc, nnz = map(int, dimensions.split())
    require((nr, nc) == (pin['rows'], pin['columns']) and 0 <= nnz <= min(MAX_COORDINATES, nr * nc), 'matrix_dimensions_changed')
    require(all(type(r) is int and 1 <= r <= nr for r in RNA_rows), 'RNA_rows_out_of_bounds')
    B, R, positive_cells, coordinates, previous = 0, 0, 0, 0, (0, 0)
    while (text := line()) is not None:
        fields = text.split()
        require(len(fields) == 3 and all(re.fullmatch('[0-9]{1,9}', t) for t in fields[:2]), 'matrix_coordinate_syntax')
        row, col = map(int, fields[:2])
        require(1 <= row <= nr and 1 <= col <= nc and (col, row) > previous, 'matrix_bounds_duplicate_or_order')
        previous = (col, row)
        value = exact_positive_integer(fields[2])
        coordinates += 1
        require(coordinates <= nnz, 'matrix_nnz_exceeded')
        if col in selected_columns and row in RNA_rows:
            R += value
            if row == target_row:
                B += value
                positive_cells += 1
    require(coordinates == nnz and size == pin['bytes'] and hasher.hexdigest() == pin['sha256'], 'matrix_frozen_hash_bytes_or_nnz_failed')
    return {'B': B, 'R': R, 'B_positive_cells': positive_cells, 'selected_cells': len(selected_columns)}


def as_decimal(value):
    with localcontext(Context(prec=PRECISION)) as ctx:
        ctx.prec = PRECISION
        if isinstance(value, Fraction):
            return Decimal(value.numerator) / Decimal(value.denominator)
        return +Decimal(value)


def decimal_string(value):
    value = as_decimal(value)
    require(value.is_finite(), 'nonfinite_serialization')
    return str(value) if value else '0'



def log_positive_fraction(value):
    """Stable ln of an exact positive rational, including arbitrarily near one.

    The atanh series is an algebraic evaluation rule, not a variance threshold
    or a change of estimand. It avoids losing a tiny exact numerator difference
    when a rational near one would round to Decimal(1).
    """
    require(isinstance(value, Fraction) and value > 0, 'positive_log_argument_required')
    if value == 1:
        return Decimal(0)
    with localcontext(Context(prec=PRECISION)):
        if Fraction(1, 2) <= value <= Fraction(3, 2):
            x = Decimal(value.numerator - value.denominator) / Decimal(value.denominator)
            z = x / (2 + x)
            squared, power, total = z * z, z, z
            for k in range(1, 2048):
                power *= squared
                updated = total + power / (2 * k + 1)
                if updated == total:
                    return 2 * total
                total = updated
            raise GateError('fixed_log_series_failed')
        return (Decimal(value.numerator) / Decimal(value.denominator)).ln()


def logarithmic_means_and_contrast(arguments):
    """Same mean-log estimand, exact product identity gives stable zero/sign.

    For n1/n0=4/4 or3/4, Delta=ln(prod(Q1)^n0/prod(Q0)^n1)/(n1*n0*ln2).
    Q=1+1e6*B/R. No additional endpoint or empirical tuning is introduced.
    """
    require(set(arguments) == set(GROUPS) and all(len(arguments[g]) in (3, 4) for g in GROUPS)
            and sum(map(len, arguments.values())) in (7, 8), 'fixed_log_contrast_sizes')
    products = {g: math.prod(arguments[g], start=Fraction(1)) for g in GROUPS}
    require(all(isinstance(q, Fraction) and 1 <= q <= 1_000_001 for group in arguments.values() for q in group), 'bounded_source_log_argument')
    n1, n0 = len(arguments['SNP']), len(arguments['noSNP'])
    with localcontext(Context(prec=PRECISION)):
        log2 = Decimal(2).ln()
        means = {g: log_positive_fraction(products[g]) / len(arguments[g]) / log2 for g in GROUPS}
        delta = log_positive_fraction(products['SNP'] ** n0 / products['noSNP'] ** n1) / (n1 * n0) / log2
    return means, delta


def score_units(B, R):
    require(type(B) is int and type(R) is int and 0 <= B <= R and R > 0, 'whole_endpoint_invalid_RNA_units')
    require(R < MAX_COORDINATES * 10 ** MAX_VALUE_DIGITS, 'RNA_accumulation_bound')
    if B == 0:
        return Decimal(0)
    with localcontext(Context(prec=PRECISION)) as ctx:
        ctx.prec = PRECISION
        score = log_positive_fraction(1 + Fraction(1_000_000 * B, R)) / Decimal(2).ln()
    require(score.is_finite() and score > 0, 'score_precision_failed')
    return score


def moment(values):
    require(len(values) >= 2 and all(isinstance(x, (Decimal, Fraction)) for x in values), 'finite_exact_scores_required')
    try:
        values = [Fraction(x) for x in values]
    except (ValueError, OverflowError) as exc:
        raise GateError('nonfinite_score') from exc
    n = len(values)
    mean = sum(values, Fraction()) / n
    variance = sum(((x - y) ** 2 for i, x in enumerate(values) for y in values[i + 1:]), Fraction()) / (n * (n - 1))
    return mean, variance


def welch(scores, *, exact_log_arguments=None):
    """Exactly 4/4 donor scores. Exact rational moments, fixed Decimal sqrt."""
    require(set(scores) == set(GROUPS) and all(len(scores[g]) == 4 for g in GROUPS), 'whole_eight_donor_gate')
    moments = {g: moment(scores[g]) for g in GROUPS}
    means = {g: moments[g][0] for g in GROUPS}
    variances = {g: moments[g][1] for g in GROUPS}
    delta = means['SNP'] - means['noSNP']
    if exact_log_arguments is not None:
        require(all(len(exact_log_arguments[g]) == 4 for g in GROUPS), 'whole_eight_log_argument_gate')
        means, delta = logarithmic_means_and_contrast(exact_log_arguments)
    variance_SE = sum(variances.values(), Fraction()) / 4
    with localcontext(Context(prec=PRECISION)) as ctx:
        ctx.prec = PRECISION
        SE = as_decimal(variance_SE).sqrt()
        result = {'delta': decimal_string(delta), 'SE': decimal_string(SE), 'df': None, 'CI95': None, 'P': None,
                  'inference_status': 'unavailable_zero_SE',
                  'groups': {g: {'n': 4, 'mean': decimal_string(means[g]), 'sample_sd': decimal_string(as_decimal(variances[g]).sqrt())} for g in GROUPS}}
        if not variance_SE:
            return result
        v1, v0 = variances['SNP'] / 4, variances['noSNP'] / 4
        df = (v1 + v0) ** 2 / (v1 ** 2 / 3 + v0 ** 2 / 3)
        require(Fraction(3) <= df <= Fraction(6), 'invalid_Welch_df')
        result['df'] = decimal_string(df)
        from scipy.stats import t
        critical = float(t.ppf(0.975, float(df)))
        if not math.isfinite(critical) or critical <= 0 or not SE.is_finite() or SE <= 0:
            result['inference_status'] = 'unavailable_nonfinite_inference'
            return result
        half = Decimal.from_float(critical) * SE
        centre = as_decimal(delta)
        result['CI95'] = [decimal_string(centre - half), decimal_string(centre + half)]
        statistic = abs(centre / SE)
        statistic_float = float(statistic)
        probability = 2 * float(t.sf(statistic_float, float(df))) if math.isfinite(statistic_float) else float('nan')
        if math.isfinite(probability) and 0 < probability <= 1:
            result['P'] = str(probability)
            result['inference_status'] = 'available_nominal_Welch'
        else:
            result['inference_status'] = 'CI_available_P_unavailable_numeric_tail'
        return result


def sign(x):
    return int(x > 0) - int(x < 0)


def leaveouts(records, *, exact_log_arguments=None):
    require(len(records) == 8 and len({r['sample'] for r in records}) == 8 and Counter(r['group'] for r in records) == {'SNP': 4, 'noSNP': 4}, 'whole_eight_LOO_gate')
    scores = {g: [Fraction(r['Y']) for r in records if r['group'] == g] for g in GROUPS}
    full = sum(scores['SNP']) / 4 - sum(scores['noSNP']) / 4
    if exact_log_arguments is not None:
        require(set(exact_log_arguments) == {r['sample'] for r in records}, 'LOO_log_argument_membership')
        _, full = logarithmic_means_and_contrast({g: [exact_log_arguments[r['sample']] for r in records if r['group'] == g] for g in GROUPS})
    private, values = [], []
    for omitted in records:
        kept = [r for r in records if r['sample'] != omitted['sample']]
        means = {g: sum((Fraction(r['Y']) for r in kept if r['group'] == g), Fraction()) / sum(r['group'] == g for r in kept) for g in GROUPS}
        delta = means['SNP'] - means['noSNP']
        if exact_log_arguments is not None:
            _, delta = logarithmic_means_and_contrast({g: [exact_log_arguments[r['sample']] for r in kept if r['group'] == g] for g in GROUPS})
        values.append(delta)
        private.append({'omitted_sample': omitted['sample'], 'delta': decimal_string(delta)})
    summary = {'count': 8, 'minimum': decimal_string(min(values)), 'maximum': decimal_string(max(values)),
        'positive_count': sum(v > 0 for v in values), 'negative_count': sum(v < 0 for v in values), 'zero_count': sum(v == 0 for v in values),
        'strict_sign_reversal_count': sum(sign(v) * sign(full) < 0 for v in values), 'full_point_sign': sign(full)}
    return private, summary


def analyze_donors(measurements, selection):
    donors = selection['donors']
    require(len(donors) == 8 and {d['sample'] for d in donors} == set(SAMPLES) and len({d['donor'] for d in donors}) == 8, 'whole_eight_identity_gate')
    require(Counter(d['group'] for d in donors) == {'SNP': 4, 'noSNP': 4}, 'whole_four_four_gate')
    require(set(measurements) == set(SAMPLES), 'whole_eight_measurement_gate')
    records, ratios, scores = [], {g: [] for g in GROUPS}, {g: [] for g in GROUPS}
    for donor in donors:
        m = measurements[donor['sample']]
        n = len(donor['selected'])
        require(type(m['B_positive_cells']) is int and type(m['selected_cells']) is int and 0 <= m['B_positive_cells'] <= n and m['selected_cells'] == n, 'whole_selected_membership_gate')
        require((m['B'] == 0) == (m['B_positive_cells'] == 0), 'B_detection_integrity')
        Y = score_units(m['B'], m['R'])
        ratios[donor['group']].append(Fraction(m['B'], m['R']))
        scores[donor['group']].append(Y)
        records.append({'sample': donor['sample'], 'donor': donor['donor'], 'group': donor['group'],
            'selected_cells': n, 'state_cell_counts': donor['state_counts'],
            'state_cell_fractions': {s: decimal_string(Fraction(c, n)) for s, c in donor['state_counts'].items()},
            'B_source_units': str(m['B']), 'RNA_source_units': str(m['R']), 'B_positive_cells': m['B_positive_cells'],
            'B_positive_cell_fraction': decimal_string(Fraction(m['B_positive_cells'], n)), 'Y': decimal_string(Y)})
    # Any impossible rounding collision fails, rather than declaring constancy.
    require(all(len(set(ratios[g])) == len(set(scores[g])) for g in GROUPS), 'score_precision_collision')
    arguments = {s: 1 + Fraction(1_000_000 * m['B'], m['R']) for s, m in measurements.items()}
    grouped_arguments = {g: [arguments[r['sample']] for r in records if r['group'] == g] for g in GROUPS}
    inference = welch(scores, exact_log_arguments=grouped_arguments)
    private_LOO, summary_LOO = leaveouts(records, exact_log_arguments=arguments)
    return {'donors': records, 'inference': inference, 'LOO_private': private_LOO, 'LOO_summary': summary_LOO,
            'zero_B_donors': sum(m['B'] == 0 for m in measurements.values()),
            'B_positive_selected_cells': sum(m['B_positive_cells'] for m in measurements.values())}



def result_provenance(plan, authorization_sha):
    return {'plan_sha256': digest(plan), 'selection_sha256': plan['selection_sha256'],
        'source_manifest_sha256': plan['source_manifest_sha256'],
        'implementation_sha256': plan['implementation'][0]['sha256'],
        'runtime_sha256': digest(plan['runtime']), 'public_schema_sha256': plan['public_schema_sha256'],
        'authorization_sha256': authorization_sha, 'candidate_sha256': CANDIDATE_SHA}


def exact_run_inventory(directory, expected):
    directory = Path(directory)
    require(directory.is_dir() and not directory.is_symlink(), 'run_directory_invalid')
    files = list(directory.iterdir())
    require({p.name for p in files} == set(expected)
            and all(p.is_file() and not p.is_symlink() for p in files), 'run_inventory_missing_extra_or_symlink')


def payload_pins(directory):
    pins = {}
    for name in RUN_PAYLOADS:
        body = (Path(directory) / name).read_bytes()
        pins[name] = {'bytes': len(body), 'sha256': sha_bytes(body)}
    return pins


def completion_record(directory, plan, authorization_sha, output):
    return {'schema_version': 1, 'status': COMPLETION_STATUS, 'output_directory': output,
        'allowed_inventory': list(RUN_INVENTORY), 'provenance': result_provenance(plan, authorization_sha),
        'payloads': payload_pins(directory)}


def write_completion(directory, plan, authorization_sha, output):
    """Last-written durable completion. A hash inventory is not a math audit."""
    directory = Path(directory)
    exact_run_inventory(directory, RUN_PAYLOADS)
    record = completion_record(directory, plan, authorization_sha, output)
    temporary = directory / 'completion.pending'
    try:
        with temporary.open('xb') as stream:
            stream.write(canonical(record))
            stream.flush()
            import os
            os.fsync(stream.fileno())
        temporary.replace(directory / 'completion.json')
    except Exception:
        temporary.unlink(missing_ok=True)
        (directory / 'completion.json').unlink(missing_ok=True)
        raise
    return record


def validate_completed_run(directory, plan, authorization_sha, output):
    """Only for an authorized result audit or synthetic fixtures, never default."""
    directory = Path(directory)
    exact_run_inventory(directory, RUN_INVENTORY)
    raw = (directory / 'completion.json').read_bytes()
    record = json.loads(raw)
    require(raw == canonical(completion_record(directory, plan, authorization_sha, output)), 'completion_payload_or_binding_mismatch')
    return record


def public_projection(private, plan, authorization_sha):
    """Construct an allowlist. Never copy donor, join, or private LOO objects."""
    inference = private['inference']
    summary = private['LOO_summary']
    result = {'schema_version': 1, 'status': 'fixed_endpoint_complete',
        'endpoint': 'ENSG00000112182/BACH2; fixed pooled C0+C1 naive CD4; eight source donors',
        'quantity': 'mean donor log2(1+1000000*released_BACH2_RNA/all36601_released_RNA_units)',
        'groups': {g: {k: inference['groups'][g][k] for k in ('n', 'mean', 'sample_sd')} for g in GROUPS},
        'contrast': {'SNP_minus_noSNP': inference['delta'], 'SE': inference['SE'], 'df': inference['df'],
            'CI95': None if inference['CI95'] is None else list(inference['CI95']), 'P_two_sided': inference['P'], 'inference_status': inference['inference_status']},
        'leave_one_donor_out': {k: summary[k] for k in PUBLIC_SCHEMA['leave_one_donor_out']},
        'QC': {'donors': 8, 'selected_cells': 16659, 'retained_cells': 55460, 'RNA_rows': 36601,
               'zero_BACH2_donors_retained': private['zero_B_donors'], 'BACH2_positive_selected_cells': private['B_positive_selected_cells'],
               'source_and_membership_gates': 'passed'},
        'provenance': result_provenance(plan, authorization_sha)}
    return validate_public_result(result)



def read_root_authorization(relative):
    """Reject count/annotation files used in the authorization argument role."""
    require(type(relative) is str and relative.startswith(AUTHORIZATION_ROOT), 'authorization_requires_root_freeze_namespace')
    name = relative[len(AUTHORIZATION_ROOT):]
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,79}\.json', name) is not None, 'authorization_filename')
    path = source_path(relative)
    require(path.is_file() and 0 < path.stat().st_size <= MAX_AUTHORIZATION_BYTES, 'authorization_file_size')
    with path.open('rb') as stream:
        raw = stream.read(MAX_AUTHORIZATION_BYTES + 1)
    require(0 < len(raw) <= MAX_AUTHORIZATION_BYTES, 'authorization_read_bound')
    return raw


def execute(plan_directory, authorization_path, authorization_sha, output):
    require(authorization_path is not None and authorization_sha is not None and output is not None, 'separate_root_authorization_required')
    # Validate the input role BEFORE reading any body. A matrix is never auth.
    raw_auth = read_root_authorization(authorization_path)
    require(sha_bytes(raw_auth) == authorization_sha, 'authorization_SHA_mismatch')
    auth = json.loads(raw_auth)
    plan, selection = validate_plan(plan_directory)
    check_authorization(auth, authorization_sha, sha_bytes(raw_auth), plan, output)
    destination = fresh_directory(output, create=True)
    try:
        measurements = {}
        for pin in plan['source_manifest']['matrix_expected_pins_NOT_read']:
            donor = next(d for d in selection['donors'] if d['sample'] == pin['sample'])
            # This is the only production matrix-open site. Never used in prepare.
            with source_path(pin['path']).open('rb') as stream:
                measurements[pin['sample']] = stream_selected_matrix(stream, pin,
                    [c['column_1based'] for c in donor['selected']], plan['RNA_rows_1based'])
        private = analyze_donors(measurements, selection)
        public = public_projection(private, plan, authorization_sha)
        (destination / 'measurements.private.json').write_bytes(canonical(private))
        (destination / 'public-aggregate.json').write_bytes(canonical(public))
        (destination / 'authorization.used.json').write_bytes(raw_auth)
        write_completion(destination, plan, authorization_sha, output)
        validate_completed_run(destination, plan, authorization_sha, output)
        return {'status': 'authorized_fixed_endpoint_complete', 'output_directory': output}
    except Exception:
        # Keep this run consumed. Never publish a partial point or n=7 fallback.
        (destination / 'completion.json').unlink(missing_ok=True)
        (destination / 'completion.pending').unlink(missing_ok=True)
        (destination / 'failure.json').write_bytes(canonical({'status': 'whole_endpoint_unavailable', 'no_repair_or_rescue': True}))
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--validate-only', action='store_true')
    mode.add_argument('--compile-plan', action='store_true')
    mode.add_argument('--execute', action='store_true')
    parser.add_argument('--plan-directory', default=DEFAULT_PLAN)
    parser.add_argument('--authorization')
    parser.add_argument('--authorization-sha256')
    parser.add_argument('--output-directory')
    args = parser.parse_args(argv)
    try:
        if args.execute:
            result = execute(args.plan_directory, args.authorization, args.authorization_sha256, args.output_directory)
        else:
            require(not any((args.authorization, args.authorization_sha256, args.output_directory)), 'execution_options_require_execute')
            if args.compile_plan:
                result = compile_plan(args.plan_directory)
            else:
                plan, _ = validate_plan(args.plan_directory)
                result = {'status': 'source_only_validation_passed', 'plan_sha256': digest(plan),
                          'matrix_bodies_opened': 0, 'measurements_executed': False, 'execution_authorized': False}
        print(json.dumps(result, sort_keys=True))
        return 0
    except (GateError, OSError, ValueError, KeyError, TypeError) as exc:
        # No raw source values or individual identifiers in public CLI errors.
        print(json.dumps({'status': 'unavailable_no_rescue', 'exception_class': type(exc).__name__}), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
