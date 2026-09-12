import io
import sys
import tarfile
import unittest
from decimal import Decimal
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from audit_finemap_archive import (audit_dataset, crosswalk, parse_configs, parse_snps,
                                   safe_members, singleton_sets)
from prepare_controlled_comparison import assess_mapping, nominate, select_loci


SNP_HEADER = 'index snp snp_prob snp_log10bf\n'
CONFIG_HEADER = 'rank config config_prob config_log10bf\n'


def log(n, maximum=1):
    return (f'Welcome to FINEMAP v1.1\n- n-causal-max : {maximum}\n'
            f'- Number of SNPs in region : {n} SNPs\n'
            '- Number of individuals in GWAS : 123 SNPs\n'
            '- Post-Pr(# of causal SNPs is k)  :\n   0 -> 0\n   1 -> 1\n')


class FinemapTests(unittest.TestCase):
    def test_multicausal_pips_are_not_normalized_into_a_credible_set(self):
        config = parse_configs(CONFIG_HEADER + '1 rs1,rs2 0.7 1\n2 rs1 0.3 1\n', {'rs1', 'rs2'})
        self.assertEqual(config['marginals']['rs1'], 1)
        self.assertEqual(sum(config['marginals'].values()), Decimal('1.7'))
        self.assertEqual(singleton_sets(config, 2), [])

    def test_singleton_reconstruction_keeps_boundary_ties_and_requires_complete_mass(self):
        config = parse_configs(CONFIG_HEADER + '1 rs1 0.90 1\n2 rs2 0.05 1\n3 rs3 0.05 1\n', {'rs1', 'rs2', 'rs3'})
        sets = singleton_sets(config, 3)
        self.assertEqual(sets[0]['size'], 3)
        self.assertEqual(sets[0]['mass'], '1.00')
        incomplete = parse_configs(CONFIG_HEADER + '1 rs1 0.6 1\n', {'rs1', 'rs2'})
        self.assertEqual(singleton_sets(incomplete, 2), [])

    def test_log_mismatch_is_reported_even_when_pips_agree(self):
        stem = 'data_release/psc_gwas_fine_mapping/test'
        contents = {stem + '.snp': SNP_HEADER + '2 rs1 1.0000 1\n1 rs2 0.7000 1\n',
                    stem + '.config': CONFIG_HEADER + '1 rs1,rs2 0.7 1\n2 rs1 0.3 1\n',
                    stem + '.log': log(2)}
        audit, rows = audit_dataset(stem, contents)
        self.assertEqual(Decimal(audit['max_pip_config_difference']), 0)
        self.assertIn('configuration_exceeds_log_max_causal_variants', audit['issues'])
        self.assertIn('configuration_posterior_k_disagrees_with_log', audit['issues'])
        self.assertEqual(rows[0]['source_index'], '2')

    def test_missing_configuration_is_not_a_zero_probability_result(self):
        stem = 'data_release/eqtl_fine_mapping/test'
        audit, rows = audit_dataset(stem, {stem + '.snp': SNP_HEADER + '1 rs1 0.0000 -inf\n',
                                          stem + '.config': '', stem + '.log': log(1, 5)})
        self.assertIsNone(audit['config_rows'])
        self.assertIsNone(audit['max_pip_config_difference'])
        self.assertEqual(rows[0]['config_summed_marginal_probability'], '')
        self.assertEqual(rows[0]['snp_log10bf_as_published'], '-inf')

    def test_malformed_or_duplicate_variants_and_unknown_config_ids_fail(self):
        for body in ['1 rs1 nan 0\n', '1 rs1 1.2 0\n', '1 rs1 0.1 0\n2 rs1 0.2 0\n']:
            with self.assertRaises(ValueError):
                parse_snps(SNP_HEADER + body)
        for body in ['1 rs1,rs9 1 0\n', '1 rs1,rs1 1 0\n', '1 rs1 0.5 0\n2 rs1 0.5 0\n']:
            with self.assertRaises(ValueError):
                parse_configs(CONFIG_HEADER + body, {'rs1', 'rs2'})

    def test_coordinate_match_does_not_certify_alleles_or_signal_membership(self):
        signal = {'source_table_data_row': '1', 'candidate_gene_as_published': 'ETS2', 'signal': '1',
                  'highest_pp_snp': 'rs1', 'chromosome': '21', 'highest_pp_position_b37': '100',
                  'posterior_probability_as_published': '.5', 'credible_set_size_as_published': '5'}
        result = crosswalk([signal], [{'dataset_id': 'gwas/PSMG1', 'variant_id_as_published': 'chr21:100',
                                       'snp_prob_as_published': '.5'}])[0]
        self.assertIn('unverified', result['match_method'])
        self.assertFalse(result['complete_published_signal_membership_available'])
        self.assertIn('not declared equivalent', result['region_link_note'])

    def test_archive_paths_and_links_cannot_escape_reader(self):
        for name, kind in [('../outside', tarfile.REGTYPE), ('/absolute', tarfile.REGTYPE),
                           ('data/link', tarfile.SYMTYPE)]:
            buffer = io.BytesIO()
            with tarfile.open(fileobj=buffer, mode='w') as archive:
                member = tarfile.TarInfo(name)
                member.type = kind
                archive.addfile(member, io.BytesIO())
            buffer.seek(0)
            with tarfile.open(fileobj=buffer) as archive:
                with self.assertRaises(ValueError):
                    safe_members(archive)

    def test_locus_selection_excludes_prior_results_and_preserves_tie_rule(self):
        def signal(label, pip, size, signal_id='1'):
            return {'candidate_gene_as_published': label, 'posterior_probability_as_published': pip,
                    'credible_set_size_as_published': size, 'signal': signal_id,
                    'chromosome': '2', 'lead_position_b37': '100'}
        rows = [signal('prior', '.4', '1'), signal('multi', '.3', '2'), signal('multi', '.2', '2', '2'),
                signal('small', '.46', '5'), signal('tie2', '.18', '12'), signal('tie1', '.20', '12'),
                signal('below', '.09', '1'), signal('high', '.99', '1')]
        rule = {'exclude_prior_pilot_regions': ['prior'], 'minimum_published_top_pip': .1,
                'maximum_published_top_pip_exclusive': .5, 'number_of_regions': 3}
        selected, _ = select_loci(rows, rule)
        self.assertEqual([r['candidate_gene_as_published'] for r in selected], ['small', 'tie1', 'tie2'])

    def test_multiallelic_and_indel_records_cannot_be_coerced_to_a_snp(self):
        snv = {'start': 10, 'end': 10, 'allele_string': 'G/T'}
        self.assertEqual(assess_mapping(snv, {**snv, 'start': 20, 'end': 20}), '')
        multiple = {**snv, 'allele_string': 'G/A/C/T'}
        indel = {**snv, 'end': 17, 'allele_string': 'TTTTTTTT/TTTTTTTTT'}
        self.assertIn('Multiallelic', assess_mapping(multiple, multiple))
        self.assertIn('Indel', assess_mapping(indel, indel))
        self.assertIn('differ', assess_mapping(snv, {**snv, 'allele_string': 'A/C'}))

    def test_reference_mismatch_fails_and_candidate_cap_keeps_a_reserve(self):
        sources = [{'dataset_id': 'gwas/test', 'variant_id_as_published': f'rs{i}', 'identifier_type': 'rsid',
                    'source_row': i, 'snp_prob_as_published': p} for i, p in [(1, '.3'), (2, '.2'), (3, '.1')]]
        def read_json(name):
            if name.endswith('-sequence.json'):
                return {'seq': 'T'}
            _, build, identifier = name.removesuffix('.json').split('-')
            return {'name': identifier, 'most_severe_consequence': 'intron_variant',
                    'mappings': [{'assembly_name': build, 'seq_region_name': '2', 'coord_system': 'chromosome',
                                  'strand': 1, 'start': 100, 'end': 100, 'allele_string': 'T/A'}]}
        rows = nominate(sources, '2', {'maximum_per_region': 2}, read_json, lambda _: 'T')
        self.assertEqual([r['role'] for r in rows], ['anchor_pending_comparators', 'anchor_pending_comparators', 'reserve'])
        with self.assertRaises(ValueError):
            nominate(sources, '2', {'maximum_per_region': 2}, read_json, lambda _: 'C')


if __name__ == '__main__':
    unittest.main()
