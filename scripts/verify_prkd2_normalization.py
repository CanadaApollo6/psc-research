"""Cross-check every queued full-allele normalization with independent bcftools."""

import json
import subprocess
from pathlib import Path

import pysam
import pysam.bcftools

from prkd2_common import ROOT, save_json, sha
from prkd2_expanded_common import load_references


def main():
    checks = json.loads((ROOT / 'work/prkd2/expanded-normalization-checks.json').read_text())
    refs = load_references()
    folder = ROOT / 'work/prkd2/normalization'
    folder.mkdir(exist_ok=True)
    counts = {}
    for build, reference in refs.items():
        selected = [r for r in checks if r['build'] == build]
        expected = {r['id']: tuple(r['expected']) for r in selected}
        assert len(expected) == len(selected)
        fasta = folder / (build + '.fa')
        fasta.write_text('>block\n' + reference.sequence + '\n')
        pysam.faidx(str(fasta))
        vcf = folder / (build + '.vcf')
        lines = ['##fileformat=VCFv4.2', f'##contig=<ID=block,length={len(reference.sequence)}>',
                 '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']
        for r in selected:
            lines.append(f"block\t{r['position'] - reference.start + 1}\t{r['id']}\t{r['ref']}\t{r['alt']}\t.\tPASS\t.")
        vcf.write_text('\n'.join(lines) + '\n')
        result = pysam.bcftools.norm('-f', str(fasta), '-c', 'e', '-Ov', str(vcf), catch_stdout=True)
        observed = {}
        for line in result.splitlines():
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            value = (int(fields[1]) + reference.start - 1, fields[3], fields[4])
            if fields[2] in observed:
                raise ValueError('Duplicate normalization output')
            observed[fields[2]] = value
        if observed != expected:
            failures = [k for k in expected if observed.get(k) != expected[k]]
            raise ValueError('Independent normalization differs: ' + str(failures[:10]))
        counts[build] = len(selected)
    sources = json.loads((ROOT / 'config/prkd2-expanded-sources.json').read_text())
    previous_reference = ROOT / 'work/ubash3a-raw-read-audit/reference.json'
    old = json.loads(previous_reference.read_text())
    assert sources['local_GRCh38_reference']['sha256'] == old['fasta_sha256']
    name = 'Homo_sapiens.GRCh37.dna.chromosome.19.fa.gz'
    published = next(line.split()[:2] for line in (ROOT / 'data/raw/prkd2-CHECKSUMS').read_text().splitlines() if line.split()[-1] == name)
    calculated = subprocess.check_output(['sum', str(ROOT / ('data/raw/prkd2-' + name))], text=True).split()[:2]
    assert published == calculated
    record = {'all_checks_pass': True, 'independent_normalization_checks': counts, 'total_normalizations': sum(counts.values()),
        'pysam': pysam.__version__, 'bcftools': pysam.__samtools_version__,
        'GRCh37_publisher_BSD_checksum_and_blocks': published,
        'GRCh38_prior_verified_decompression_matches': True,
        'GRCh38_prior_reference_record_sha256': sha(previous_reference.read_bytes()),
        'check_input_sha256': sha((ROOT / 'work/prkd2/expanded-normalization-checks.json').read_bytes()),
        'expanded_source_manifest_sha256': sha((ROOT / 'config/prkd2-expanded-sources.json').read_bytes()),
        'verifier_sha256': sha(Path(__file__).read_bytes())}
    save_json(ROOT / 'reports/prkd2-normalization-verification.json', record)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
