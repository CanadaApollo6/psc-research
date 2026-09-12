"""Exact sequence identity and missingness helpers for the frozen follow-up.

Coordinates are one based. Unordered GWAS alleles are not reference alleles.
An rsID match is never used to establish the identity of an edit.
"""

from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from coloc_common import reciprocal_base


@dataclass(frozen=True)
class Reference:
    chromosome: str
    start: int
    sequence: str

    def get(self, position, length=1):
        offset = position - self.start
        if offset < 0 or length < 1 or offset + length > len(self.sequence):
            raise ValueError('Edit or normalization reaches reference block boundary')
        return self.sequence[offset:offset + length]

    def haplotype(self, position, ref, alt):
        if self.get(position, len(ref)) != ref:
            raise ValueError('REF does not match reference sequence')
        offset = position - self.start
        return self.sequence[:offset] + alt + self.sequence[offset + len(ref):]


def normalize(reference, position, ref, alt):
    """Left align and minimally represent a nonsymbolic biallelic VCF edit.

    Trim common suffixes, extending to the left before an allele becomes
    empty. This handles repeat rotations as well as ordinary padding.
    Independently cross-checked with bcftools norm for the real targets.
    """
    if not ref or not alt or ref == alt or set(ref + alt) - set('ACGT'):
        raise ValueError('Expected two distinct explicit DNA alleles')
    original = reference.haplotype(position, ref, alt)
    while ref[-1] == alt[-1]:
        if min(len(ref), len(alt)) == 1:
            previous = reference.get(position - 1)
            position -= 1
            ref, alt = previous + ref, previous + alt
        ref, alt = ref[:-1], alt[:-1]
    while min(len(ref), len(alt)) > 1 and ref[0] == alt[0]:
        position += 1
        ref, alt = ref[1:], alt[1:]
    if original != reference.haplotype(position, ref, alt):
        raise ValueError('Normalization changed alternate haplotype')
    return position, ref, alt


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def reciprocal_edit(reference, destination, forward, reverse, position, ref, alt, flank=100):
    """Require every REF and flank base to map uniquely and reciprocally."""
    source_start, length = position - flank, len(ref) + 2 * flank
    sequence = reference.get(source_start, length)
    hits = [reciprocal_base(forward, reverse, reference.chromosome, p - 1)
            for p in range(source_start, source_start + length)]
    if any(h is None for h in hits):
        raise ValueError('Nonunique or nonreciprocal full-allele mapping')
    strand = hits[0][2]
    step = 1 if strand == '+' else -1
    if any(h[0] != destination.chromosome or h[2] != strand or h[1] != hits[0][1] + step * i
           for i, h in enumerate(hits)):
        raise ValueError('Noncontiguous allele or flank mapping')
    dest_start = min(h[1] for h in hits) + 1
    expected = sequence if strand == '+' else reverse_complement(sequence)
    if destination.get(dest_start, length) != expected:
        raise ValueError('Mapped reference flanks differ between builds')
    mapped_position = min(h[1] for h in hits[flank:flank + len(ref)]) + 1
    mapped_ref, mapped_alt = (ref, alt) if strand == '+' else (reverse_complement(ref), reverse_complement(alt))
    # Both source and destination edits must have a valid REF, not just an anchor.
    reference.haplotype(position, ref, alt)
    destination.haplotype(mapped_position, mapped_ref, mapped_alt)
    return (mapped_position, mapped_ref, mapped_alt), length, strand


def unordered_edit_candidates(reference, position, allele0, allele1):
    """Enumerate valid REF orientations, retaining ambiguity for the caller.

    A shorter repeat allele can also match a reference prefix. Its insertion
    interpretation is retained, never silently selected as the true GWAS REF.
    A supplied exact deletion haplotype can disambiguate the interpretations.
    """
    candidates = []
    for ref, alt in [(allele0, allele1), (allele1, allele0)]:
        if set(ref + alt) - set('ACGT') or not ref or not alt or ref == alt:
            continue
        if reference.get(position, len(ref)) == ref:
            candidates.append(normalize(reference, position, ref, alt))
    return sorted(set(candidates))


def exact_qtl_status(rows, position, ref, alt, gene):
    at_position = [r for r in rows if int(r['position']) == position]
    exact = [r for r in at_position if r['ref'] == ref and r['alt'] == alt]
    measured = [r for r in exact if r['gene_id'].split('.')[0] == gene]
    if measured:
        status = 'exact_allele_measured_for_target_gene'
    elif exact:
        status = 'exact_allele_measured_only_for_other_genes'
    elif at_position:
        status = 'position_measured_only_with_different_alleles'
    else:
        status = 'position_absent_across_all_genes'
    return status, at_position, exact, measured


def component_weights(logbf):
    values = np.asarray(logbf, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('Missing or infinite component evidence')
    return np.exp(values - logsumexp(values, axis=0))


def coverage_parts(weights, baseline, recovered):
    weights = np.asarray(weights)
    baseline, recovered = np.asarray(baseline, dtype=bool), np.asarray(recovered, dtype=bool)
    if weights.ndim != 1 or not (weights.shape == baseline.shape == recovered.shape):
        raise ValueError('Coverage masks do not match the full signal universe')
    if not np.isfinite(weights).all() or (weights < 0).any() or not np.isclose(weights.sum(), 1):
        raise ValueError('Use full normalized signal weights')
    if (baseline & recovered).any():
        raise ValueError('Recovered evidence was already counted')
    return float(weights[baseline].sum()), float(weights[recovered].sum()), float(weights[baseline | recovered].sum())
