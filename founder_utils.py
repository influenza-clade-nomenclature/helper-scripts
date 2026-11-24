"""
Shared helpers for lineage founder sequences.

Each lineage yaml can carry a `founder` block that describes the founder sequence
of the lineage as nucleotide differences relative to the founder of its parent
(or relative to the reference sequence for root lineages):

    founder:
      method: representative-consensus
      provenance: "..."
      mutations:
      - {position: 473, state: T}

Positions are 1-based in the coordinates of the repository reference sequence,
same as the `nuc` locus in `defining_mutations`. Deletions are encoded as state `-`.
"""
import glob
import itertools
import re

import yaml
from Bio import SeqIO
from Bio.Seq import Seq

NUCS = set('ACGT-')


def load_lineages(input_dirs, include_revoked=False):
    """Read all lineage yamls. Returns dict name -> yaml data (with '_file' added)."""
    lineages = {}
    for d in input_dirs:
        for fname in sorted(glob.glob(d + '/*.yml')):
            with open(fname) as fh:
                data = yaml.safe_load(fh)
            if data.get('revoked', False) and not include_revoked:
                continue
            data['_file'] = fname
            lineages[data['name']] = data
    return lineages


def parent_of(lineages, name):
    p = lineages[name].get('parent', 'none')
    return None if p in (None, 'none', '') else p


def preorder(lineages):
    """Lineage names ordered such that parents precede their children."""
    children = {}
    roots = []
    for name in lineages:
        p = parent_of(lineages, name)
        if p is None or p not in lineages:
            roots.append(name)
        else:
            children.setdefault(p, []).append(name)
    order, stack = [], sorted(roots, reverse=True)
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(sorted(children.get(n, []), reverse=True))
    return order


def load_reference(fasta):
    rec = next(SeqIO.parse(fasta, 'fasta'))
    return str(rec.seq).upper()


def load_annotation(gff):
    """Parse a minimal gff3. Returns dict feature name -> (start, end), 1-based inclusive."""
    features = {}
    for line in open(gff):
        if line.startswith('#') or not line.strip():
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9:
            # gff blocks copied from documentation are often separated by spaces
            fields = line.split(None, 8)
        if len(fields) < 9:
            raise ValueError(f"Cannot parse annotation line: {line.rstrip()}")
        m = re.search(r'(?:^|;)\s*(?:gene_name|name|Name)="?([^";]+)"?', fields[8])
        if m:
            features[m.group(1)] = (int(fields[3]), int(fields[4]))
    return features


def apply_mutations(seq, mutations):
    seq = list(seq)
    for m in mutations:
        seq[m['position'] - 1] = m['state']
    return ''.join(seq)


def diff(seq, background):
    return [{'position': i + 1, 'state': b}
            for i, (a, b) in enumerate(zip(background, seq)) if a != b]


def assemble_founders(lineages, reference):
    """
    Build founder sequences for all lineages with a founder block.
    Lineages whose founder (or that of an ancestor) is missing are skipped.
    Returns dict name -> sequence and list of names lacking founders.
    """
    founders, missing = {}, []
    for name in preorder(lineages):
        block = lineages[name].get('founder')
        p = parent_of(lineages, name)
        background = reference if p is None else founders.get(p)
        if block is None or background is None:
            missing.append(name)
            continue
        founders[name] = apply_mutations(background, block.get('mutations') or [])
    return founders, missing


def translate(seq, start, end):
    codons = [seq[i:i + 3] for i in range(start - 1, end - 2, 3)]
    aa = []
    for c in codons:
        if c == '---':
            aa.append('-')
        elif set(c) <= set('ACGT'):
            aa.append(str(Seq(c).translate()))
        else:
            aa.append('X')
    return ''.join(aa)


def left_align_gaps(seq, reference):
    """
    Shift internal deletions as far towards the 5' end as the reference allows.

    A deletion in a repeat can be placed at several equivalent positions, and
    aligners do not all make the same choice. Left alignment is the convention
    that matches how deletions are usually named, e.g. the deletion of B/Vic HA1
    codons 162-164 is placed at nucleotides 562-570 rather than further 3'.
    Gap runs at either end are missing coverage, not deletions, and are left alone.
    """
    s = list(seq)
    i = 0
    while i < len(s):
        if s[i] != '-':
            i += 1
            continue
        j = i
        while j < len(s) and s[j] == '-':
            j += 1
        if i > 0 and j < len(s):  # internal gap run [i, j)
            while (i > 0 and s[i - 1] != '-'
                   and s[i - 1] == reference[i - 1] == reference[j - 1]):
                s[j - 1], s[i - 1] = s[i - 1], '-'
                i -= 1
                j -= 1
        i = j
    return ''.join(s)


def children_of(lineages, name):
    return [n for n in lineages if parent_of(lineages, n) == name]


def descendants_of(lineages, name):
    out, stack = [], [name]
    while stack:
        for child in children_of(lineages, stack.pop()):
            out.append(child)
            stack.append(child)
    return out


def _codon_for_aa(current_codon, target_aa, preferred=None):
    """Pick a codon encoding target_aa, preferring `preferred`, else the fewest changes."""
    candidates = [''.join(c) for c in itertools.product('ACGT', repeat=3)
                  if str(Seq(''.join(c)).translate()) == target_aa]
    if not candidates:
        return None
    if preferred in candidates:
        return preferred
    return min(candidates, key=lambda c: (sum(a != b for a, b in zip(c, current_codon)), c))


def enforce_defining_mutations(lineages, founders, annotation):
    """
    Make every founder carry its own defining mutations.

    Ancestral reconstruction sometimes places a defining mutation one node below the
    lineage root, leaving the founder without it. For amino acid loci the codon used by
    the lineage's descendants is preferred, otherwise the codon closest to the current one.
    Returns (repairs, unfixable) as lists of human readable strings.
    """
    repairs, unfixable = [], []
    for name in preorder(lineages):
        if name not in founders:
            continue
        for m in lineages[name].get('defining_mutations') or []:
            locus, pos, state = m['locus'], m['position'], str(m['state'])
            seq = founders[name]
            if locus == 'nuc':
                if seq[pos - 1] == state:
                    continue
                founders[name] = seq[:pos - 1] + state + seq[pos:]
                repairs.append(f"{name}: nuc:{pos} {seq[pos - 1]} -> {state}")
                continue
            if locus not in annotation:
                unfixable.append(f"{name}: unknown locus {locus}")
                continue
            start, end = annotation[locus]
            if translate(seq, start, end)[pos - 1] == state:
                continue
            i = start - 1 + (pos - 1) * 3
            current = seq[i:i + 3]
            kids = [founders[d] for d in descendants_of(lineages, name) if d in founders]
            kid_codons = {k[i:i + 3] for k in kids}
            preferred = kid_codons.pop() if len(kid_codons) == 1 else None
            codon = _codon_for_aa(current, state, preferred)
            if codon is None:
                unfixable.append(f"{name}: no codon encodes {locus}:{pos} {state}")
                continue
            founders[name] = seq[:i] + codon + seq[i + 3:]
            repairs.append(f"{name}: {locus}:{pos} {current} -> {codon} "
                           f"({translate(seq, start, end)[pos - 1]} -> {state})"
                           + (" [codon of descendants]" if codon == preferred else ""))
    return repairs, unfixable


def check_defining_mutations(name, seq, lineage, annotation):
    """Return list of defining mutations not carried by the founder sequence."""
    problems = []
    for m in lineage.get('defining_mutations') or []:
        locus, pos, state = m['locus'], m['position'], str(m['state'])
        if locus == 'nuc':
            observed = seq[pos - 1]
        elif locus in annotation:
            observed = translate(seq, *annotation[locus])[pos - 1]
        else:
            problems.append(f"{name}: unknown locus {locus}")
            continue
        if observed != state:
            problems.append(f"{name}: {locus}:{pos} expected {state}, founder has {observed}")
    return problems


def write_founder_block(yaml_file, method, mutations, provenance=None):
    """Replace (or append) the `founder` block of a lineage yaml, leaving the rest untouched."""
    lines = open(yaml_file).read().splitlines()
    out, skipping = [], False
    for line in lines:
        if line.startswith('founder:'):
            skipping = True
            continue
        if skipping and (line.startswith(' ') or line.startswith('-')):
            continue
        skipping = False
        out.append(line)
    out.append('founder:')
    out.append(f'  method: {method}')
    if provenance:
        out.append(f'  provenance: "{provenance}"')
    if mutations:
        out.append('  mutations:')
        out.extend(f"  - {{position: {m['position']}, state: '{m['state']}'}}" for m in mutations)
    else:
        out.append('  mutations: []')
    with open(yaml_file, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
