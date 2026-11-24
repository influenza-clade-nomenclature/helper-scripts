"""
Compute lineage founders from representative sequences.

At each position, the founder takes the state shared by all representatives (with
an unambiguous state at that position) if at least `--min-representatives` of them
are informative and their shared state differs from the parent founder. Otherwise
the founder inherits the state of the parent founder. Requiring unanimity rather
than a majority keeps mutations that arose within the lineage out of the founder,
provided the representatives sample the basal diversity of the lineage.

Representative sequences have to be aligned to the repository reference, e.g.
    nextclade run --gap-alignment-side right --input-ref config/reference.fasta \
        --input-annotation config/annotation.gff --output-fasta aligned.fasta reps.fasta
The gap alignment side matters: nextclade defaults to `left`, which places a deletion in a
repeat further 3' than the convention used for the founders.
They are matched to the yaml entries by isolate name (ignoring spaces and case)
or by accession.
"""
import sys

from Bio import SeqIO

from founder_utils import (load_lineages, load_reference, load_annotation, assemble_founders,
                           preorder, parent_of, diff, apply_mutations, check_defining_mutations,
                           left_align_gaps, write_founder_block)


def normalize(name):
    return name.replace(' ', '').lower()


def consensus_founder(rep_seqs, parent_founder, min_reps):
    founder = list(parent_founder)
    unresolved = []
    ranges = [aligned_range(s) for s in rep_seqs]
    for i, parent_state in enumerate(parent_founder):
        covered = [s for s, (lo, hi) in zip(rep_seqs, ranges) if lo <= i <= hi]
        # internal gaps can be real deletions, but only adopt one if no representative has a
        # base here; gaps outside a sequence's aligned range are trimming and never informative
        states = [s[i] for s in covered if s[i] in 'ACGT'] or [s[i] for s in covered if s[i] == '-']
        if len(states) < min_reps:
            continue
        distinct = set(states)
        if len(distinct) == 1:
            founder[i] = states[0]
        else:
            # representatives disagree, so this mutation arose within the lineage as far as
            # we can tell: keep the parent state, but report it as it may indicate that the
            # representatives do not sample the basal diversity of the lineage
            unresolved.append((i + 1, ''.join(sorted(distinct)), parent_state))
    return ''.join(founder), unresolved


def aligned_range(seq):
    """First and last index of a non-gap character; positions outside are untrimmed ends."""
    bases = [i for i, c in enumerate(seq) if c != '-']
    return (bases[0], bases[-1]) if bases else (0, -1)


def coverage(seq):
    return sum(c in 'ACGT' for c in seq) / len(seq)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input-dir', nargs='+', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--annotation', required=True)
    parser.add_argument('--sequences', required=True, help='representative sequences aligned to the reference')
    parser.add_argument('--lineages', nargs='+', help='lineages to compute (default: all without founder)')
    parser.add_argument('--min-representatives', type=int, default=2)
    parser.add_argument('--overwrite', action='store_true', help='replace existing founder blocks')
    parser.add_argument('--min-coverage', type=float, default=0.9,
                        help='warn about representatives with a smaller fraction of unambiguous bases')
    parser.add_argument('--force', action='store_true',
                        help='write founders even if they lack defining mutations')
    parser.add_argument('--dry-run', action='store_true', help='report but do not modify yamls')
    args = parser.parse_args()

    lineages = load_lineages(args.input_dir)
    reference = load_reference(args.reference)
    annotation = load_annotation(args.annotation)
    founders, _ = assemble_founders(lineages, reference)

    sequences, shifted = {}, []
    for r in SeqIO.parse(args.sequences, 'fasta'):
        seq = str(r.seq).upper()
        if len(seq) != len(reference):
            print(f"Skipping {r.id}: not aligned to reference", file=sys.stderr)
            continue
        aligned = left_align_gaps(seq, reference)
        if aligned != seq:
            shifted.append(r.id)
        sequences[normalize(r.id)] = aligned

    if shifted:
        print(f"WARNING {len(shifted)} of {len(sequences)} representatives had a deletion placed "
              f"differently from the founders and were adjusted. Align them with\n"
              f"  nextclade run --gap-alignment-side right --input-ref ... --input-annotation ...\n"
              f"to place deletions in repeats the same way, e.g. B/Vic HA1 162-164 at nucleotides 562-570.",
              file=sys.stderr)

    targets = set(args.lineages) if args.lineages else {n for n in lineages if 'founder' not in lineages[n]}
    for name in preorder(lineages):
        if name not in targets:
            continue
        if 'founder' in lineages[name] and not args.overwrite:
            print(f"{name}: founder exists, use --overwrite to replace", file=sys.stderr)
            continue
        p = parent_of(lineages, name)
        background = reference if p is None else founders.get(p)
        if background is None:
            print(f"{name}: SKIPPED, parent {p} has no founder", file=sys.stderr)
            continue

        reps = lineages[name].get('representatives') or []
        rep_seqs, found = [], []
        for r in reps:
            for key in (r.get('isolate'), r.get('accession')):
                if key and normalize(str(key)) in sequences:
                    seq = sequences[normalize(str(key))]
                    rep_seqs.append(seq)
                    found.append(r.get('isolate') or str(r.get('accession')))
                    if coverage(seq) < args.min_coverage:
                        print(f"  WARNING {name}: representative {found[-1]} has coverage {coverage(seq):.2f}",
                              file=sys.stderr)
                    break
        if len(rep_seqs) < args.min_representatives:
            print(f"{name}: SKIPPED, {len(rep_seqs)} of {len(reps)} representatives found "
                  f"(need {args.min_representatives})", file=sys.stderr)
            continue

        founder, unresolved = consensus_founder(rep_seqs, background, args.min_representatives)
        mutations = diff(founder, background)
        problems = check_defining_mutations(name, founder, lineages[name], annotation)
        print(f"{name}: {len(mutations)} mutations relative to {p or 'reference'} "
              f"from {len(rep_seqs)} representatives", file=sys.stderr)
        for pos, states, kept in unresolved:
            print(f"  position {pos}: representatives disagree ({states}), kept parent state {kept}",
                  file=sys.stderr)
        for problem in problems:
            print(f"  INCONSISTENT {problem}", file=sys.stderr)

        if problems and not args.force:
            print(f"  {name}: NOT WRITTEN, founder lacks defining mutations (use --force to override)",
                  file=sys.stderr)
            continue
        if not args.dry_run:
            write_founder_block(lineages[name]['_file'], 'representative-consensus', mutations,
                                'consensus of ' + ', '.join(found))
        # make the new founder available to descendants processed later in this run
        founders[name] = apply_mutations(background, mutations)
