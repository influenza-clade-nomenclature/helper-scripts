"""
Compute lineage founders from representative sequences.

At each position, the founder takes the state shared by all representatives (with
an unambiguous state at that position) if at least `--min-representatives` of them
are informative and their shared state differs from the parent founder. Otherwise
the founder inherits the state of the parent founder. Requiring unanimity rather
than a majority keeps mutations that arose within the lineage out of the founder,
provided the representatives sample the basal diversity of the lineage.

Representative sequences can be provided in two ways:
  --all-sequences: any fasta (e.g. a GISAID download) with strain names without spaces at
      the start of the header, optionally followed by `|`-separated fields. Representatives
      of the target lineages are extracted with seqkit and aligned with nextclade.
  --sequences: representatives already aligned to the repository reference, e.g.
      nextclade run --gap-alignment-side right --input-ref config/reference.fasta \
          --input-annotation config/annotation.gff --output-fasta aligned.fasta reps.fasta
The gap alignment side matters: nextclade defaults to `left`, which places a deletion in a
repeat further 3' than the convention used for the founders.
Sequences are matched to the yaml entries by isolate name (ignoring spaces and case)
or by accession, using any `|`-separated field of the header.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

from Bio import SeqIO

from founder_utils import (load_lineages, load_reference, load_annotation, assemble_founders,
                           preorder, parent_of, diff, apply_mutations, check_defining_mutations,
                           left_align_gaps, write_founder_block)


def normalize(name):
    return name.replace(' ', '').lower()


def header_keys(header):
    """Normalized names a sequence can be matched by: each `|`-separated field, without accession version."""
    keys = set()
    for field in header.split('|'):
        key = normalize(field)
        if key:
            keys.add(key)
            keys.add(re.sub(r'\.\d+$', '', key))
    return keys


def go_regex_escape(text):
    # seqkit uses Go regular expressions, which reject some escapes produced by re.escape
    return re.sub(r'([\\.+*?()|\[\]{}^$])', r'\\\1', text)


def extract_and_align(all_sequences, keys, reference_fasta, annotation_gff, workdir):
    """Extract sequences matching isolate names or accessions with seqkit and align them with nextclade."""
    for tool in ('seqkit', 'nextclade'):
        if shutil.which(tool) is None:
            sys.exit(f"{tool} not found, it is required for --all-sequences")
    patterns = os.path.join(workdir, 'patterns.txt')
    with open(patterns, 'w') as fh:
        for key in sorted(keys):
            fh.write(f"(^|\\|){go_regex_escape(key)}(\\.\\d+)?(\\||$)\n")
    extracted = os.path.join(workdir, 'representatives.fasta')
    aligned = os.path.join(workdir, 'representatives.aligned.fasta')
    # without spaces in strain names, the sequence ID seqkit matches against is the full header
    subprocess.run(['seqkit', 'grep', '--quiet', '-r', '-i', '-f', patterns, '-o', extracted, all_sequences],
                   check=True)
    n_extracted = sum(1 for _ in SeqIO.parse(extracted, 'fasta'))
    print(f"Extracted {n_extracted} sequences matching {len(keys)} isolate names and accessions",
          file=sys.stderr)
    if n_extracted == 0:
        return aligned, n_extracted
    # sequences that do not align, e.g. other segments of the same isolate, are dropped by nextclade
    subprocess.run(['nextclade', 'run', '--quiet', '--gap-alignment-side', 'right',
                    '--input-ref', reference_fasta, '--input-annotation', annotation_gff,
                    '--output-fasta', aligned, extracted], check=True)
    return aligned, n_extracted


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
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--all-sequences', help='unaligned fasta containing the representatives among other sequences')
    source.add_argument('--sequences', help='representative sequences aligned to the reference')
    parser.add_argument('--output-sequences', help='with --all-sequences: save the aligned representatives here')
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

    targets = set(args.lineages) if args.lineages else {n for n in lineages if 'founder' not in lineages[n]}
    unknown = targets - set(lineages)
    if unknown:
        sys.exit(f"Unknown lineages: {', '.join(sorted(unknown))}")

    with tempfile.TemporaryDirectory() as workdir:
        if args.all_sequences:
            keys = {normalize(str(r[k])) for n in targets for r in lineages[n].get('representatives') or []
                    for k in ('isolate', 'accession') if r.get(k)}
            aligned_fasta, n_extracted = extract_and_align(args.all_sequences, keys, args.reference,
                                                           args.annotation, workdir)
            if args.output_sequences and n_extracted:
                shutil.copy(aligned_fasta, args.output_sequences)
            records = list(SeqIO.parse(aligned_fasta, 'fasta')) if n_extracted else []
        else:
            records = list(SeqIO.parse(args.sequences, 'fasta'))

    sequences, shifted, duplicates = {}, [], set()
    for r in records:
        seq = str(r.seq).upper()
        if len(seq) != len(reference):
            print(f"Skipping {r.id}: not aligned to reference", file=sys.stderr)
            continue
        aligned = left_align_gaps(seq, reference)
        if aligned != seq:
            shifted.append(r.id)
        for key in header_keys(r.id):
            # several sequences of the same isolate (e.g. different passages): keep the best covered
            if key in sequences:
                duplicates.add(key)
                if coverage(sequences[key]) >= coverage(aligned):
                    continue
            sequences[key] = aligned

    if shifted:
        print(f"WARNING {len(shifted)} of {len(records)} representatives had a deletion placed "
              f"differently from the founders and were adjusted. Align them with\n"
              f"  nextclade run --gap-alignment-side right --input-ref ... --input-annotation ...\n"
              f"to place deletions in repeats the same way, e.g. B/Vic HA1 162-164 at nucleotides 562-570.",
              file=sys.stderr)

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
            label = r.get('isolate') or str(r.get('accession'))
            for key in (r.get('isolate'), r.get('accession')):
                if key and normalize(str(key)) in sequences:
                    seq = sequences[normalize(str(key))]
                    rep_seqs.append(seq)
                    found.append(label)
                    if normalize(str(key)) in duplicates:
                        print(f"  NOTE {name}: several sequences match representative {label}, "
                              f"using the best covered", file=sys.stderr)
                    if coverage(seq) < args.min_coverage:
                        print(f"  WARNING {name}: representative {label} has coverage {coverage(seq):.2f}",
                              file=sys.stderr)
                    break
            else:
                print(f"  WARNING {name}: representative {label} not found", file=sys.stderr)
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
