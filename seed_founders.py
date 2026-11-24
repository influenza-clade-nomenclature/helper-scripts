"""
Seed `founder` blocks of lineage yamls from an existing set of full-length founder
sequences (e.g. ancestral sequences inferred by the nextclade dataset workflow).
The sequences need to be aligned to the repository reference and named by lineage.
Each lineage's founder is stored as the difference to the founder of its parent.
"""
import sys

from Bio import SeqIO

from founder_utils import (load_lineages, load_reference, load_annotation, preorder, parent_of, diff,
                           enforce_defining_mutations, write_founder_block)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', nargs='+', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--founders', required=True, help='aligned fasta with sequences named by lineage')
    parser.add_argument('--annotation', help='gff3, required for --enforce-defining-mutations')
    parser.add_argument('--enforce-defining-mutations', action='store_true',
                        help='repair founders that do not carry their own defining mutations')
    parser.add_argument('--method', default='ancestral-inference')
    parser.add_argument('--provenance', default=None)
    parser.add_argument('--overwrite', action='store_true', help='replace existing founder blocks')
    args = parser.parse_args()

    lineages = load_lineages(args.input_dir)
    reference = load_reference(args.reference)
    founders = {r.id: str(r.seq).upper() for r in SeqIO.parse(args.founders, 'fasta')}
    for name, seq in founders.items():
        if len(seq) != len(reference):
            sys.exit(f"Founder {name} has length {len(seq)}, reference has length {len(reference)}")

    if args.enforce_defining_mutations:
        if not args.annotation:
            sys.exit("--enforce-defining-mutations requires --annotation")
        repairs, unfixable = enforce_defining_mutations(lineages, founders, load_annotation(args.annotation))
        for r in repairs:
            print(f"REPAIRED {r}", file=sys.stderr)
        for u in unfixable:
            print(f"UNFIXABLE {u}", file=sys.stderr)

    for name in preorder(lineages):
        if name not in founders:
            print(f"No founder sequence for {name}", file=sys.stderr)
            continue
        if 'founder' in lineages[name] and not args.overwrite:
            continue
        p = parent_of(lineages, name)
        if p is None:
            background = reference
        elif p in founders:
            background = founders[p]
        else:
            print(f"Skipping {name}: no founder sequence for parent {p}", file=sys.stderr)
            continue
        # ambiguous states in the seed are treated as identical to the background
        mutations = [m for m in diff(founders[name], background) if m['state'] in 'ACGT-']
        write_founder_block(lineages[name]['_file'], args.method, mutations, args.provenance)
        print(f"{name}: {len(mutations)} mutations relative to {p or 'reference'}", file=sys.stderr)
