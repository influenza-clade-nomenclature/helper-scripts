"""
Assemble founder sequences from the `founder` blocks of the lineage yamls,
validate them against the defining mutations, and write auto-generated outputs:
an aligned fasta of founders and a table of representative strains.
"""
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from founder_utils import (load_lineages, load_reference, load_annotation, assemble_founders,
                           check_defining_mutations, left_align_gaps, parent_of)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', nargs='+', required=True)
    parser.add_argument('--reference', required=True, help='reference fasta defining nucleotide coordinates')
    parser.add_argument('--annotation', required=True, help='gff3 defining the coordinates of the protein loci')
    parser.add_argument('--output-fasta', help='aligned founder sequences')
    parser.add_argument('--output-representatives', help='tsv of representative strains per lineage')
    parser.add_argument('--strict', action='store_true',
                        help='exit with an error if founders are missing or inconsistent with defining mutations')
    args = parser.parse_args()

    lineages = load_lineages(args.input_dir)
    reference = load_reference(args.reference)
    annotation = load_annotation(args.annotation)

    founders, missing = assemble_founders(lineages, reference)

    problems = []
    for name, seq in founders.items():
        problems.extend(check_defining_mutations(name, seq, lineages[name], annotation))
        p = parent_of(lineages, name)
        if p and founders.get(p) == seq:
            problems.append(f"{name}: founder identical to founder of parent {p}")
        if left_align_gaps(seq, reference) != seq:
            problems.append(f"{name}: deletion not placed as far 5' as the reference allows, "
                            f"align representatives with nextclade --gap-alignment-side right")

    for name in missing:
        print(f"MISSING founder: {name}", file=sys.stderr)
    for p in problems:
        print(f"INCONSISTENT {p}", file=sys.stderr)
    print(f"{len(founders)} founders assembled, {len(missing)} missing, {len(problems)} problems", file=sys.stderr)

    if args.output_fasta:
        SeqIO.write([SeqRecord(Seq(founders[n]), id=n, description='') for n in sorted(founders)],
                    args.output_fasta, 'fasta')

    if args.output_representatives:
        with open(args.output_representatives, 'w') as fh:
            fh.write('\t'.join(['lineage', 'isolate', 'source', 'accession']) + '\n')
            for name in sorted(lineages):
                for r in lineages[name].get('representatives') or []:
                    fh.write('\t'.join([name, r.get('isolate', ''), r.get('source', ''), str(r.get('accession', ''))]) + '\n')

    if args.strict and (missing or problems):
        sys.exit(1)
