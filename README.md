# Helper scripts for the influenza clade nomenclature repositories

These scripts are shared by the nomenclature repositories of the
[influenza-clade-nomenclature](https://github.com/influenza-clade-nomenclature) organization
(e.g. `seasonal_A-H3N2_HA`). They expect to be checked out next to these repositories and are
run from within a nomenclature repository as `python ../helper-scripts/<script>.py`.
Requirements: `biopython` and `pyyaml`; `compute_founders.py --all-sequences` additionally
requires `seqkit` and `nextclade` on the path.

| script | purpose |
|---|---|
| `construct_tsv.py` | tsv files, alias table and newick tree of the lineage hierarchy for augur/nextclade |
| `generate_markdown_summary.py` | markdown summary of the lineage definitions |
| `compute_founders.py` | founder sequences of new lineages from their representative sequences |
| `assemble_founders.py` | assemble and validate all founder sequences, list representatives |
| `seed_founders.py` | one-off: seed founders from ancestral sequences (e.g. of the nextclade dataset workflow) |
| `founder_utils.py` | shared functions for the founder scripts |

## Adding a new lineage

All commands are run from the root of the nomenclature repository.

### 1. Write the yaml file
Create `subclades/<name>.yml` (or `clades/<name>.yml`), e.g.
```
name: K.4
unaliased_name: A.2.1.1.2.2.1.2.1.3.1.1.2.4.1.4
parent: K
representatives:
- isolate: A/Wisconsin/82/2026
  source: gisaid
  accession: EPI1234567
- isolate: A/SouthAustralia/2527007683/2025
  source: gisaid
  accession: EPI4739151
defining_mutations:
- locus: HA1
  position: 158
  state: D
clade: none
```
- `parent` has to be an existing lineage. A lineage that receives a new letter alias additionally
  carries `alias_of: <full name>`, and the alias has to be added to `config/aliases.json`.
- Defining mutations use the coordinates of `config/reference.fasta` (`locus: nuc`) or of the
  proteins in `config/annotation.gff`. They need not be exhaustive, but every one of them has to be
  carried by the founder (see below).
- Choose at least two, better three or more, **representatives**. They determine the founder
  sequence, so they should
  - be (nearly) complete sequences of the segment,
  - sample the basal diversity of the lineage, i.e. come from different sublineages close to
    the root of the lineage rather than from one twig,
  - not carry the defining mutations of descendant lineages.

  Isolate names are matched ignoring spaces and case.

### 2. Compute the founder sequence
The founder is the founder of the parent plus all differences that are shared by every
representative covering a position. Positions where representatives disagree keep the state of
the parent, positions covered by fewer than `--min-representatives` (default 2) representatives
as well.

`compute_founders.py --all-sequences` takes a fasta containing the representatives among any number
of other sequences, e.g. `data/h3n2/ha.fasta` of a
[nextstrain/seasonal-flu](https://github.com/nextstrain/seasonal-flu) build. Headers have to start
with the strain name without spaces, optionally followed by `|`-separated fields (accessions are
matched against these fields as well). The script extracts the representatives with `seqkit grep`,
aligns them with `nextclade --gap-alignment-side right` against `config/reference.fasta`, and computes
the founder. Other segments of the same isolate do not align and are dropped automatically; if an
isolate occurs several times, the best covered sequence is used.

Check the report first with `--dry-run`:
```
python ../helper-scripts/compute_founders.py --input-dir subclades \
    --reference config/reference.fasta --annotation config/annotation.gff \
    --all-sequences ../../seasonal-flu/data/h3n2/ha.fasta --lineages K.4 --dry-run
```
The report lists for each lineage
- the number of founder mutations relative to the parent founder,
- representatives that were not found or have low coverage (`--min-coverage`, default 0.9),
- positions at which the representatives disagree,
- `INCONSISTENT` defining mutations that the founder does not carry.

Representatives that are missing, poorly covered, or disagree with the others at many positions
are worth replacing. If a defining mutation is not carried by the founder, either a representative
lacks it (replace it) or the defining mutation is wrong. Founders with inconsistent defining
mutations are not written unless `--force` is given.

If the report looks good, run the same command without `--dry-run`. This adds a `founder` block to
the yaml:
```
founder:
  method: representative-consensus
  provenance: "consensus of A/Wisconsin/82/2026, A/SouthAustralia/2527007683/2025"
  mutations:
  - {position: 473, state: 'T'}
```
Lineages that already have a founder are skipped unless `--overwrite` is given. When several
nested lineages are added at once, list them all in `--lineages`: parents are processed before
their children. Use `--output-sequences` to keep the aligned representatives, or
`--sequences aligned.fasta` instead of `--all-sequences` to provide representatives that are
already aligned to the reference.

### 3. Validate and update the auto-generated files
```
python ../helper-scripts/assemble_founders.py --input-dir subclades \
    --reference config/reference.fasta --annotation config/annotation.gff --strict \
    --output-fasta .auto-generated/founders.fasta --output-representatives .auto-generated/representatives.tsv
```
`--strict` fails if a founder is missing, lacks one of its defining mutations, is identical to the
founder of its parent, or places a deletion differently from the convention. Then regenerate the
remaining auto-generated files (markdown summary, tsv files, newick) with the commands listed in the
README of the nomenclature repository and commit the yaml together with the auto-generated files.
