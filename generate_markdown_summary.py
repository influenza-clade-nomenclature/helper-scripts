#! python
#%%
import yaml
import glob
from urllib.parse import quote

def generate_lineage_md(subclade, lineage, segment):
    lines = []
    revoked = subclade.get('revoked', False)
    if revoked:
        lines.append(f"## ~~{subclade['name']}~~ (revoked)")
    else:
        lines.append(f"## {subclade['name']}")

    if 'alias_of' in subclade:
        lines.append(f" * alias of: {subclade['alias_of']}")
    lines.append(f" * parent: [{subclade['parent']}](#{subclade['parent'].replace('.', '')})")
    if 'comment' in subclade and subclade['comment']:
        lines.append(f" * comment: {subclade['comment']}")
    snp_str = ', '.join(f"{x['locus']}:{x['position']}{x['state']}" for x in subclade['defining_mutations'])
    lines.append(f" * defining mutations or substitutions: {snp_str}")
    if "clade" in subclade and subclade['clade'] != "none":
        lines.append(f" * clade: {subclade['clade']}")

    ref_seqs = []
    for x in subclade['representatives']:
        nextstrain_link = f"[View on Nextstrain](https://nextstrain.org/seasonal-flu/{lineage}/{segment}/6y?c=subclade&s={quote(x['isolate'])})"
        if x['source']=='genbank' and 'accession' in x:
            accession_link = f"[{x['accession']}](https://www.ncbi.nlm.nih.gov/nuccore/{x['accession']})"
        elif x['source']=='gisaid':
            accession_link = x['accession']

        if 'other_accession' in x:
            other_accession = f", {x['other_accession']}"
        else: other_accession=''
        ref_seqs.append(f"{x.get('isolate', x['accession'])} ({accession_link}{other_accession}) {nextstrain_link}")

    if len(ref_seqs)==1:
        lines.append(f" * representative sequence: {ref_seqs[0]}")
    elif len(ref_seqs)>1:
        lines.append(" * representative sequences:")
        for r in ref_seqs:
            lines.append(f"   - {r}")
    lines.append(f" * [View on Nextstrain](https://nextstrain.org/seasonal-flu/{lineage}/{segment}/6y?branchLabel=Subclade&c=subclade&label=Subclade:{subclade['name']})")
    return '\n'.join(lines) + '\n'

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', required=True)
    parser.add_argument('--lineage')
    parser.add_argument('--segment')
    args = parser.parse_args()
    repo_name = f"{args.lineage}_{args.segment}"

    subclades = []
    # Iterate through all lineage definition files
    for yaml_file in sorted(glob.glob(f"{args.input_dir}/*.yml")):
        with open(yaml_file, 'r') as stream:
            yaml_data = yaml.safe_load(stream)
        subclades.append(yaml_data)

    subclades.sort(key=lambda x:x['name'])
    clade_lineage_map = [(x['clade'], x['name'], x['unaliased_name'])
                         for x in subclades if 'clade' in x and x['clade'] != 'none' and (not x.get('revoked', False))]
    # Write to json file
    with open('.auto-generated/subclades.md', 'w') as outfile:
        outfile.write("# Summary of designated subclades\n")

        for subclade in subclades:
            print("output clade", subclade['name'])
            outfile.write(generate_lineage_md(subclade, args.lineage, args.segment) + '\n')

        if len(clade_lineage_map):
            # write table of clade -- subclade correspondence
            outfile.write("# Clade -- subclade correspondence\n")
            outfile.write("|*Clade*|*Subclade*|*full subclade name*|\n")
            outfile.write("|-------------|---------|----------------------|\n")
            for clade, lineage, unaliased_name in clade_lineage_map:
                outfile.write(f"|{clade}|[{lineage}](#{lineage.replace('.','')})|{unaliased_name}|\n")


    with open('subclades.tex', 'w') as latexoutfile:
        latexoutfile.write("Subclade & Clade & full subclade name\\\\\\hline\n")

        for subclade in subclades:
            latexoutfile.write(f"{subclade['name']} & {subclade.get('clade','')} & {subclade['unaliased_name']}\\\\\n")

    with open('subclades.tsv', 'w') as tsvoutfile:
        tsvoutfile.write("Subclade\tClade\tfull subclade name\n")

        for subclade in subclades:
            tsvoutfile.write(f"{subclade['name']}\t{subclade.get('clade','')}\t{subclade['unaliased_name']}\n")

