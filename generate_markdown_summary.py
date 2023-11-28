#! python
#%%
import yaml
import glob

def generate_lineage_md(subclade, lineage, segment):
    lines = []
    lines.append(f"## {subclade['name']}")
    lines.append(f" * parent: [{subclade['parent']}](#{subclade['parent'].replace('.', '')})")
    snp_str = ', '.join(f"{x['locus']}:{x['position']}{x['state']}" for x in subclade['defining_mutations'])
    lines.append(f" * defining mutations or substitutions: {snp_str}")
    if "clade" in subclade and subclade['clade'] != "none":
        lines.append(f" * clade: {subclade['clade']}")

    ref_seqs = []
    for x in subclade['representatives']:
        nextstrain_link = f"[View on Nextstrain](https://nextstrain.org/flu/seasonal/{lineage}/{segment}/6y?c=subclade&s={x['isolate']})"
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
        lines.append(f" * representative sequences:")
        for r in ref_seqs:
            lines.append(f"   - {r}")
    lines.append(f" * [View on Nextstrain](https://nextstrain.org/flu/seasonal/{lineage}/{segment}/6y?branchLabel=Subclade&c=subclade&label=Subclade:{subclade['name']})")
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
                         for x in subclades if 'clade' in x and x['clade'] != 'none']
    # Write to json file
    with open('.auto-generated/subclades.md', 'w') as outfile:
        outfile.write("# Summary of designated subclades\n")

        for subclade in subclades:
            print("output clade", subclade['name'])
            outfile.write(generate_lineage_md(subclade, args.lineage, args.segment) + '\n')

        if len(clade_lineage_map):
            # write table of clade -- subclade correspondence
            outfile.write("# Clade -- subclade correspondence\n")
            outfile.write(f"|*Clade*|*Subclade*|*full subclade name*|\n")
            outfile.write(f"|-------------|---------|----------------------|\n")
            for clade, lineage, unaliased_name in clade_lineage_map:
                outfile.write(f"|{clade}|[{lineage}](#{lineage.replace('.','')})|{unaliased_name}|\n")

