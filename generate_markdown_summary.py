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
    if "clade" in lineage and subclade['clade'] != "none":
        lines.append(f" * clade: {subclade['clade']}")

    ref_seqs = [f"[{x.get('isolate', x['accession'])}](https://www.ncbi.nlm.nih.gov/nuccore/{x['accession']})" for x in subclade['representatives'] if x['source']=='genbank']
    if len(ref_seqs)==1:
        lines.append(f" * representative sequence: {ref_seqs[0]}")
    elif len(ref_seqs)>1:
        lines.append(f" * representative sequences:")
        for r in ref_seqs:
            lines.append(f"   - {r}")
    lines.append(f" * [View in Nextstrain](https://nextstrain.org/flu/seasonal/{lineage}/{segment}/6y?branchLabel=Subclade&c=subclade&label=Subclade:{subclade['name']})")
    return '\n'.join(lines) + '\n'

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--lineage')
    parser.add_argument('--segment')
    args = parser.parse_args()
    repo_name = f"{args.lineage}_{args.segment}"

    subclades = []
    # Iterate through all lineage definition files
    for yaml_file in sorted(glob.glob(f"subclades/*.yml")):
        with open(yaml_file, 'r') as stream:
            yaml_data = yaml.safe_load(stream)
        subclades.append(yaml_data)

    subclades.sort(key=lambda x:x['name'])
    # Write to json file
    with open('.auto-generated/subclades.md', 'w') as outfile:
        outfile.write("# Summary of designated subclades\n")

        for subclade in subclades:
            outfile.write(generate_lineage_md(subclade, args.lineage, args.segment) + '\n')

