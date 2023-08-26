import yaml, glob

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir')
    parser.add_argument('--flat-output', action='store_true', default=False)
    parser.add_argument('--output-tsv')
    args = parser.parse_args()
    yml_files = glob.glob(args.input_dir+'/*yml')

    clades = {}
    for yfile in yml_files:
        with open(yfile, 'r') as stream:
            yaml_data = yaml.safe_load(stream)
            clades[yaml_data['name']] = {'parent': yaml_data['parent'],
                                        'defining_muts':{(x['locus'], x['position']):x['state']
                                        for x in yaml_data['defining_mutations']}}


    tsv_file = open(args.output_tsv, 'w')
    sep = '\t'
    tsv_file.write(sep.join(['clade','gene','site','alt'])+'\n')

    if args.flat_output:
        all_muts = {}
        for c in sorted(clades.keys()):
            if clades[c].get('parent', 'none')=='none':
                all_muts[c] = clades[c]['defining_muts']
            else:
                all_muts[c] = {k:v for k,v in all_muts[clades[c]['parent']].items()}
            for (locus, position), state in clades[c]['defining_muts'].items():
                all_muts[c][(locus, position)] = state

            for (locus, position), state in all_muts[c].items():
                tsv_file.write(sep.join([c, locus, str(position),state])+'\n')

    else:
        for c in sorted(clades.keys()):
            if clades[c].get('parent', 'none')!='none':
                tsv_file.write(sep.join([c,'clade', clades[c]['parent'],''])+'\n')
            for (locus, position), state in clades[c]['defining_muts'].items():
                tsv_file.write(sep.join([c, locus, str(position),state])+'\n')

    tsv_file.close()