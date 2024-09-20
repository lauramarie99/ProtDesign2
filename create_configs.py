import yaml, argparse, random, utils, copy


def create_random_contig(contig):
    new_contig_sections = []
    sections = contig.split("/")
    for section in sections:
        if section[0].isalpha():
            new_contig_sections.append(section)
        else:
            lb = int(section.split('-')[0])
            ub = int(section.split('-')[1])
        if lb != ub:
            random_number = random.randint(lb,ub)
            new_contig_sections.append(str(random_number) + "-" + str(random_number))
        else:
            new_contig_sections.append(str(lb) + "-" + str(lb))
    new_contig = '/'.join(new_contig_sections)
    return new_contig


# Get config
parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True)
parser.add_argument('--num_contigs', type=int, required=True)
parser.add_argument('--outdir', type=str, required=True)
args = parser.parse_args()
with open(args.config, "r") as file:
    config_args = yaml.safe_load(file)
name = args.config.split("/")[-1][:-4]

general_contig = config_args["diffusion"]["contigs"]

for n in range(args.num_contigs):
    new_contig = create_random_contig(general_contig)
    new_config_args = copy.deepcopy(config_args)
    new_config_args["diffusion"]["contigs"] = new_contig
    new_config_args["diffusion"]["name"] = f"{config_args['diffusion']['name']}{n+1}"

    if 'cst_file' in config_args['seqdesign']:
        model_motif, ref_motif = utils.get_motifs(new_contig)
        model_motif = [(resi[1:]+resi[0]) for resi in model_motif]
        print(model_motif)
        ref_motif = [(resi[1:]+resi[0]) for resi in ref_motif]
        print(ref_motif)
        with open(config_args['seqdesign']['cst_file']) as infile:
            data = infile.readlines()
        new_data = []
        for line in data:
            new_line = line
            for i, resi in enumerate(ref_motif):
                if resi in line:
                    new_line = new_line.replace(resi, model_motif[i])
            new_data.append(new_line)
        with open(f"{args.outdir}/constraints_{n+1}.cst", "w") as outfile:
            outfile.writelines(new_data)

        new_config_args["seqdesign"]["cst_file"] = f"{args.outdir}/constraints_{n+1}.cst"

    with open(f"{args.outdir}/{name}_{n+1}.yml", "w") as outfile:
        yaml.dump(new_config_args, outfile)


