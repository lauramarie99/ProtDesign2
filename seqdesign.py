import utils, glob, yaml, argparse, json, os, config

# Creates a json file given a dictionary
def create_json(outpath, dict):
    with open(outpath, 'w') as file:
        json.dump(dict, file)

# Creates the json input files for sequence design
def create_input_files(path, name, contig_str, outdir):
    pdb_files = glob.glob(f"{path}/{name}*.pdb")
    pdb_dict = {path: "" for path in pdb_files}
    fixed_resi = []
    index = 0
    for block in contig_str.split("/"):
        if block[0].isalpha():
            lb = int(block.split("-")[0][1:])
            ub = int(block.split("-")[1])
            range = ub-lb+1
            i = 0
            while i < range:
                index += 1
                fixed_resi.append(block[0] + str(index))
                i += 1
        else:
            index = index + int(block.split("-")[0])
        fixed_resi_str = " ".join(fixed_resi)
        resi_dict = {path: fixed_resi_str for path in pdb_files}
    
    os.makedirs(outdir, exist_ok=True)
    create_json(f"{outdir}/pdb_ids.json", pdb_dict)
    create_json(f"{outdir}/fix_residues_multi.json", resi_dict)

# Postprocessing of fasta files
def postprocessing(outdir, name):
    fasta_seq = []
    for file in glob.glob(f"{outdir}/*.fa"):
        with open(file, "r") as infile:
            lines = infile.readlines()
        save_lines = lines[2:]
        for i in range(0,len(save_lines)-1,2):
            save_lines[i] = save_lines[i].split(",")[0] + "_n" + save_lines[i].split(",")[1][4:] + "\n"
        save_lines[-1] = save_lines[-1] + "\n"
        fasta_seq.extend(save_lines)

    # Write the remaining entries to the output file
    with open(f"{outdir}/{name}.fa", "w") as outfile:
        outfile.writelines(fasta_seq)

"""
MAIN
"""

# Read given config
parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True)
args = parser.parse_args()
args = yaml.safe_load(open(args.config))
args_seqdesign = args["seqdesign"]
contigs_str = args["diffusion"]["contigs"]
num_seqs = args_seqdesign["num_seqs"]
model_type = args_seqdesign["model_type"]
path = args["diffusion"]["path"]
name = args["diffusion"]["name"]
full_path = f"{path}/{name}"
outdir = f"{full_path}/SeqDesign"

# Create input files
create_input_files(f"{full_path}/Diffusion", name, contigs_str, f"{outdir}/inputs")

# Run sequence design
opts = [f"--model_type {model_type}",
        f"--out_folder {outdir}/outputs",
        f"--batch_size {num_seqs}",
        f"--pdb_path_multi {outdir}/inputs/pdb_ids.json",
        f"--fixed_residues_multi {outdir}/inputs/fix_residues_multi.json",
        f"--zero_indexed 1"]
if "seed" in args_seqdesign: opts.append(f"--seed {args_seqdesign['seed']}")
if "temperature" in args_seqdesign: opts.append(f"--temperature {args_seqdesign['temperature']}")
opts = ' '.join(opts)

print("running sequence design...")
cmd = f"cd {config.LIGANDMPNN_PATH} && python3.9 run.py {opts}"
print(cmd)
utils.run(cmd)

# Postprocessing of resulting fasta files
postprocessing(f"{outdir}/outputs/seqs", name)