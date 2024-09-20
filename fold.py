import utils, argparse, yaml, os, shutil, glob, json
from statistics import mean

def get_mean_plddt(dir, motif):
    json_file = glob.glob(f"{dir}/*rank_001*.json")[0]
    with open(json_file, "r") as f:
        data = json.load(f) 
    mean_plddt = round(mean(data["plddt"]), 2)
    motif_plddt = [data["plddt"][i] for i in motif]
    mean_motif_plddt = round(mean(motif_plddt),2)
    return mean_plddt, mean_motif_plddt

def get_scores(dir, design_id, ref_path, diff_path, model_motif, ref_motif):
    score_dict = {}
    model_motif_index_list = [int(resi[1:]) for resi in model_motif]
    model_path = glob.glob(f"{dir}/*rank_001*.pdb")[0]
    score_dict["id"] = design_id
    score_dict["mean-plddt"],score_dict["mean-motif-plddt"] = get_mean_plddt(dir, model_motif_index_list)
    score_dict["ca-rmsd"] = utils.get_ca_rmsd(model_path=model_path, ref_path=diff_path)
    score_dict["motif-ca-rmsd"] = utils.get_motif_ca_rmsd(model_path=model_path, ref_path=ref_path, model_motif=model_motif, ref_motif=ref_motif)
    score_dict["motif-all-atom-rmsd"] = utils.get_motif_all_atom_rmsd(model_path=model_path, ref_path=ref_path, model_motif=model_motif, ref_motif=ref_motif)
    return score_dict


def create_score_file(outdir, args_diffusion):
    model_motif_list, ref_motif_list = utils.get_motifs(args_diffusion["contigs"])
    subfolders = [f for f in glob.glob(f"{outdir}/*") if os.path.isdir(f)]
    print(subfolders)
    combined_scores = {}
    for folder in subfolders:
        design_id = folder.split("/")[-1]
        diff_id = design_id.split("_")[0] + "_" + design_id.split("_")[1] + ".pdb"
        diff_path = f"{outdir}/../Diffusion/{diff_id}"
        print(design_id)
        scores = get_scores(dir=folder, 
                            design_id=design_id, 
                            ref_path=args_diffusion["pdb"], 
                            diff_path=diff_path, 
                            model_motif=model_motif_list, 
                            ref_motif=ref_motif_list)
        combined_scores[design_id] = scores
    print(combined_scores)
    with open(f"{outdir}/scores.json", "w") as output_file:
        json.dump(combined_scores, output_file, indent=4)

# Postprocessing: Save resulting files for each diffused model in one separate directory
def postprocessing(outdir, args_diffusion):
    for filename in os.listdir(outdir):
        if filename.endswith(".done.txt") or filename.endswith(".a3m"):
            file_path = os.path.join(outdir, filename)
            os.remove(file_path)
            continue
        elif filename.endswith("bibtex") or filename.endswith("config.json") or filename.endswith("log.txt"):
            continue
        else:
            id = filename.split("_")[0] + "_" + filename.split("_")[1] + "_" + filename.split("_")[2] + "_" + filename.split("_")[3]
            os.makedirs(os.path.join(outdir, id), exist_ok=True)
            shutil.move(os.path.join(outdir, filename), os.path.join(outdir, id))
    create_score_file(outdir=outdir, args_diffusion=args_diffusion)

"""
MAIN
"""

# Read given config
parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True)
args = parser.parse_args()
args = yaml.safe_load(open(args.config))
args_folding = args["folding"]
msa_mode = args_folding["msa_mode"]
num_models = args_folding["num_models"]
recycles = args_folding["num_recycles"]
path = args["diffusion"]["path"]
name = args["diffusion"]["name"]
full_path = f"{path}/{name}"
outdir = f"{full_path}/Folding"
input = f"{full_path}/SeqDesign/{name}.fa"

# Run ColabFold
opts = [f"--msa-mode {msa_mode}",
        f"--num-models {num_models}",
        f"--num-recycle {recycles}",
        f"{input} {outdir}"]
opts = ' '.join(opts)

# print("running ColabFold...")
# cmd = f"colabfold_batch {opts}"
# print(cmd)
# utils.run(cmd)

# Postprocessing of resulting files
postprocessing(outdir=outdir, args_diffusion=args["diffusion"])
