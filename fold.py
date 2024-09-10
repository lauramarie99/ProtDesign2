import utils, argparse, yaml, os, shutil

# Postprocessing: Save resulting files for each diffused model in one separate directory
def postprocessing(outdir):
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
input = f"{full_path}/SeqDesign/outputs/seqs/{name}.fa"

# Run ColabFold
opts = [f"--msa-mode {msa_mode}",
        f"--num-models {num_models}",
        f"--num-recycle {recycles}",
        f"{input} {outdir}"]
opts = ' '.join(opts)

print("running ColabFold...")
cmd = f"colabfold_batch {opts}"
print(cmd)
utils.run(cmd)

# Postprocessing of resulting files
postprocessing(outdir)
