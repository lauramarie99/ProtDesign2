import random, string, os, yaml, argparse, utils, config

# Run diffusion
def run_diffusion(type, contigs, name, path,
                  pdb=None, 
                  iterations=50,
                  num_designs=10,
                  guide_scale=1,
                  guide_potentials="",
                  substrate="2PE",
                  ckpt_override_path="null",
                  enzyme_design=False,
                  partial_diffusion=False,
                  noise_scale=1,
                  deterministic=False):
    """
    This function runs a diffusion simulation using provided input parameters, 
    applies contigs processing, and generates the final PDB structures.

    Args:
    contigs (str): Input contigs string to define the fixed and free portions.
    name (str): Experiment name.
    path (str): The output directory path for generated results.
    pdb (str, optional): The PDB file path. Defaults to None.
    iterations (int, optional): Number of diffusion iterations. Defaults to 50.
    num_designs (int, optional): Number of designs to generate. Defaults to 10.
    guide_scale (float): Scaling factor for guiding potentials. Defaults to 1.
    guide_potentials (str): The guiding potentials string. Defaults to an empty string.
    substrate (str): The substrate design. Defaults to "2PE".
    ckpt_override_path (str): The path of the checkpoint file. Defaults to "null".
    enzyme_design (bool, optional): If True, generates substrate pockets by adding guiding potential. Defaults to False.
    noise_scale (int, optional): Change noise_scale_ca and noise_scale_frame.
    deterministic (bool, optional): Deterministic initialization.
    partial_diffusion (bool, optional): Carry out partial_diffusion
    
    Returns:
    tuple: The updated contigs list and the number of symmetry-equivalent copies.
    """

    # Make output directory
    full_path = f"{path}/{name}/Diffusion"
    os.makedirs(full_path, exist_ok=True)
    output_prefix = f"{full_path}/{name}"

    # Add general options
    opts = [f"inference.output_prefix={output_prefix}", 
            f"inference.num_designs={num_designs}",
            f"denoiser.noise_scale_ca={noise_scale}",
            f"denoiser.noise_scale_frame={noise_scale}",
            f"inference.ckpt_override_path={ckpt_override_path}"]

    copies = 1
    
    # Store input PDB in outdir
    if pdb:
        pdb_filename = f"{full_path}/input.pdb"
        os.system(f"cp {pdb} {pdb_filename}")
        opts.append(f"inference.input_pdb={pdb_filename}")

    # Add contig to options
    opts.append(f"'contigmap.contigs=[{contigs}]'")

    # Add enzyme_design related options if enzyme_design is True
    if enzyme_design:
        opts.append(f"potentials.guide_scale={guide_scale}")
        opts.append(f"'potentials.guiding_potentials=[\"{guide_potentials}\"]'")
        opts.append(f"potentials.substrate={substrate}")

    # Add number of diffusion steps
    if partial_diffusion:
        opts.append(f"diffuser.partial_T={iterations}")
    else:
        opts.append(f"diffuser.T={iterations}")

    if deterministic:
        opts.append(f"inference.deterministic=True")

    # Print different parameters for diagnostic purposes
    print("output:", full_path)
    print("contigs:", contigs)

    # Create the command with options to run the inference script
    opts_str = " ".join(opts)
    cmd = f"cd {config.RFDIFFUSION_PATH} && python3.9 run_inference.py {opts_str}"
    print(cmd)
    # Run the command using a helper function "run"
    utils.run(cmd)
    return contigs, copies

# Run diffusion
def run_diffusion_aa(type, contigs, name, path,
                    pdb=None, 
                    iterations=50,
                    num_designs=10,
                    noise_scale=1,
                    deterministic=False,
                    substrate="2PE"):
    """
    This function runs a diffusion-all-atom simulation using provided input parameters, 
    applies contigs processing, and generates the final PDB structures.

    Args:
    contigs (str): Input contigs string to define the fixed and free portions.
    name (str): Experiment name.
    path (str): The output directory path for generated results.
    pdb (str, optional): The PDB file path. Defaults to None.
    iterations (int, optional): Number of diffusion iterations. Defaults to 50.
    num_designs (int, optional): Number of designs to generate. Defaults to 10.
    substrate (str): The substrate to build a pocket around. Defaults to "2PE".
    noise_scale (int, optional): Change noise_scale_ca and noise_scale_frame.
    deterministic (bool, optional): Deterministic initialization.
    
    Returns:
    tuple: The updated contigs list and the number of symmetry-equivalent copies.
    """

    # Make output directory
    full_path = f"{path}/{name}/Diffusion"
    os.makedirs(full_path, exist_ok=True)
    output_prefix = f"{full_path}/{name}"
    contigs = contigs.replace("/", ",")
    copies = 1

    # Add general options
    opts = [f"inference.output_prefix={output_prefix}", 
            f"inference.num_designs={num_designs}",
            f"denoiser.noise_scale_ca={noise_scale}",
            f"denoiser.noise_scale_frame={noise_scale}",
            f"diffuser.T={iterations}",
            f"inference.ligand={substrate}"]

    opts.append(f"contigmap.contigs=[\\'{contigs}\\']")

    pdb_filename = f"{full_path}/input.pdb"
    os.system(f"cp {pdb} {pdb_filename}")
    opts.append(f"inference.input_pdb={pdb_filename}")

    if deterministic:
        opts.append(f"inference.deterministic=True")

    print("output:", full_path)
    print("contigs:", contigs)

    # Create the command with options to run the inference script
    opts_str = " ".join(opts)
    cmd = f"cd {config.RFDIFFUSIONAA_PATH} && python run_inference.py {opts_str}"
    print(cmd)
    # Run the command using a helper function "run"
    utils.run(cmd)

    return contigs, copies


# Read given config
parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True)
args = parser.parse_args()
config_file = args.config
args = yaml.safe_load(open(config_file))
args_diffusion = args["diffusion"]

# Check if output directory already exists
name = args_diffusion["name"]
path = args_diffusion["path"]
if os.path.exists(f"{path}/{name}/Diffusion/{name}_0.pdb"):
  args_diffusion["name"] = name = args_diffusion["name"] + "_" + ''.join(random.choices(string.ascii_lowercase + string.digits, k=5))

# Get diffusion arguments
for k,v in args_diffusion.items():
  if isinstance(v,str):
    args_diffusion[k] = v.replace("'","").replace('"','')

# Run diffusion
if args_diffusion["type"] == "all-atom":
     contigs, copies = run_diffusion_aa(**args_diffusion)
else:
    contigs, copies = run_diffusion(**args_diffusion)

# Copy config to results directory
os.system(f"cp {config_file} {path}/{name}")

# Print output contigs
print("the final contigs are:")
print(contigs, copies)