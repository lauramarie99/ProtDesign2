import random, string, os, yaml, argparse, utils, config, glob

# Add ligand and sidechains of motif to designs
def postprocessing(contigs, name, path, pdb, substrate):
    full_path = f"{path}/{name}/Diffusion"
    design_motif, ref_motif, redesigned_residues = utils.get_motifs(contigs)
    ref_path = pdb
    for design_path in glob.glob(f"{full_path}/{name}*.pdb"):
        if "/0" in contigs:
            # Remove second chain
            utils.remove_chain_from_pdb(design_path=design_path,
                                        chain_to_remove="B")
        utils.add_sidechain_and_ligand_coordinates(design_path=design_path,
                                                   ref_path=ref_path,
                                                   design_motif=design_motif,
                                                   ref_motif=ref_motif,
                                                   ligand_name=substrate)

# Run RFdiffusion
def run_diffusion(type, 
                  contigs, 
                  name, 
                  path,
                  pdb=None, 
                  iterations=50,
                  num_designs=10,
                  enzyme_design=False,
                  guide_scale=1,
                  guide_potentials="",
                  substrate="",
                  ckpt_override_path=None,
                  noise_scale=1,
                  deterministic=False,
                  partial_diffusion=False):
    """
    This function runs a diffusion simulation using provided input parameters and generates the final PDB structures.

    Args:
    contigs (str): Input contigs string to define the fixed and free portions.
    name (str): Experiment name.
    path (str): The output directory path for generated results.
    pdb (str, optional): The PDB file path. Defaults to None.
    iterations (int, optional): Number of diffusion iterations. Defaults to 50.
    num_designs (int, optional): Number of designs to generate. Defaults to 10.
    enzyme_design (bool, optional): If True, generates substrate pockets by adding guiding potential. Defaults to False.
    guide_scale (float): Scaling factor for guiding potentials. Defaults to 1.
    guide_potentials (str): The guiding potentials string. Defaults to an empty string.
    substrate (str): Substrate for the guiding potential. Defaults to "".
    ckpt_override_path (str, optional): The path of the checkpoint file. Defaults to None.
    noise_scale (int, optional): Change noise_scale_ca and noise_scale_frame. Defaults to 1.
    deterministic (bool, optional): Deterministic initialization. Defaults to False.
    partial_diffusion (bool, optional): Carry out partial_diffusion. Defaults to False.
    
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
            f"'contigmap.contigs=[{contigs}]'"]

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

    # Deterministic run
    if deterministic:
        opts.append(f"inference.deterministic=True")

    # Store input PDB in outdir
    if pdb:
        pdb_filename = f"{full_path}/input.pdb"
        os.system(f"cp {pdb} {pdb_filename}")
        opts.append(f"inference.input_pdb={pdb_filename}")
    if ckpt_override_path:
        opts.append(f"inference.ckpt_override_path={ckpt_override_path}")

    # Print different parameters for diagnostic purposes
    print("output:", full_path)
    print("contigs:", contigs)

    # Run the inference script
    opts_str = " ".join(opts)
    cmd = f"cd {config.RFDIFFUSION_PATH} && python3.9 run_inference.py {opts_str}"
    print(cmd)
    utils.run(cmd)
    if pdb and substrate:
        postprocessing(contigs, name, path, pdb, substrate)

# Run RFdiffusion all-atom
def run_diffusion_aa(type, 
                     contigs, 
                     name, 
                     path,
                     pdb,
                     substrate, 
                     iterations=50,
                     num_designs=10,
                     noise_scale=1,
                     deterministic=False):
    """
    This function runs a diffusion-all-atom simulation using provided input parameters and generates the final PDB structures.

    Args:
    contigs (str): Input contigs string to define the fixed and free portions.
    name (str): Experiment name.
    path (str): The output directory path for generated results.
    pdb (str): The PDB file path.
    substrate (str): The substrate to build a pocket around.
    iterations (int, optional): Number of diffusion iterations. Defaults to 50.
    num_designs (int, optional): Number of designs to generate. Defaults to 10.
    noise_scale (int, optional): Change noise_scale_ca and noise_scale_frame. Defaults to 1.
    deterministic (bool, optional): Deterministic initialization. Defaults to False.
    
    """

    # Make output directory
    full_path = f"{path}/{name}/Diffusion"
    os.makedirs(full_path, exist_ok=True)
    output_prefix = f"{full_path}/{name}"

    # Preprocess contig string
    contigs = contigs.replace("/", ",")

    # Add general options
    opts = [f"inference.output_prefix={output_prefix}", 
            f"inference.num_designs={num_designs}",
            f"denoiser.noise_scale_ca={noise_scale}",
            f"denoiser.noise_scale_frame={noise_scale}",
            f"diffuser.T={iterations}",
            f"inference.ligand={substrate}",
            f"contigmap.contigs=[\\'{contigs}\\']"]

    # Deterministic run
    if deterministic:
        opts.append(f"inference.deterministic=True")
    
    # Store input PDB in outdir
    pdb_filename = f"{full_path}/input.pdb"
    os.system(f"cp {pdb} {pdb_filename}")
    opts.append(f"inference.input_pdb={pdb_filename}")
    
    # Print different parameters for diagnostic purposes
    print("output:", full_path)
    print("contigs:", contigs)

    # Run the inference script
    opts_str = " ".join(opts)
    cmd = f"cd {config.RFDIFFUSIONAA_PATH} && python run_inference.py {opts_str}"
    print(cmd)
    utils.run(cmd)


"""
MAIN
"""

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

if not os.path.exists(f"{path}/{name}/Diffusion/{name}_0.pdb"):

    # Get diffusion arguments
    for k,v in args_diffusion.items():
        if isinstance(v,str):
            args_diffusion[k] = v.replace("'","").replace('"','')

    # Run diffusion
    if args_diffusion["type"] == "all-atom":
        run_diffusion_aa(**args_diffusion)
    else:
        run_diffusion(**args_diffusion)

    # Copy config file to results directory
    os.system(f"cp {config_file} {path}/{name}")