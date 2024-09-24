import utils, glob, yaml, argparse, json, os, config, time
import pyrosetta as pyr
from Bio.PDB import PDBParser, Superimposer
from concurrent.futures import ProcessPoolExecutor, as_completed
start = time.time()

# Creates a json file given a dictionary
def create_json(outpath, dict):
    with open(outpath, 'w') as file:
        json.dump(dict, file)


# Filter PDB files and create the json input files for sequence design
def preprocessing(inpath, name, outpath, contig_str, ref_path, threshold):
    pdb_files = glob.glob(f"{inpath}/{name}*.pdb")                                  # Get all pdb file paths
    model_motif, ref_motif, redesigned_residues = utils.get_motifs(contig_str)      # Get motifs
    filtered_pdb_files = []                                                         # Filter pdb files based on Motif Ca-RMSD
    for path in pdb_files:
        rmsd = utils.get_motif_ca_rmsd(model_path=path,ref_path=ref_path,model_motif=model_motif,ref_motif=ref_motif)
        if rmsd <= threshold:
            filtered_pdb_files.append(path)
    print("Filtered PDB files: ", filtered_pdb_files)
    pdb_dict = {path: "" for path in filtered_pdb_files}                            # Create dictionaries with paths and motif
    fixed_resi_str = " ".join(model_motif)
    fixed_resi_dict = {path: fixed_resi_str for path in filtered_pdb_files}
    redesigned_resi_str = " ".join(redesigned_residues)
    redesigned_resi_dict = {path: redesigned_resi_str for path in filtered_pdb_files}
    os.makedirs(outpath, exist_ok=True)
    create_json(f"{outpath}/pdb_ids.json", pdb_dict)                                # Create json input files
    create_json(f"{outpath}/fix_residues_multi.json", fixed_resi_dict)
    create_json(f"{outpath}/redesigned_residues_multi.json", redesigned_resi_dict)
    return model_motif, ref_motif


# Postprocessing of fasta files
def postprocessing(outpath, name, relax_round):
    fasta_seq = []
    fasta_files = glob.glob(f"{outpath}/*.fa")
    print(fasta_files)
    for file in fasta_files:
        with open(file, "r") as infile:
            lines = infile.readlines()
        save_lines = lines[2:]
        for i in range(0,len(save_lines)-1,2):
            id = save_lines[i].split(",")[0]
            if relax_round > 0:
                save_lines[i] = id + "\n"
            else:
                save_lines[i] = id + "_n" + save_lines[i].split(",")[1][4:] + "_c" + str(relax_round) + "\n"
        save_lines[-1] = save_lines[-1] + "\n"
        fasta_seq.extend(save_lines)

    # Write the remaining entries to the output file
    with open(f"{outpath}/{name}_c{str(relax_round)}.fa", "w") as outfile:
        outfile.writelines(fasta_seq)

# Sequence design
def design(inpath, outpath, args_seqdesign, args_diffusion, relax_round):
    # Preprocessing
    model_motif, ref_motif = preprocessing(inpath=inpath, 
                                           name=args_diffusion["name"], 
                                           contig_str=args_diffusion["contigs"], 
                                           outpath=f"{outpath}/inputs", 
                                           ref_path=args_diffusion["pdb"], 
                                           threshold=int(args_seqdesign["rmsd_cutoff"]))

    # Run sequence design
    opts = [f"--model_type {args_seqdesign['model_type']}",
            f"--out_folder {outpath}/outputs",
            f"--number_of_batches {args_seqdesign['num_seqs']}",
            f"--pdb_path_multi {outpath}/inputs/pdb_ids.json",
            f"--fixed_residues_multi {outpath}/inputs/fix_residues_multi.json",
            f"--redesigned_residues_multi {outpath}/inputs/redesigned_residues_multi.json",
            f"--zero_indexed 1",
            "--pack_side_chains 1",
            "--number_of_packs_per_design 1"]
    if args_diffusion["type"] != "all-atom": opts.append("--repack_everything 1")
    if "seed" in args_seqdesign: opts.append(f"--seed {args_seqdesign['seed']}")
    if "temperature" in args_seqdesign: opts.append(f"--temperature {args_seqdesign['temperature']}")
    opts = ' '.join(opts)
    print("running sequence design...")
    cmd = f"cd {config.LIGANDMPNN_PATH} && python3.9 run.py {opts}"
    print(cmd)
    utils.run(cmd)

    # Postprocessing
    postprocessing(outpath=f"{outpath}/outputs/seqs", name=args_diffusion["name"], relax_round=relax_round)
    return outpath, model_motif, ref_motif

# Create specific cst file for input pdb
def create_cst_file(model_motif, ref_motif, old_cst_file, pdb_file, ligand, outdir, name):
    model_motif = [(resi[1:]+resi[0]) for resi in model_motif]
    ref_motif = [(resi[1:]+resi[0]) for resi in ref_motif]
    ligand_index = utils.get_ligand_index(pdb_file=pdb_file, ligand=ligand)
    with open(old_cst_file) as infile:
        data = infile.readlines()
    new_data = []
    for line in data:
        new_line = line
        blocks = line.split(" ")
        block1 = blocks[2].strip()
        block2 = blocks[4].strip()
        for resi in [block1, block2]:
            if resi in ref_motif:
                index = ref_motif.index(resi)
                new_line = new_line.replace(resi, model_motif[index])
            else:
                new_line = new_line.replace(resi, f"{ligand_index}B")
        new_data.append(new_line)
    filename = f"{outdir}/{name}.cst"
    with open(filename, "w") as outfile:
        outfile.writelines(new_data)
    return filename


# Rosetta FastRelax
def relax(pdb, id, outpath, model_motif, ref_motif, ligand, params=None, cst=None):

    # Create cst file
    if cst is not None:
        cst_file = create_cst_file(model_motif=model_motif, 
                                   ref_motif=ref_motif, 
                                   old_cst_file=cst, 
                                   pdb_file=pdb,
                                   ligand=ligand,
                                   outdir=f"{outpath}/relaxed", 
                                   name="constraints")

    # Initialization
    extra_res_fa = f"-extra_res_fa {params}" if params is not None else ""
    constraints = f"-constraints:cst_fa_file {cst}" if cst is not None else ""
    pyr.init(f"-ignore_zero_occupancy false -ex1 -ex2 {extra_res_fa} {constraints}")

    # Define scorefunction
    scorefxn = pyr.get_fa_scorefxn()                                                                    # Default full-atom energy terms
    scorefxn.set_weight(pyr.rosetta.core.scoring.score_type_from_name("atom_pair_constraint"), 1)       # Add constraint weight to score function    
    scorefxn.set_weight(pyr.rosetta.core.scoring.score_type_from_name("coordinate_constraint"), 1)
    # Read input pdb file
    pose = pyr.pose_from_file(pdb)
    pose2 = pose.clone()

    # Add constraints
    if cst is not None:
        constraint_mover = pyr.rosetta.protocols.constraint_movers.ConstraintSetMover()
        constraint_mover.constraint_file(cst_file)
        constraint_mover.apply(pose2)


    # Specify task operations
    tf = pyr.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyr.rosetta.core.pack.task.operation.InitializeFromCommandline())     # Use options specified in the init() function
    tf.push_back(pyr.rosetta.core.pack.task.operation.RestrictToRepacking())           # Avoid design of residues
    tf.push_back(pyr.rosetta.core.pack.task.operation.IncludeCurrent())                # Includes the current rotamers to the rotamers sets   
    tf.push_back(pyr.rosetta.core.pack.task.operation.NoRepackDisulfides())            # Prevent repacking of disulfid bond side-chains   
    residue_selector = pyr.rosetta.core.select.residue_selector.ResidueIndexSelector() # Prevent repacking of catalytic residues 
    for resi in model_motif:
        residue = int(resi[1:].strip())
        residue_selector.append_index(residue)
    residue_selector.apply(pose2)
    chain_selector = pyr.rosetta.core.select.residue_selector.ChainSelector("B")
    chain_selector.apply(pose2)
    no_repacking_selector = pyr.rosetta.core.select.residue_selector.OrResidueSelector()
    no_repacking_selector.add_residue_selector(residue_selector)
    no_repacking_selector.add_residue_selector(chain_selector)
    no_repacking_selector.apply(pose2)
    prevent_repacking_rlt = pyr.rosetta.core.pack.task.operation.PreventRepackingRLT()
    prevent_repacking_op = pyr.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, no_repacking_selector, False)
    tf.push_back(prevent_repacking_op)
    print(tf.create_task_and_apply_taskoperations(pose2))

    # Specify which torsions to minimize
    mmf = pyr.rosetta.core.select.movemap.MoveMapFactory()
    mmf.add_chi_action(pyr.rosetta.core.select.movemap.mm_disable, no_repacking_selector)
    mmf.add_bb_action(pyr.rosetta.core.select.movemap.mm_disable, no_repacking_selector)

    #mm = pyr.rosetta.core.kinematics.MoveMap()
    #mm.set_jump(True)
    #mm.set_chi(True)
    #mm.set_bb(True)

    # Perform FastRelax
    fastRelax = pyr.rosetta.protocols.relax.FastRelax()
    fastRelax.constrain_relax_to_start_coords(True)
    fastRelax.set_scorefxn(scorefxn)
    fastRelax.set_movemap_factory(mmf)
    #fastRelax.set_movemap(mm)
    fastRelax.set_task_factory(tf)
    before = scorefxn.score(pose2)
    print(scorefxn.show(pose2))
    fastRelax.apply(pose2)
    after = scorefxn.score(pose2)
    pose2.dump_pdb(f"{outpath}/relaxed/{id}.pdb")
    print("Before: ", before)
    print("After: ", after)
    print("Delta: ", after-before)
    print(scorefxn.show(pose2))
    return outpath
        
# Call relax function
def relax_pdb(pdb, relax_round, outpath, args_seqdesign, model_motif, ref_motif, ligand):
    id = pdb.split("/")[-1].split("_packed")[0]
    if "_c" in id:
        index = id.index("_c")
        design_id = id[:index+2] + str(relax_round) + id[index+3:]
    else:
        design_id = id + "_n" + pdb.split("/")[-1].split("packed_")[1][:-6] + "_c" + str(relax_round)
    outpath = relax(pdb=pdb, 
                    id=design_id, 
                    outpath=outpath, 
                    **({"params": args_seqdesign["params_file"]} if "params_file" in args_seqdesign else {}),
                    **({"cst": args_seqdesign["cst_file"]} if "cst_file" in args_seqdesign else {}),
                    ligand=ligand,
                    model_motif=model_motif,
                    ref_motif=ref_motif)


# Perform cycles of Rosetta FastRelax and Sequence Design
def relax_and_design(inpath, outpath, args_seqdesign, args_diffusion, relax_round, model_motif, ref_motif):
    pdb_files = glob.glob(f"{inpath}/*.pdb")                                                            # Extract all pdb files
    os.makedirs(os.path.join(outpath, "relaxed"), exist_ok=True)                                        # Create outdir
    
    # Relaxation (in batches)
    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(relax_pdb, pdb, relax_round, outpath, args_seqdesign, model_motif, ref_motif, args_diffusion["substrate"]) for pdb in pdb_files]
        for future in as_completed(futures):
            future.result()

    # Proceed with the design function after all relax calls are complete
    outpath, model_motif, ref_motif = design(inpath=f"{outpath}/relaxed", 
                                             outpath=outpath, 
                                             args_seqdesign=args_seqdesign, 
                                             args_diffusion=args_diffusion,
                                             relax_round=relax_round)
    return outpath

def concat_fasta_files(outpath, name):
    all_lines = []
    fasta_files = glob.glob(f"{outpath}/Recycle-*/outputs/seqs/{name}_c*.fa")
    print(fasta_files)
    for fasta in fasta_files:
        with open(fasta, "r") as infile:
            lines = infile.readlines()
            all_lines.extend(lines)
    
    # Write all entries to the output file
    with open(f"{outpath}/{name}.fa", "w") as outfile:
        outfile.writelines(all_lines)


"""
MAIN
"""

def main():
    # Read given config
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True)
    args = parser.parse_args()
    args = yaml.safe_load(open(args.config))
    args_seqdesign = args["seqdesign"]
    args_diffusion = args["diffusion"]

    # Define outpath
    full_path = f"{args_diffusion['path']}/{args_diffusion['name']}"

    # Perform sequence design
    outpath, model_motif, ref_motif = design(inpath=f"{full_path}/Diffusion",
                                outpath=f"{full_path}/SeqDesign/Recycle-0",
                                args_seqdesign=args_seqdesign,
                                args_diffusion=args_diffusion,
                                relax_round=0)

    # Perform cycles of Rosetta FastRelax and Sequence Design
    if args_seqdesign["relax_design_cycles"]:
        args_seqdesign["num_seqs"] = 1
        for c in range(int(args_seqdesign["relax_design_cycles"])):
            relax_round = c+1
            outpath = relax_and_design(inpath=f"{outpath}/outputs/packed", 
                                    outpath=f"{full_path}/SeqDesign/Recycle-{str(relax_round)}",
                                    args_seqdesign=args_seqdesign, 
                                    args_diffusion=args_diffusion,
                                    model_motif=model_motif,
                                    ref_motif=ref_motif,
                                    relax_round=relax_round)
            
    concat_fasta_files(outpath=f"{full_path}/SeqDesign", name=args_diffusion["name"])

    # Print time needed
    end = time.time()
    print("Time needed: ", round(end-start, 2), "sec")


if __name__ == "__main__":
    main()