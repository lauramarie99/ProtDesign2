import utils, glob, yaml, argparse, json, os, config, time
import pyrosetta as pyr
from Bio.PDB import PDBParser, Superimposer
from concurrent.futures import ProcessPoolExecutor, as_completed
start = time.time()

# Creates a json file given a dictionary
def create_json(outpath, dict):
    with open(outpath, 'w') as file:
        json.dump(dict, file)

# Filter pdbs based on Ca motif RMSD
def filter_pdb(model_path, ref_path, model_motif, ref_motif, threshold):
    ref_chain = ref_motif[0][0]
    parser=PDBParser(QUIET=True)
    model_structure = parser.get_structure('model', model_path)[0]['A']
    ref_structure = parser.get_structure('ref', ref_path)[0][ref_chain]
    ref_atoms = []
    model_atoms = []
    for resi in ref_motif:
        ref_atom_list = [atom for atom in ref_structure[int(resi[1:])].get_atoms() if atom.get_id() == 'CA']
        ref_atoms.extend(ref_atom_list)
    for resi in model_motif:
        model_atom_list = [atom for atom in model_structure[int(resi[1:])].get_atoms() if atom.get_id() == 'CA']
        model_atoms.extend(model_atom_list)
    super_imposer = Superimposer()
    super_imposer.set_atoms(model_atoms, ref_atoms)
    print(model_path, super_imposer.rms)
    if super_imposer.rms <= threshold:
        return True
    return False


# Get fixed motif from contig string
def get_motifs(contig_str):
    model_motif = []
    ref_motif = []
    index = 0
    for block in contig_str.split("/"):
        if block[0].isalpha():
            lb = int(block.split("-")[0][1:])
            ub = int(block.split("-")[1])
            range = ub-lb+1
            i = 0
            while i < range:
                index += 1
                ref_motif.append(block[0] + str(lb+i))
                model_motif.append(block[0] + str(index))
                i += 1
        else:
            index = index + int(block.split("-")[0])
    return model_motif, ref_motif


# Filter PDB files and create the json input files for sequence design
def preprocessing(inpath, name, outpath, contig_str, ref_path, threshold):
    pdb_files = glob.glob(f"{inpath}/{name}*.pdb")                                  # Get all pdb file paths
    model_motif, ref_motif = get_motifs(contig_str)                                 # Get motifs
    filtered_pdb_files = []                                                         # Filter pdb files based on Motif Ca-RMSD
    for path in pdb_files:
        if filter_pdb(model_path=path,ref_path=ref_path,model_motif=model_motif,ref_motif=ref_motif,threshold=threshold):
            filtered_pdb_files.append(path)
    print("Filtered PDB files: ", filtered_pdb_files)
    pdb_dict = {path: "" for path in filtered_pdb_files}                            # Create dictionaries with paths and motif
    fixed_resi_str = " ".join(model_motif)
    resi_dict = {path: fixed_resi_str for path in filtered_pdb_files}
    
    os.makedirs(outpath, exist_ok=True)
    create_json(f"{outpath}/pdb_ids.json", pdb_dict)                                # Create json input files
    create_json(f"{outpath}/fix_residues_multi.json", resi_dict)
    return model_motif


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
            if "n" in id:
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
    model_motif = preprocessing(inpath=inpath, 
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
    return outpath, model_motif


# Rosetta FastRelax
def relax(pdb, id, outpath, fixed_resi, params=None, cst=None):
    # Initialization
    extra_res_fa = f"-extra_res_fa {params}" if params is not None else ""
    constraints = f"-constraints:cst_fa_file {cst}" if cst is not None else ""
    pyr.init(f"-ignore_zero_occupancy false -ex1 -ex2 {extra_res_fa} {constraints}")

    # Define scorefunction
    scorefxn = pyr.get_fa_scorefxn()                                                                    # Default full-atom energy terms
    scorefxn.set_weight(pyr.rosetta.core.scoring.score_type_from_name("atom_pair_constraint"), 1)       # Add constraint weight to score function    

    # Read input pdb file
    pose = pyr.pose_from_file(pdb)
    pose2 = pose.clone()

    # Add constraints
    if cst is not None:
        constraint_mover = pyr.rosetta.protocols.constraint_movers.ConstraintSetMover()
        constraint_mover.constraint_file(cst)
        constraint_mover.apply(pose2)

    # Specify task operations
    tf = pyr.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyr.rosetta.core.pack.task.operation.InitializeFromCommandline())     # Use options specified in the init() function
    tf.push_back(pyr.rosetta.core.pack.task.operation.RestrictToRepacking())           # Avoid design of residues
    tf.push_back(pyr.rosetta.core.pack.task.operation.IncludeCurrent())                # Includes the current rotamers to the rotamers sets   
    tf.push_back(pyr.rosetta.core.pack.task.operation.NoRepackDisulfides())            # Prevent repacking of disulfid bond side-chains   
    residue_selector = pyr.rosetta.core.select.residue_selector.ResidueIndexSelector() # Prevent repacking of catalytic residues 
    for resi in fixed_resi:
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

    # Perform FastRelax
    fastRelax = pyr.rosetta.protocols.relax.FastRelax()
    fastRelax.constrain_relax_to_start_coords(True)
    fastRelax.set_scorefxn(scorefxn)
    fastRelax.set_movemap_factory(mmf)
    fastRelax.set_task_factory(tf)
    #fastRelax.constrain_coords(True)
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
def relax_pdb(pdb, relax_round, outpath, args_seqdesign, fixed_resi):
    id = pdb.split("/")[-1].split("_packed")[0]
    if "c" in id:
        index = id.index("c")
        design_id = id[:index+1] + str(relax_round) + id[index+2:]
    else:
        design_id = id + "_n" + pdb.split("/")[-1].split("packed_")[1][:-6] + "_c" + str(relax_round)
    outpath = relax(pdb=pdb, 
                    id=design_id, 
                    outpath=outpath, 
                    **({"params": args_seqdesign["params_file"]} if "params_file" in args_seqdesign else {}),
                    **({"cst": args_seqdesign["cst_file"]} if "cst_file" in args_seqdesign else {}),
                    fixed_resi=fixed_resi)


# Perform cycles of Rosetta FastRelax and Sequence Design
def relax_and_design(inpath, outpath, args_seqdesign, args_diffusion, relax_round, fixed_resi):
    pdb_files = glob.glob(f"{inpath}/*.pdb")                                                            # Extract all pdb files
    os.makedirs(os.path.join(outpath, "relaxed"), exist_ok=True)                                        # Create outdir
    
    # Relaxation (in batches)
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(relax_pdb, pdb, relax_round, outpath, args_seqdesign, fixed_resi) for pdb in pdb_files]
        for future in as_completed(futures):
            future.result()

    # Proceed with the design function after all relax calls are complete
    outpath, model_motif = design(inpath=f"{outpath}/relaxed", 
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
    outpath, model_motif = design(inpath=f"{full_path}/Diffusion",
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
                                    fixed_resi=model_motif,
                                    relax_round=relax_round)
            
    concat_fasta_files(outpath=f"{full_path}/SeqDesign", name=args_diffusion["name"])

    # Print time needed
    end = time.time()
    print("Time needed: ", round(end-start, 2), "sec")


if __name__ == "__main__":
    main()