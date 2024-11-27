import subprocess, copy
from Bio.PDB import PDBParser, Superimposer, PDBIO, Chain
from Bio import SeqIO
import numpy as np

# Key to sort atoms in pdb file
def atom_sort_key(atom):
    return atom.id

# Run subprocess
def run(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, text=True)
    while True:
        line = process.stdout.readline()
        if line:
            print(f"Output: {line.strip()}")
        else: 
            break    
    return_code = process.wait()

# Get fixed motif from contig string
def get_motifs(contig_str):
    design_motif = []
    ref_motif = []
    all_residues = []
    idx = 0
    for block in contig_str.split(" ")[0].split("/"):
        if block == "0":
            # Ignore chain breaks
            continue
        if block[0].isalpha():
            lb = int(block.split("-")[0][1:])
            ub = int(block.split("-")[1])
            diff = ub-lb+1
            i = 0
            while i < diff:
                idx += 1
                ref_motif.append(block[0] + str(lb+i))
                design_motif.append("A" + str(idx))
                i += 1
        else:
            idx = idx + int(block.split("-")[0])
    for i in range(idx):
        all_residues.append(f"A{i+1}")
    redesigned_residues = list(set(all_residues) - set(design_motif))
    return design_motif, ref_motif, redesigned_residues


# Get index of ligand residue
def get_ligand_index(pdb_file, ligand):
    parser=PDBParser(QUIET=True)
    chain = parser.get_structure('design', pdb_file)[0]["B"]
    for resi in chain.get_residues():
        if resi.get_resname() == ligand:
            return resi.get_id()[1]
    return None         

# Get Ca RMSD
def get_ca_rmsd(design_path, ref_path):
    parser=PDBParser(QUIET=True)
    design_structure = parser.get_structure('design', design_path)[0]['A']
    ref_structure = parser.get_structure('ref', ref_path)[0]['A']
    design_ca_atoms = [atom for atom in design_structure.get_atoms() if atom.get_id() == "CA"]
    ref_ca_atoms = [atom for atom in ref_structure.get_atoms() if atom.get_id() == "CA"]
    super_imposer = Superimposer()
    super_imposer.set_atoms(design_ca_atoms, ref_ca_atoms)
    return round(super_imposer.rms, 2)

# Get Ca motif RMSD
def get_motif_ca_rmsd(design_path, ref_path, design_motif, ref_motif):
    parser=PDBParser(QUIET=True)
    design_structure = parser.get_structure('design', design_path)[0]['A']
    ref_model = parser.get_structure('ref', ref_path)[0]
    ref_atoms = []
    design_atoms = []
    for resi in ref_motif:
        ref_chain = resi[0]
        ref_atom_list = [atom for atom in ref_model[ref_chain][int(resi[1:])].get_atoms() if atom.get_id() == 'CA']
        ref_atoms.extend(ref_atom_list)
    for resi in design_motif:
        design_atom_list = [atom for atom in design_structure[int(resi[1:])].get_atoms() if atom.get_id() == 'CA']
        design_atoms.extend(design_atom_list)
    super_imposer = Superimposer()
    super_imposer.set_atoms(design_atoms, ref_atoms)
    return round(super_imposer.rms, 2)

# Superimpose based on all-atom motif
def superimpose_motif_all_atom(design_model, ref_model, design_motif, ref_motif, apply_change):
    design_structure = design_model['A']
    ref_atoms = []
    design_atoms = []
    for resi in ref_motif:
        ref_chain = resi[0]
        ref_atom_list = [atom for atom in ref_model[ref_chain][int(resi[1:])].get_atoms() if atom.element != 'H']
        ref_atom_list.sort(key=atom_sort_key)
        ref_atoms.extend(ref_atom_list)
    for resi in design_motif:
        design_atom_list = [atom for atom in design_structure[int(resi[1:])].get_atoms() if atom.element != 'H']
        design_atom_list.sort(key=atom_sort_key)
        design_atoms.extend(design_atom_list)
    super_imposer = Superimposer()
    super_imposer.set_atoms(design_atoms, ref_atoms)
    if not apply_change:
        return round(super_imposer.rms, 2)
    else:
        super_imposer.apply(ref_model.get_atoms())
        return ref_model

# Get all-atom motif RMSD
def get_motif_all_atom_rmsd(design_path, ref_path, design_motif, ref_motif):
    parser=PDBParser(QUIET=True)
    design_model = parser.get_structure('design', design_path)[0]
    ref_model = parser.get_structure('ref', ref_path)[0]
    return superimpose_motif_all_atom(design_model=design_model,
                                      ref_model=ref_model,
                                      design_motif=design_motif,
                                      ref_motif=ref_motif,
                                      apply_change=False)

# Add sidechain coordinates for given residue using a reference structure
def add_sidechain_coordinates(design_model, ref_model, design_resi, ref_resi):
    ref_chain = ref_resi[0]
    bb_atoms = ["CA", "N", "C", "O"]
    design_atoms = [atom for atom in design_model["A"][int(design_resi[1:])].get_atoms() if atom.get_id() in bb_atoms]
    ref_atoms = [atom for atom in ref_model[ref_chain][int(ref_resi[1:])].get_atoms() if atom.get_id() in bb_atoms]
    design_atoms.sort(key=atom_sort_key)
    ref_atoms.sort(key=atom_sort_key)
    super_imposer = Superimposer()
    super_imposer.set_atoms(design_atoms, ref_atoms)
    sidechain_atoms = [atom for atom in ref_model[ref_chain][int(ref_resi[1:])].get_atoms() if (atom.get_id() not in bb_atoms) & (atom.element != 'H')]
    super_imposer.apply(sidechain_atoms)
    for atom in sidechain_atoms:
        design_model["A"][int(design_resi[1:])].add(atom)
    return design_model

# Get ligand residue from model given the name of the ligand
def get_ligand_residue(model, ligand_name):
    print("Search for ligand: ", ligand_name)
    for chain in model:
        print(chain.id)
        for residue in chain:
            if residue.resname == ligand_name:
                print(f"Ligand {ligand_name} found in Chain {chain.id}, Residue {residue.id}")
                return residue
    print(f"Ligand {ligand_name} not found in the PDB file.")
    return None

# Get ligand residue from model path given the name of the ligand
def get_ligand_residue_from_path(path, ligand_name):
    parser=PDBParser(QUIET=True)
    model = parser.get_structure('structure', path)[0]
    return get_ligand_residue(model=model, ligand_name=ligand_name)

# Add coordinates of ligand to design
def add_ligand_coordinates(design_model, ref_model, design_motif, ref_motif, ligand_name):
    ref_model = superimpose_motif_all_atom(design_model=design_model, 
                                           ref_model=ref_model, 
                                           design_motif=design_motif, 
                                           ref_motif=ref_motif, 
                                           apply_change=True)
    ligand_residue = get_ligand_residue(model=ref_model, ligand_name=ligand_name)
    new_chain = Chain.Chain("B")
    design_model.add(new_chain)
    design_model["B"].add(ligand_residue)
    return design_model

# Add sidechain and ligand coordinates to design
def add_sidechain_and_ligand_coordinates(design_path, ref_path, design_motif, ref_motif, ligand_name):
    io = PDBIO()
    parser=PDBParser(QUIET=True)
    design_model = parser.get_structure('design', design_path)[0]
    ref_model = parser.get_structure('ref', ref_path)[0]
    for i, design_resi in enumerate(design_motif):
        design_model = add_sidechain_coordinates(design_model=design_model, 
                                                 ref_model=copy.deepcopy(ref_model), 
                                                 design_resi=design_resi, 
                                                 ref_resi=ref_motif[i])
    design_model = add_ligand_coordinates(design_model=design_model,
                                          ref_model=copy.deepcopy(ref_model),
                                          design_motif=design_motif,
                                          ref_motif=ref_motif,
                                          ligand_name=ligand_name)
    io.set_structure(design_model)
    io.save(design_path)

# Remove chain from pdb file
def remove_chain_from_pdb(design_path, chain_to_remove):
    parser=PDBParser(QUIET=True)
    design_model = parser.get_structure('design', design_path)[0]
    new_model = design_model.copy()
    for chain in list(new_model):  # Use list() to avoid modifying the structure while iterating
        print(chain.id)
        if chain.id == chain_to_remove:
            print("Chain found")
            new_model.detach_child(chain_to_remove)  # Remove the chain
            break

    # Save the new structure to a PDB file
    io = PDBIO()
    io.set_structure(new_model)
    io.save(design_path)

# Get sequences from fasta_file
def get_seqs_from_fasta(fasta_file):
    seq_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_dict[record.id] = str(record.seq)
    return seq_dict

# Returns distance from the protein C-alpha to the closest ligand atom
def residue_dist_to_ligand(protein_residue, ligand_residue, atomtype='CA'):
    dist = []
    for atom in ligand_residue:
        if atomtype in protein_residue:
            vector  = protein_residue[atomtype].coord - atom.coord
            dist.append(np.sqrt(np.sum(vector * vector)))
    if len(dist) > 0:
        return min(dist)
    else:
        return None
    
# Get protein atoms located within a given distance from ligand
def get_close_protein_atoms(ligand, distance, model, atom_list):
    close_atoms = []
    chains = model.child_dict
    for c in chains:
        for protein_res in chains[c].child_list:
            if not protein_res.resname == ligand.resname:
                for atom in atom_list:
                    dist = residue_dist_to_ligand(protein_res, ligand, atom)
                    if dist and dist < distance :
                        close_atoms.append((protein_res.resname, protein_res.id[1], atom, dist))
    return close_atoms

# Get Ca atoms located within a given distance from ligand
def get_close_ca_atoms(ligand, distance, design_path):
    parser=PDBParser(QUIET=True)
    design_model = parser.get_structure('design', design_path)[0]
    close_ca_atoms = get_close_protein_atoms(ligand=ligand, distance=distance, model=model, atom_list=["CA"])
    return close_ca_atoms


# Returns backbone atoms located within a given distance from ligand
def get_close_backbone_atoms(ligand, distance, design_path):
    parser=PDBParser(QUIET=True)
    design_model = parser.get_structure('design', design_path)[0]
    close_backbone_atoms = get_close_protein_atoms(ligand=ligand, distance=distance, model=design_model, atom_list=["CA","C","O","N"])
    return close_backbone_atoms


