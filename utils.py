import subprocess
from Bio.PDB import PDBParser, Superimposer

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

# Get index of ligand residue
def get_ligand_index(pdb_file, ligand):
    parser=PDBParser(QUIET=True)
    chain = parser.get_structure('model', pdb_file)[0]["B"]
    for resi in chain.get_residues():
        if resi.get_resname() == ligand:
            return resi.get_id()[1]
    return None         


# Get Ca RMSD
def get_ca_rmsd(model_path, ref_path):
    parser=PDBParser(QUIET=True)
    model_structure = parser.get_structure('model', model_path)[0]['A']
    ref_structure = parser.get_structure('ref', ref_path)[0]['A']
    model_ca_atoms = [atom for atom in model_structure.get_atoms() if atom.get_id() == "CA"]
    ref_ca_atoms = [atom for atom in ref_structure.get_atoms() if atom.get_id() == "CA"]
    super_imposer = Superimposer()
    super_imposer.set_atoms(model_ca_atoms, ref_ca_atoms)
    return round(super_imposer.rms, 2)


# Get Ca motif RMSD
def get_motif_ca_rmsd(model_path, ref_path, model_motif, ref_motif):
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
    return round(super_imposer.rms, 2)

# Get all-atom motif RMSD
def get_motif_all_atom_rmsd(model_path, ref_path, model_motif, ref_motif):
    ref_chain = ref_motif[0][0]
    parser=PDBParser(QUIET=True)
    model_structure = parser.get_structure('model', model_path)[0]['A']
    ref_structure = parser.get_structure('ref', ref_path)[0][ref_chain]
    ref_atoms = []
    model_atoms = []
    for resi in ref_motif:
        ref_atom_list = [atom for atom in ref_structure[int(resi[1:])].get_atoms() if atom.element != 'H']
        ref_atom_list.sort(key=atom_sort_key)
        ref_atoms.extend(ref_atom_list)
    for resi in model_motif:
        model_atom_list = [atom for atom in model_structure[int(resi[1:])].get_atoms() if atom.element != 'H']
        model_atom_list.sort(key=atom_sort_key)
        model_atoms.extend(model_atom_list)
    super_imposer = Superimposer()
    super_imposer.set_atoms(model_atoms, ref_atoms)
    return round(super_imposer.rms, 2)








