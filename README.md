# ProtDesign2
Protein design with RFdiffusion(all-atom) and ProteinMPNN/LigandMPNN/RosettaFastRelax.\
Create new proteins such as enzymes or small-molecule binders.

## Overview
- diffuse.py: Diffusion (RFdiffusion/RFdiffusion all-atom)
- seqdesign.py: SequenceDesign (ProteinMPNN/LigandMPNN and Rosetta FastRelax)
- fold.py: Forward Folding (ColabFold)
- run_cluster.py: Slurm script generation and job submission to cluster

## Setup

- Clone this repository
- Setup RFdiffusion, LigandMPNN/ProteinMPNN and ColabFold

### RFdiffusion
- Use RFdiffusion_dockerfile to build a singularity/docker container
- Clone official repo
- Get weights
```
cd RFdiffusion
mkdir models && cd models
wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
```

### RFdiffusion all-atom
- Clone official repo
- Download the official container
```
wget http://files.ipd.uw.edu/pub/RF-All-Atom/containers/rf_se3_diffusion.sif
```
- Get weights
```
cd rf_diffusion_all_atom
wget http://files.ipd.uw.edu/pub/RF-All-Atom/weights/RFDiffusionAA_paper_weights.pt
```
- Initialize git submodules
```
git submodule init
git submodule update
```

### LigandMPNN/ProteinMPNN
- Clone official repo
- Use LigandMPNN_Rosetta_dockerfile to build a singularity/docker container
- Get weights
```
cd LigandMPNN
bash get_model_params.sh "./model_params"
```

### ColabFold
- Get official container
```
singularity pull docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2
```
- Get weights
```
mkdir cache
singularity run -B /local/path/to/cache:/cache \
  colabfold_1.5.5-cuda12.2.2.sif \
  python -m colabfold.download
```

### General
- Adapt paths in config.py


## Get started
For the diffusion and validation (ProteinMPNN + AF) steps, only one single config file is used.

### Diffusion
- type: "base" (RFdiffusion) or "all-atom" (RFdiffusion all-atom)
- ckpt_override_path: Override RFdiffusion model path (e.g. Active_Site model)
- contigs: Contig string (specify always ranges, e.g. 16-16 instead of 16!)
- enzyme_design: Set true if you want to use an external potential
- guide_potentials: External potential to use (only used if enzyme_design = true)
- guide_scale: Scale factor for guide potential (only used if enzyme_design = true)
- pdb: Input structure (the structure where the fixed residues are taken from)
- substrate: Substrate name (must be contained in input structure pdb)
- iterations: Number of RFdiffusion steps
- name: Experiment name (don't use special characters)
- noise scale: RFdiffusion noise scale
- num_designs: Number of designs to generate
- path: Directory where to store results

### Sequence Design
- model_type: Sequence design model (e.g. "protein_mpnn")
- num_seqs: Number of sequences to generate
- rmsd_cutoff: Threshold for motif-Ca-RMSD (if a backbone has a motif-Ca-RMSD > rmsd_threshold, it will not be considered for sequence design)
- clash_cutoff: Distance threshold for clashes (if a design has backbone atoms within the given distance, it will not be considered for sequence design)
- relax_design_cycles: Number of FastRelax+SeqDesign cycles, 0 will only carry out sequence design without FastRelax
- params_file: Params file for Rosetta FastRelax (only used if relax_design_cycles > 0)
- cst_file: CST file for Rosetta FastRelax (only used if relax_design_cycles > 0)

### Folding
- num_recycles: Number of AF2 recycles
- num_models: Number of AF2 models
- msa_mode: single_sequence

### Run diffusion
```
# RFdiffusion
singularity exec --nv /media/data/Container/RFdiffusionAA/rf_se3_diffusion.sif python3.9 diffuse.py --config config.yml
# RFdiffusion all-atom
singularity exec --nv /media/data/Container/RFdiffusion/rfdiffusion.sif python3.9 diffuse.py --config config.yml        
```

### Run sequence design
```
singularity exec --nv /media/data/Container/LigandMPNN/ligandmpnn_rosetta.sif python3.9 seqdesign.py --config config.yml
```

### Run forward folding
```
singularity exec --nv -B /home/iwe28/ProteinFolding:/cache /media/data/Container/ColabFold/colabfold_1.5.5-cuda12.2.2.sif python3 fold.py --config config.yml
```

## Large scale studies
For generation of many config files based on a general config file, the script create_configs.py in the folder configs can be used.
To automatically generate slurm scripts and submit the jobs, the script run_cluster.py can be used.

## Notes
- Avoid chain breaks in the contig strings
- This pipeline is not suitable for protein binder design
