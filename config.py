# Define paths to repositories
RFDIFFUSION_PATH = "RFdiffusion"
RFDIFFUSIONAA_PATH = "rf_diffusion_all_atom"
LIGANDMPNN_PATH = "LigandMPNN"

# To run on cluster
REPO_PATH = ""                              # Path to ProtDesign2 Repo
RFDIFFUSION_CONTAINER = ""                  # Paths to singularity container
RFDIFFUSIONAA_CONTAINER = ""
SEQDESIGN_CONTAINER = ""
FOLDING_CONTAINER = ""
EMAIL = ""                                  # Email for slurm notifications
EMAIL_TYPE = "ALL"
PARTITION = "paula"                         # Partition to use for runs
GPU = "a30:1"                               # GPU to use
COLABFOLD_WEIGHTS = ""                      # Path to colabfold weights

