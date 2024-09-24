import subprocess, glob, os, argparse, config, yaml

# Create a slurm script
def create_slurm_script(repo_path, 
                        container, 
                        config_file, 
                        script, 
                        out_path, 
                        name, 
                        jobname, 
                        time, 
                        mem, 
                        cpus, 
                        gpu, 
                        partition,
                        email,
                        emailType,
                        dependency=None):
    with open(f'{out_path}/{name}.slurm', 'w') as slurmFile:
        slurmFile.writelines([
                "#!/bin/bash\n",
                f"#SBATCH --job-name={jobname}\n",
                f"#SBATCH --output={out_path}/{name}.out\n",
                f"#SBATCH --error={out_path}/{name}.err\n",
                f"#SBATCH --time={time}\n",
                f"#SBATCH --mem={mem}\n",
                f"#SBATCH --cpus-per-task={cpus}\n",
                f"#SBATCH --gres=gpu:{gpu}\n",
                f"#SBATCH --partition={partition}\n",
                f"#SBATCH --mail-user={email}\n",
                f"#SBATCH --mail-type={emailType}\n"
        ])
        if dependency:
            slurmFile.writelines([
                f"#SBATCH --dependency=afterok:{dependency}\n"
            ])
        slurmFile.writelines([
                '# define CONTAINER\n',
                f'CONTAINER={container}\n',
                '# define SCRIPT or program to call inside the container\n',
                f'SCRIPT="{script} --config {config_file}"\n',
                f'cd {repo_path}\n',
                f'singularity exec --nv -B {config.COLABFOLD_WEIGHTS}:/cache --cleanenv $CONTAINER $SCRIPT\n'      
        ])

# Run single slurm script, returns job id and error message
def run_slurm_script(name, cwd):
    slurm_command = f'/usr/bin/sbatch {name}.slurm'
    process = subprocess.Popen(slurm_command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               universal_newlines=True,
                               shell=True,
                               cwd=cwd)
    output, error = process.communicate()
    if error == '':
        output = output.split()[-1]
    return output, error

# Run all slurm scripts, returns dictionary with job ids and dictionary with error messages
def run_all_slurm_scripts(slurm_path):
    slurm_files = glob.glob(
        f'{slurm_path}/*.slurm')
    job_ids = {}
    errors = {}
    for slurm_file in slurm_files:
        name = slurm_file.split('/')[-1].split('.')[0]
        output, error = run_slurm_script(name=name, cwd=slurm_path)
        job_ids[name] = output
        errors[name] = error
    return job_ids, errors

    
# Start diffusion jobs, returns dictionary with job ids and dictionary with error messages   
def run_diffusion(repo_path, config_path, out_path, diffusion_cmd):
    config_files = glob.glob(f'{config_path}/*.yml')
    for config_file in config_files:
        args = yaml.safe_load(open(config_file))
        if args["diffusion"]["type"] == "all-atom":
            container = config.RFDIFFUSIONAA_CONTAINER
        else:
            container = config.RFDIFFUSION_CONTAINER
        name = args["diffusion"]["name"]
        create_slurm_script(repo_path=repo_path, 
                            container=container,
                            config_file=config_file, 
                            script=diffusion_cmd, 
                            out_path=out_path, 
                            name=name,
                            jobname=f"diff-{name}", 
                            time="01:00:00", 
                            mem="4000", 
                            cpus=1, 
                            gpu=config.GPU, 
                            partition=config.PARTITION,
                            email=config.EMAIL,
                            emailType=config.EMAIL_TYPE)
    job_ids, errors = run_all_slurm_scripts(slurm_path=out_path)
    return job_ids, errors


# Start validation job if diffusion is done, returns dictionary with job ids and dictionary with error messages 
def run_validation(repo_path, config_path, out_path, dep_job_ids, container, cmd, stage, time, cpus):
    for name,job_id in dep_job_ids.items():
        config_file = f"{config_path}/{name}.yml"
        create_slurm_script(repo_path=repo_path, 
                            container=container,
                            config_file=config_file, 
                            script=cmd, 
                            out_path=out_path, 
                            name=name,
                            jobname=f"{stage}-{name}", 
                            time=time, 
                            mem="5000", 
                            cpus=cpus, 
                            gpu=config.GPU,
                            partition=config.PARTITION, 
                            email=config.EMAIL,
                            emailType=config.EMAIL_TYPE,
                            dependency=job_id)
    job_ids, errors = run_all_slurm_scripts(slurm_path=out_path)
    return job_ids, errors


"""
MAIN
"""

# Global variables
argParser = argparse.ArgumentParser()
argParser.add_argument('-i','--input')                                      # Input path that contains Config folder
args = argParser.parse_args()

config_path = f"{args.input}/Configs"                                       # Location config files
out_path = f"{args.input}/Slurm"                                            # Location slurm files

# Run diffusion and validation
os.makedirs(f"{out_path}/Diffusion", exist_ok=True)
os.makedirs(f"{out_path}/SeqDesign", exist_ok=True)
os.makedirs(f"{out_path}/Folding", exist_ok=True)

diffusion_job_ids, diffusion_errors = run_diffusion(repo_path=config.REPO_PATH, 
                                                    config_path=config_path,
                                                    out_path=f"{out_path}/Diffusion",
                                                    diffusion_cmd="python3.9 diffuse.py")
print("Diffusion jobs submitted")

seqdesign_job_ids, seqdesign_errors = run_validation(repo_path=config.REPO_PATH,
                                                     config_path=config_path,
                                                     out_path=f"{out_path}/SeqDesign",
                                                     dep_job_ids=diffusion_job_ids,
                                                     container=config.SEQDESIGN_CONTAINER,
                                                     cmd="python3.9 seqdesign.py",
                                                     stage="seqdesign",
                                                     time="03:00:00",
                                                     cpus=9)
print("Seqdesign jobs submitted")

folding_job_ids, folding_errors = run_validation(repo_path=config.REPO_PATH,
                                                       config_path=config_path,
                                                       out_path=f"{out_path}/Folding",
                                                       dep_job_ids=seqdesign_job_ids,
                                                       container=config.FOLDING_CONTAINER,
                                                       cmd="python fold.py",
                                                       stage="folding",
                                                       time="01:00:00",
                                                       cpus=1)
print("Colabfold jobs submitted")




