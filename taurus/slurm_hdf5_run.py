''' Generate and submit SLURM jobs for running GNE (get_nebular_emission) '''
import os
import gne.gne_slurm as su

verbose = True
nvol = 64

submit_jobs = True  # False to only generate scripts
check_all_jobs = True
clean = False

# Optional: user-defined suffix for job names
# If None, suffix is derived from cutcols/mincuts/maxcuts in param_file
job_suffix = None 

# Galform in taurus
hpc = 'taurus'
sam = 'Shark'

simulations = {
    "Galform": {
        "script": "run_gne_SU1.py",
        "runs": [
            ('SU1', [109, 104, 98, 90, 87, 128, 96, 78]),
            ('SU1', [128, 90, 87, 96, 78]),
            ('SU2', [90]),
            ('UNIT1GPC_fnl0', [128, 90]),
            ('UNIT1GPC_fnl100', [127, 89]),
            #('SU1', [109, 104, 98, 90, 96, 78]),
            #('SU2', [109, 104, 98, 90, 87]),
            #('UNIT1GPC_fnl0', [98, 109, 87, 90, 104]),
            #('UNIT1GPC_fnl100', [108, 103, 97, 89, 86]),
        ]
    },
    "Shark": {
        "script": "run_gne_shark.py",
        "runs": [
            ('SU1', [128, 90, 87, 96, 78]),
            ('SU2', [128, 90]),
            ('UNIT1GPC_fnl0', [128, 90]),
            ('UNIT1GPC_fnl100', [127, 89]),
        ]
    }
}


# Parameter file to use as base
# The catalogue path, subvols and snapshot will be modified
param_file = os.path.join(os.getcwd(),simulations[sam]["script"])

# Select which runs to process
runs = simulations[sam]["runs"]
if hpc=='taurus':
    root = '/home2/vgonzalez/Data' 

logdir =  os.path.join(os.getcwd(),'logs')
    
# Submit, check or clean
if clean:
    su.clean_all_jobs(runs,root,sam,param_file,nvol,only_show=True,
                      logdir=logdir,job_suffix=job_suffix)
elif check_all_jobs:
    results = su.check_all_jobs(runs, root, sam, param_file, nvol,
                                logdir=logdir,job_suffix=job_suffix,verbose=True)
else:    
    job_count = 0
    for sim, snaps in runs:
        simpath = os.path.join(root,sam,sim)
        for snap in snaps:
            # Generate SLURM script
            script_path, job_name = su.create_slurm_script(
                hpc, param_file, simpath, snap, nvol,
                logdir=logdir,job_suffix=job_suffix,
                verbose=verbose
            )
            if verbose: 
                print(f'  Created script: {script_path}')
                
            # Submit the job
            if submit_jobs:
                su.submit_slurm_job(script_path, job_name)
                job_count += 1
    
    if submit_jobs and verbose:
        print(f'Total jobs submitted: {job_count}')

