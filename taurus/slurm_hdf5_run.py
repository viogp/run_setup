''' Generate and submit SLURM jobs for running GNE (get_nebular_emission) '''
import os
import gne.gne_slurm as su

verbose = True
nvol = 64

submit_jobs = True  # False to only generate scripts
check_all_jobs = False
clean = False

# Optional: user-defined suffix for job names
# If None, suffix is derived from cutcols/mincuts/maxcuts in param_file
job_suffix = None 

# Galform in taurus
hpc = 'taurus'
sam = 'Galform'

simulations = {
    "Galform": {
        "script": "run_gne_galform.py",
        "runs": [
            ('UNIT1GPC_fnl0', [90, 87, 86, 81, 74, 65, 98, 104, 109]),
            ('UNIT1GPC_fnl100', [89, 86, 85, 80, 73, 64, 97, 103, 108]),
            ('SU1', [90, 87, 86, 81, 74, 65, 96, 98, 104, 109]),
            ('SU2', [90, 87, 86, 81, 74, 65, 96, 98, 104, 109]),
            #('UNIT1GPC_fnl0', [98, 109, 87, 90, 104]),
            #('UNIT1GPC_fnl100', [108, 103, 97, 89, 86]),
        ]
    },
    "Shark": {
        "script": "run_gne_shark.py",
        "runs": [
            ('UNIT1GPC_fnl0', [90, 87, 86, 81, 74, 65, 98, 104, 109]),
            ('UNIT1GPC_fnl100', [89, 86, 85, 80, 73, 64, 97, 103, 108]),
            ('SU1', [90, 87, 86, 81, 74, 65, 96, 98, 104, 109]),
            ('SU2', [90, 87, 86, 81, 74, 65, 96, 98, 104, 109]),
            # ('SU1', [65])
        ]
    }
}


# Parameter file to use as base
# The catalogue path, subvols and snapshot will be modified
param_file = os.path.join(os.getcwd(),simulations[sam]["script"])

# Select which runs to process
runs = simulations[sam]["runs"]
if hpc=='taurus':
    root = '/data21/users/vgonzalez/Data' 

logdir =  os.path.join(os.getcwd(),'logs')
    
# Submit, check or clean
if clean:
    su.clean_all_jobs(runs,root,sam,param_file,nvol,only_show=True,
                      logdir=logdir,job_suffix=job_suffix)
elif check_all_jobs:
    # results = su.check_all_jobs(runs, root, sam, param_file, str(nvol),
    #                             logdir=logdir,job_suffix=job_suffix,verbose=True)
    for sim, snaps in runs:
        simpath = os.path.join(root,sam,sim)
        for snap in snaps:
            results = su.check_all_jobs(sam, snap, logdir=logdir,job_suffix=job_suffix,verbose=True)
else:    
    job_count = 0
    for sim, snaps in runs:
        simpath = os.path.join(root,sam,sim)
        for snap in snaps:
            # Generate SLURM script
            script_path = su.create_slurm_script(
                hpc, param_file, simpath, sam, snap, str(nvol),
                logdir=logdir,job_suffix=job_suffix,
                verbose=verbose
            )
            if verbose: 
                print(f'  Created script: {script_path}')
                
            # Submit the job
            if submit_jobs:
                su.submit_slurm_job(script_path, verbose=verbose)
                job_count += 1
    
    if submit_jobs and verbose:
        print(f'Total jobs submitted: {job_count}')

