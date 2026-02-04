''' Generate and submit SLURM jobs for running GNE (get_nebular_emission) '''
from gne.gne_slurm as su

verbose = True
nvol = 2

submit_jobs = False  # False to only generate scripts
check_all_jobs = True
clean = False

# Parameter file to use as base (will modify path, subvols and snapshot in root)
param_files = 'run_gne_SU1.py'

# Optional: user-defined suffix for job names
# If None, suffix is derived from cutcols/mincuts/maxcuts in param_file
job_suffix = None 

# Galform in taurus
hpc = 'taurus'
sam = 'Galform'
taurus_galform = [
    #('SU1', [109, 104, 98, 90, 87, 128, 96, 78]),
    ('SU1', [87, 128]),
    #('SU1', [109, 104, 98, 90, 96, 78]),
    ('SU2', [109, 104, 98, 90, 87]),
    ('UNIT1GPC_fnl0', [98, 109, 87, 90, 104]),
    ('UNIT1GPC_fnl100', [108, 103, 97, 89, 86]),
]

# Select which runs to process
runs = taurus_galform
if hpc=='taurus':
    root = '/home2/vgonzalez/Data' 

# Submit, check or clean
if clean:
    su.clean_all_jobs(runs, only_show=True)
elif check_all_jobs:
    results = su.check_all_jobs(runs,verbose=True)
else:    
    job_count = 0
    for sim, snaps in runs:
        simpath = os.path.join(root,sam,sim)
        print(simpath)
#        for snlap in snaps:
#            # Generate SLURM script
#            script_path, job_name = su.create_slurm_script(
#                hpc, param_file, simpath, snap, nvol,
#                verbose=verbose, job_suffix=job_suffix
#            )
#            if verbose: 
#                print(f'  Created script: {script_path}')
#                
#            # Submit the job
#            if submit_jobs:
#                su.submit_slurm_job(script_path, job_name)
#                job_count += 1
#    
#    if submit_jobs and verbose:
#        print(f'Total jobs submitted: {job_count}')
#
