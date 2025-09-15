import os



# dirs = ['/scratch/b/babul/renierht/output_directory/']*2 + ['/project/b/babul/renierht/']*4
# dirs = ['/scratch/b/babul/renierht/output_directory/', '/project/b/babul/renierht/']
# dirs = ['/project/b/babul/renierht/']
# dirs = ['/scratch/b/babul/renierht/output_directory/']
dirs = ['/project/b/babul/aspadawe/sims/']

# sims_ = ['Simba_N512L50', 'Kobayashi_N512L50', 'Simba_NoAGN_N512L50',
#          'Simba_NoFB_N512L50', 'Simba_NoJet_N512L50', 'Kobayashi_NoDust_N512L50']
# sims_ = ['Kobayashi_N512L50', 'Kobayashi_NoDust_N512L50']
# sims_ = ['Simba_NoJet_N512L50']
# sims_ = ['Kobayashi_N512L50']
sims_ = ['Simba_zdiff_N512L50']

# codes_ = ['simba', 'simba-c'] + ['simba']*3 + ['simba-c']
# codes_ = ['simba-c', 'simba-c']
codes_ = ['simba']
# codes_ = ['simba-c']

# save_codes_ = ['simba', 'simba-c'] + ['simba']*3 + ['simba-c']
# save_codes_ = ['simba-c', 'simba-c']
# save_codes_ = ['simba']
# save_codes_ = ['simba-c']
save_codes_ = ['simba_zdiff/profiles/']#, 'obsidian/profiles/]

num_sims = len(sims_)




# wall_times_ = ['2:59:0']*num_sims
wall_times_ = ['8:0:0']*num_sims

nthreads_ = [40] * num_sims

# filters_ = ['Cylinder'] * num_sims
filters_ = ['Sphere'] * num_sims

profile_types_ = ['Profile'] * num_sims

# ndims_ = [2] * num_sims
ndims_ = [3] * num_sims

weight_bys_ = ['mass'] * num_sims

xscales_ = ['R500'] * num_sims
# xscales_ = ['Physical'] * num_sims

temp_cuts_ = ['5e5 K'] * num_sims

nh_cuts_ = ['0.13 cm**-3'] * num_sims

n_lgals_ = [3] * num_sims

log_m500_cuts_ = [13.0] * num_sims







# snapshot_numbers = ['151', '125', '105', '078']
# snapshot_numbers = ['151']
snapshot_numbers = ['169']

redshift_strings = ['0.000', '0.490', '0.993', '2.025']

redshift_vals = [0, 0.49, 0.993, 2.025]

prop_file_redshifts = ['0', '05', '1', '2']


codes = []
sims = []
snapshots = []
ahf_basenames = []
prop_files = []
save_files = []

redshifts = []

wall_times = []
nthreads = []
filters = []
profile_types = []
ndims = []
weight_bys = []
xscales = []
temp_cuts = []
nh_cuts = []
n_lgals = []
log_m500_cuts = []

# for dir_, sim, code, wall_time in zip(dirs, sims, codes_):
for jj in range(num_sims):
    for snap_num, z_str, z_val, prop_z in zip(snapshot_numbers, redshift_strings, redshift_vals, prop_file_redshifts):
        codes.append(codes_[jj])
        sims.append(sims_[jj])
        
#         print(os.path.join(dirs[jj], sims_[jj], 'snapshot_%s.hdf5' % snap_num))
        snapshots.append(os.path.join(dirs[jj], sims_[jj], 'snapshot_%s.hdf5' % snap_num))
        
        ahf_basenames.append(os.path.join(dirs[jj], sims_[jj], 'snapshot_%s.hdf5.z%s.AHF_' % (snap_num, z_str)))
        
        if 'halo_z%s_NW.hdf5' % prop_z in os.listdir(os.path.join(dirs[jj], sims_[jj], 'halo_results')):
            prop_files.append(os.path.join(dirs[jj], sims_[jj], 'halo_results', 'halo_z%s_NW.hdf5' % prop_z))
        else:
            prop_files.append(os.path.join(dirs[jj], sims_[jj], 'halo_results', 'halo_z%s.hdf5' % prop_z))
            
        name = 'snapshot_%s_z%s_%s_igrm_%sd_%s_profiles_Nlgal_%s_windsremoved_xscale_%s_intermediatebins' % (snap_num, z_str, sims_[jj], ndims_[jj], filters_[jj], n_lgals_[jj], xscales_[jj])
        save_files.append(os.path.join('/scratch/b/babul/aspadawe/analysis/simba/', save_codes_[jj], 'final', name))
        
        redshifts.append(z_val)
        wall_times.append(wall_times_[jj])
        nthreads.append(nthreads_[jj])
        filters.append(filters_[jj])
        profile_types.append(profile_types_[jj])
        ndims.append(ndims_[jj])
        weight_bys.append(weight_bys_[jj])
        xscales.append(xscales_[jj])
        temp_cuts.append(temp_cuts_[jj])
        nh_cuts.append(nh_cuts_[jj])
        n_lgals.append(n_lgals_[jj])
        log_m500_cuts.append(log_m500_cuts_[jj])
        
        
        
num_jobs = len(codes)
# num_jobs = 1

        
        
# snapshots = ['/scratch/b/babul/renierht/output_directory/Simba_N512L50/snapshot_151.hdf5', '/project/b/babul/renierht/Simba_NoAGN_N512L50/snapshot_151.hdf5', '/project/b/babul/renierht/Simba_NoFB_N512L50/snapshot_151.hdf5', '/project/b/babul/renierht/Simba_NoJet_N512L50/snapshot_151.hdf5', '/scratch/b/babul/renierht/output_directory/Kobayashi_N512L50/snapshot_151.hdf5', '/project/b/babul/renierht/Kobayashi_NoDust_N512L50/snapshot_151.hdf5']

# ahf_basenames = ['/scratch/b/babul/renierht/output_directory/Simba_N512L50/snapshot_151.hdf5.z0.000.AHF_', '/project/b/babul/renierht/Simba_NoAGN_N512L50/snapshot_151.hdf5.z0.000.AHF_', '/project/b/babul/renierht/Simba_NoFB_N512L50/snapshot_151.hdf5.z0.000.AHF_', '/project/b/babul/renierht/Simba_NoJet_N512L50/snapshot_151.hdf5.z0.000.AHF_', '/scratch/b/babul/renierht/output_directory/Kobayashi_N512L50/snapshot_151.hdf5.z0.000.AHF_', '/project/b/babul/renierht/Kobayashi_NoDust_N512L50/snapshot_151.hdf5.z0.000.AHF_']

# prop_files = ['/scratch/b/babul/renierht/output_directory/Simba_N512L50/halo_results/halo_z0_NW.hdf5', '/project/b/babul/renierht/Simba_NoAGN_N512L50/halo_results/halo_z0_NW.hdf5', '/project/b/babul/renierht/Simba_NoFB_N512L50/halo_results/halo_z0_NW.hdf5', '/project/b/babul/renierht/Simba_NoJet_N512L50/halo_results/halo_z0_NW.hdf5', '/scratch/b/babul/renierht/output_directory/Kobayashi_N512L50/halo_results/halo_z0_NW.hdf5', '/project/b/babul/renierht/Kobayashi_NoDust_N512L50/halo_results/halo_z0_NW.hdf5']

# # save_files = ['/scratch/b/babul/aspadawe/analysis/simba/simba/final/snapshot_151_Simba_N512L50_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins', '/scratch/b/babul/aspadawe/analysis/simba/simba/final/snapshot_151_Simba_NoAGN_N512L50_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins', '/scratch/b/babul/aspadawe/analysis/simba/simba/final/snapshot_151_Simba_NoFB_N512L50_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins', '/scratch/b/babul/aspadawe/analysis/simba/simba/final/snapshot_151_Simba_NoJet_N512L50_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins', '/scratch/b/babul/aspadawe/analysis/simba/simba-c/final/snapshot_151_Simba-C_Kobayashi_Dust_N512L50_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins', '/scratch/b/babul/aspadawe/analysis/simba/simba-c/final/snapshot_151_Simba-C_NoDust_N512L50_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins']
# save_files = ['/scratch/b/babul/aspadawe/analysis/simba/simba/final/snapshot_151_Simba_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins', '/scratch/b/babul/aspadawe/analysis/simba/simba/final/snapshot_151_Simba_NoAGN_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins', '/scratch/b/babul/aspadawe/analysis/simba/simba/final/snapshot_151_Simba_NoFB_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins', '/scratch/b/babul/aspadawe/analysis/simba/simba/final/snapshot_151_Simba_NoJet_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins', '/scratch/b/babul/aspadawe/analysis/simba/simba-c/final/snapshot_151_Simba-C_Kobayashi_Dust_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins', '/scratch/b/babul/aspadawe/analysis/simba/simba-c/final/snapshot_151_Simba-C_NoDust_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins']


# wall_times = ['2:59:0']*num_jobs

# nthreads = [40] * num_jobs

# # filters = ['Cylinder'] * num_jobs
# filters = ['Sphere'] * num_jobs

# profile_types = ['Profile'] * num_jobs

# # ndims = [2] * num_jobs
# ndims = [3] * num_jobs

# weight_bys = ['mass'] * num_jobs

# xscales = ['R500'] * num_jobs

# temp_cuts = ['5e5 K'] * num_jobs

# nh_cuts = ['0.13 cm**-3'] * num_jobs

# n_lgals = [3] * num_jobs

# log_m500_cuts = [13.0] * num_jobs



with open('job_gen_nd_profiles_igrm_template.sh', 'r') as f:
    job_data = f.read()

job_dir = './job_scripts_for_profiles'
if not os.path.isdir(job_dir):
    os.makedirs(job_dir, mode=0o755, exist_ok=True)


## Create job file, substitute in parameters, submit, (delete job file?)
new_jobs = []
for ii in range(num_jobs):
    job_file = 'job_gen_nd_profiles_igrm_%s.sh' % ii
    param_dict = {
        'CODE':codes[ii],
        'REDSHIFT':redshifts[ii],
        'WALL_TIME':wall_times[ii],
        'SNAPSHOT':snapshots[ii],
        'AHF_BASENAME':ahf_basenames[ii],
        'NTHREADS':nthreads[ii],
        'PROP_FILE':prop_files[ii],
        'SAVE_FILE':save_files[ii],
        'FILTER':filters[ii],
        'PROFILE_TYPE':profile_types[ii],
        'NDIM':ndims[ii],
        'WEIGHT_BY':weight_bys[ii],
        'XSCALE':xscales[ii],
        'TEMP_CUT':temp_cuts[ii],
        'NH_CUT':nh_cuts[ii],
        'N_LGAL':n_lgals[ii],
        'LOG_M500_CUT':log_m500_cuts[ii],
    }
    
    ## Create job_name
    job_name = '%s-z%s-%sd_%s_%s_%s' % (sims[ii], param_dict['REDSHIFT'], param_dict['NDIM'], 
                                     param_dict['FILTER'], param_dict['XSCALE'], param_dict['PROFILE_TYPE'])
#     job_name = param_dict['SAVE_FILE']
    param_dict['JOB_NAME'] = job_name
    print('job_name:', job_name)
    
    new_jobs.append(job_data)
    for param, val in param_dict.items():
        print('param:', param, ':', val)
        new_jobs[-1] = new_jobs[-1].replace(param, str(val))
        
    with open(os.path.join(job_dir, job_file), 'w') as f:
        f.write(new_jobs[-1])
    
    os.system('cd %s && sbatch %s' % (job_dir, job_file))
    
    print()