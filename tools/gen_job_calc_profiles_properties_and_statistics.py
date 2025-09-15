import os
# import numpy as np
# from glob import glob
# from pathlib import Path
# import argparse as ap


# dirs_ = ['/scratch/b/babul/aspadawe/analysis/simba/simba/final']*4 + ['/scratch/b/babul/aspadawe/analysis/simba/simba-c/final']*2
dirs_ = ['/scratch/b/babul/aspadawe/analysis/simba/simba/final', '/scratch/b/babul/aspadawe/analysis/simba/simba-c/final'] + ['/scratch/b/babul/aspadawe/analysis/simba/simba/final']*3 + ['/scratch/b/babul/aspadawe/analysis/simba/simba-c/final']
# dirs_ = ['/scratch/b/babul/aspadawe/analysis/simba/simba-c/final']
dirs_ = ['/scratch/b/babul/aspadawe/analysis/simba/simba/final', '/scratch/b/babul/aspadawe/analysis/simba/simba-c/final']

# sims_ = ['Simba_N512L50', 'Kobayashi_N512L50', 'Simba_NoAGN_N512L50',
#          'Simba_NoFB_N512L50', 'Simba_NoJet_N512L50', 'Kobayashi_NoDust_N512L50']
sims_ = ['Simba', 'Kobayashi', 'Simba_NoAGN', 'Simba_NoFB', 'Simba_NoJet', 'Kobayashi_NoDust']
# sims_ = ['Kobayashi']
sims_ = ['Simba', 'Simba-C_Kobayashi_Dust']

# codes_ = ['simba', 'simba-c'] + ['simba']*3 + ['simba-c']

num_sims = len(sims_)



filters_ = ['Cylinder'] * num_sims
# filters_ = ['Sphere'] * num_sims

ndims_ = [2] * num_sims
# ndims_ = [3] * num_sims

xscales_ = ['R500'] * num_sims
# xscales_ = ['Physical'] * num_sims

temp_cuts_ = ['5e5 K'] * num_sims

nh_cuts_ = ['0.13 cm**-3'] * num_sims

n_lgals_ = [3] * num_sims



wall_times_ = ['10:0:0']*num_sims
# wall_times_ = ['2:59:0']*num_sims

# suffixes_ = ['median_no_errors-metal_props-CC_classifications-v2']*num_sims
# suffixes_ = ['median_no_errors-thermo_props-CC_classifications-v2']*num_sims
# suffixes_ = ['median_no_errors-metal_props-CC_classifications']*num_sims
# suffixes_ = ['median_no_errors-thermo_props-CC_classifications']*num_sims
# suffixes_ = ['median_no_errors-all_props-no_log-CC_classifications']*num_sims
# suffixes_ = ['median_no_errors-all_props-CC_classifications-v2']*num_sims
suffixes_ = ['median_no_errors-all_props']*num_sims
suffixes_ = ['median_no_errors-all_props-with_median']*num_sims
# suffixes_ = ['median_no_errors-all_props-full_annulus_res']*num_sims
# suffixes_ = ['new_xscale_Physical-median_no_errors-all_props']*num_sims
# suffixes_ = ['new_xscale_Physical-median_no_errors-all_props-full_annulus_res']*num_sims
# suffixes_ = ['median_no_errors-all_props-relaxed_res_cuts']*num_sims
# suffixes_ = ['median_no_errors-all_props-full_annulus_res-relaxed_res_cuts']*num_sims
suffixes_ = ['median_no_errors-all_props-quartiles']*num_sims

codes_ = ['Simba', 'Simba-C'] + ['Simba']*3 + ['Simba-C']
# codes_ = ['Simba-C']
codes_ = ['Simba', 'Simba-C']

sizes_ = ['N1024L100']*num_sims
# sizes_ = ['N512L50']*num_sims

change_xscales_ = [False]*num_sims
# change_xscales_ = [True]*num_sims

new_xscales_ = ['None']*num_sims
# new_xscales_ = ['Physical']*num_sims

avg_types_ = ['median']*num_sims

use_errorses_ = [False]*num_sims



res_ndims_ = [3]*num_sims
# res_ndims_ = [2]*num_sims

res_filters_ = ['Sphere']*num_sims
# res_filters_ = ['Cylinder']*num_sims

profiles_are_2d_but_res_file_is_3ds_ = [True]*num_sims
# profiles_are_2d_but_res_file_is_3ds_ = [False]*num_sims

use_3d_for_2d_reses_ = [True]*num_sims
# use_3d_for_2d_reses_ = [False]*num_sims



include_resolution_cuts_ = [True]*num_sims

# min_resolved_counts_per_halos_ = [16]*num_sims
min_resolved_counts_per_halos_ = [64]*num_sims

min_resolved_counts_per_halo_cores_ = [16]*num_sims

# min_resolved_counts_per_bins_ = [16]*num_sims
min_resolved_counts_per_bins_ = [64]*num_sims

min_total_haloses_ = [0]*num_sims
# min_total_haloses_ = [1]*num_sims

min_resolved_haloses_ = [0]*num_sims
# min_resolved_haloses_ = [1]*num_sims
min_resolved_haloses_ = [5]*num_sims

# min_frac_resolved_haloses_ = [0.0]*num_sims
min_frac_resolved_haloses_ = [0.16]*num_sims

# calc_thermo_propses_ = [False]*num_sims
calc_thermo_propses_ = [True]*num_sims

# calc_metal_propses_ = [False]*num_sims
calc_metal_propses_ = [True]*num_sims

# calc_log_propses_ = [False]*num_sims
calc_log_propses_ = [True]*num_sims

classify_haloses_ = [False]*num_sims
# classify_haloses_ = [True]*num_sims



names_ = [sim+'_'+size+'_igrm_'+str(ndim)+'d_'+filt+'_profiles_Nlgal_'+str(nlgal)+'_windsremoved_xscale_'+xscale+'_intermediatebins_temp_cut='+temp_cut+'_nh_cut='+nh_cut for sim,size,ndim,filt,nlgal,xscale,temp_cut,nh_cut in zip(sims_,sizes_,ndims_,filters_,n_lgals_,xscales_,temp_cuts_,nh_cuts_)]

res_names_= [sim+'_'+size+'_igrm_'+str(ndim)+'d_'+filt+'_profiles_Nlgal_'+str(nlgal)+'_windsremoved_xscale_'+xscale+'_intermediatebins_temp_cut='+temp_cut+'_nh_cut='+nh_cut for sim,size,ndim,filt,nlgal,xscale,temp_cut,nh_cut in zip(sims_,sizes_,res_ndims_,res_filters_,n_lgals_,xscales_,temp_cuts_,nh_cuts_)]

three_d_names_ = [sim+'_'+size+'_igrm_3d_Sphere_profiles_Nlgal_'+str(nlgal)+'_windsremoved_xscale_'+xscale+'_intermediatebins_temp_cut='+temp_cut+'_nh_cut='+nh_cut for sim,size,nlgal,xscale,temp_cut,nh_cut in zip(sims_,sizes_,n_lgals_,xscales_,temp_cuts_,nh_cuts_)]

two_d_names_ = [sim+'_'+size+'_igrm_2d_Cylinder_profiles_Nlgal_'+str(nlgal)+'_windsremoved_xscale_'+xscale+'_intermediatebins_temp_cut='+temp_cut+'_nh_cut='+nh_cut for sim,size,nlgal,xscale,temp_cut,nh_cut in zip(sims_,sizes_,n_lgals_,xscales_,temp_cuts_,nh_cuts_)]






snapshot_numbers = ['151', '125', '105', '078']
snapshot_numbers = ['151']

redshift_strings = ['0.000', '0.490', '0.993', '2.025']
redshift_strings = ['0.000']







dirs = []
wall_times = []
snapshots = []
names = []
suffixes = []
codes = []
sizes = []
ndims = []
change_xscales = []    
new_xscales = []
avg_types = []
use_errorses = []
res_files = []
profiles_are_2d_but_res_file_is_3ds = []
use_3d_for_2d_reses = []
three_d_files = []
two_d_files = []    
include_resolution_cuts = []
min_resolved_counts_per_halos = []
min_resolved_counts_per_halo_cores = []
min_resolved_counts_per_bins = []
min_total_haloses = []
min_resolved_haloses = []
min_frac_resolved_haloses = []
calc_thermo_propses = []
calc_metal_propses = []
calc_log_propses = []
classify_haloses = []

for jj in range(num_sims):
    for snap_num, z_str in zip (snapshot_numbers, redshift_strings):
        dirs.append(dirs_[jj])
        wall_times.append(wall_times_[jj])
        snapshots.append('snapshot_'+snap_num+'_z'+z_str)
        names.append(names_[jj])
        suffixes.append(suffixes_[jj])
        codes.append(codes_[jj])
        sizes.append(sizes_[jj])
        ndims.append(ndims_[jj])
        change_xscales.append(change_xscales_[jj])
        new_xscales.append(new_xscales_[jj])
        avg_types.append(avg_types_[jj])
        use_errorses.append(use_errorses_[jj])
        res_files.append(res_names_[jj])
        profiles_are_2d_but_res_file_is_3ds.append(profiles_are_2d_but_res_file_is_3ds_[jj])
        use_3d_for_2d_reses.append(use_3d_for_2d_reses_[jj])
        three_d_files.append(three_d_names_[jj])
        two_d_files.append(two_d_names_[jj])
        include_resolution_cuts.append(include_resolution_cuts_[jj])
        min_resolved_counts_per_halos.append(min_resolved_counts_per_halos_[jj])
        min_resolved_counts_per_halo_cores.append(min_resolved_counts_per_halo_cores_[jj])
        min_resolved_counts_per_bins.append(min_resolved_counts_per_bins_[jj])
        min_total_haloses.append(min_total_haloses_[jj])
        min_resolved_haloses.append(min_resolved_haloses_[jj])
        min_frac_resolved_haloses.append(min_frac_resolved_haloses_[jj])
        calc_thermo_propses.append(calc_thermo_propses_[jj])
        calc_metal_propses.append(calc_metal_propses_[jj])
        calc_log_propses.append(calc_log_propses_[jj])
        classify_haloses.append(classify_haloses_[jj])
        
        
        
num_jobs = len(codes)











# # dirs = ['/scratch/b/babul/aspadawe/analysis/simba/simba/final']*8 + ['/scratch/b/babul/aspadawe/analysis/simba/simba-c/final']*8
# # dirs = ['/scratch/b/babul/aspadawe/analysis/simba/simba/final']*4 + ['/scratch/b/babul/aspadawe/analysis/simba/simba-c/final']*4
# dirs = ['/scratch/b/babul/aspadawe/analysis/simba/simba/final']*2 + ['/scratch/b/babul/aspadawe/analysis/simba/simba-c/final']*2

# num_jobs = len(dirs)
# # num_jobs = 1

# wall_time = ['10:0:0']*num_jobs

# snapshots = ['snapshot_151']*num_jobs

# # names = ['Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*4 + ['Simba_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*4 + ['Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*4 + ['Simba-C_Kobayashi_Dust_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*4
# # names = ['Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*2 + ['Simba_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*2 + ['Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*2 + ['Simba-C_Kobayashi_Dust_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*2
# # names = ['Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3', 'Simba_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3', 'Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3', 'Simba-C_Kobayashi_Dust_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']
# names = ['Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3', 'Simba_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3', 'Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3', 'Simba-C_Kobayashi_Dust_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']

# # suffixes = [
# #     'mean_with_errors-all_props_no_classifications-v1',
# #     'mean_with_errors-all_props_no_classifications-v2',
# #     'median_with_errors-all_props_no_classifications-v1',
# #     'median_with_errors-all_props_no_classifications-v2',
# #     'mean_with_errors-metal_props_no_log_no_classifications-v1',
# #     'mean_with_errors-metal_props_no_log_no_classifications-v2',
# #     'median_with_errors-metal_props_no_log_no_classifications-v1',
# #     'median_with_errors-metal_props_no_log_no_classifications-v2',
# #     'mean_with_errors-all_props_no_classifications-v1',
# #     'mean_with_errors-all_props_no_classifications-v2',
# #     'median_with_errors-all_props_no_classifications-v1',
# #     'median_with_errors-all_props_no_classifications-v2',
# #     'mean_with_errors-metal_props_no_log_no_classifications-v1',
# #     'mean_with_errors-metal_props_no_log_no_classifications-v2',
# #     'median_with_errors-metal_props_no_log_no_classifications-v1',
# #     'median_with_errors-metal_props_no_log_no_classifications-v2',]
# # suffixes = [
# #     'mean_no_errors-all_props_no_classifications-v1',
# #     'mean_no_errors-all_props_no_classifications-v2',
# #     'median_no_errors-all_props_no_classifications-v1',
# #     'median_no_errors-all_props_no_classifications-v2',
# #     'mean_no_errors-metal_props_no_log_no_classifications-v1',
# #     'mean_no_errors-metal_props_no_log_no_classifications-v2',
# #     'median_no_errors-metal_props_no_log_no_classifications-v1',
# #     'median_no_errors-metal_props_no_log_no_classifications-v2',
# #     'mean_no_errors-all_props_no_classifications-v1',
# #     'mean_no_errors-all_props_no_classifications-v2',
# #     'median_no_errors-all_props_no_classifications-v1',
# #     'median_no_errors-all_props_no_classifications-v2',
# #     'mean_no_errors-metal_props_no_log_no_classifications-v1',
# #     'mean_no_errors-metal_props_no_log_no_classifications-v2',
# #     'median_no_errors-metal_props_no_log_no_classifications-v1',
# #     'median_no_errors-metal_props_no_log_no_classifications-v2',]
# # suffixes = [
# #     'mean_with_errors-all_props_no_classifications-v0',
# #     'median_with_errors-all_props_no_classifications-v0',
# #     'mean_with_errors-metal_props_no_log_no_classifications-v0',
# #     'median_with_errors-metal_props_no_log_no_classifications-v0',
# #     'mean_with_errors-all_props_no_classifications-v0',
# #     'median_with_errors-all_props_no_classifications-v0',
# #     'mean_with_errors-metal_props_no_log_no_classifications-v0',
# #     'median_with_errors-metal_props_no_log_no_classifications-v0',]
# # suffixes = [
# #     'mean_no_errors-all_props_no_classifications-v0',
# #     'median_no_errors-all_props_no_classifications-v0',
# #     'mean_no_errors-metal_props_no_log_no_classifications-v0',
# #     'median_no_errors-metal_props_no_log_no_classifications-v0',
# #     'mean_no_errors-all_props_no_classifications-v0',
# #     'median_no_errors-all_props_no_classifications-v0',
# #     'mean_no_errors-metal_props_no_log_no_classifications-v0',
# #     'median_no_errors-metal_props_no_log_no_classifications-v0',]
# # suffixes = [
# #     'mean_no_errors-all_props_no_classifications-v3',
# #     'median_no_errors-all_props_no_classifications-v3',
# #     'mean_no_errors-metal_props_no_log_no_classifications-v3',
# #     'median_no_errors-metal_props_no_log_no_classifications-v3',
# #     'mean_no_errors-all_props_no_classifications-v3',
# #     'median_no_errors-all_props_no_classifications-v3',
# #     'mean_no_errors-metal_props_no_log_no_classifications-v3',
# #     'median_no_errors-metal_props_no_log_no_classifications-v3',]
# # suffixes = [
# #     'median_no_errors-all_props_no_classifications-v3',
# #     'median_no_errors-metal_props_no_log_no_classifications-v3',
# #     'median_no_errors-all_props_no_classifications-v3',
# #     'median_no_errors-metal_props_no_log_no_classifications-v3',
# # ]
# # suffixes = [
# #     'median_no_errors-all_props_CC_classifications-v3',
# #     'median_no_errors-metal_props_no_log_CC_classifications-v3',
# #     'median_no_errors-all_props_CC_classifications-v3',
# #     'median_no_errors-metal_props_no_log_CC_classifications-v3',
# # ]
# suffixes = [
#     'median_no_errors-all_props_CC_classifications-v3',
#     'median_no_errors-all_props_no_log_CC_classifications-v3',
#     'median_no_errors-all_props_CC_classifications-v3',
#     'median_no_errors-all_props_no_log_CC_classifications-v3',
# ]

# # codes = ['Simba']*8 + ['Simba-C']*8
# # codes = ['Simba']*4 + ['Simba-C']*4
# codes = ['Simba']*2 + ['Simba-C']*2
    
# sizes = ['N1024L100']*num_jobs

# # ndims = [3,3,3,3,2,2,2,2,3,3,3,3,2,2,2,2]
# # ndims = [3,3,2,2,3,3,2,2]
# ndims = [3,2,3,2]
    
# change_xscales = [False]*num_jobs
    
# new_xscales = ['None']*num_jobs

# # avg_types = ['mean','mean','median','median','mean','mean','median','median',
# #              'mean','mean','median','median','mean','mean','median','median',]
# # avg_types = ['mean','median','mean','median',
# #              'mean','median','mean','median',]
# avg_types = ['median']*num_jobs
    
# use_errorses=[False]*num_jobs

# # res_files = ['Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*8 + ['Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*8
# # res_files = ['Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*4 + ['Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*4
# res_files = ['Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*2 + ['Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*2

# # profiles_are_2d_but_res_file_is_3ds = [False]*4 + [True]*4 + [False]*4 + [True]*4
# # profiles_are_2d_but_res_file_is_3ds = [False]*2 + [True]*2 + [False]*2 + [True]*2
# profiles_are_2d_but_res_file_is_3ds = [False, True, False, True]#*2

# use_3d_for_2d_reses = [True]*num_jobs

# three_d_files = ['Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*2 + ['Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*2

# two_d_files = ['Simba_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*2 + ['Simba-C_Kobayashi_Dust_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3']*2
    
# include_resolution_cuts = [True]*num_jobs

# # min_resolved_counts_per_halos = [64,0,64,0,64,0,64,0,
# #                                  64,0,64,0,64,0,64,0,]
# # min_resolved_counts_per_halos = [0,0,0,0,0,0,0,0,]
# min_resolved_counts_per_halos = [64]*num_jobs
    
# min_resolved_counts_per_bins = [64]*num_jobs

# # min_total_haloses = [0,5,0,5,0,5,0,5,
# #                      0,5,0,5,0,5,0,5,]
# # min_total_haloses = [0,0,0,0,0,0,0,0,]
# min_total_haloses = [0]*num_jobs

# # min_resolved_haloses = [5,0,5,0,5,0,5,0,
# #                         5,0,5,0,5,0,5,0,]
# # min_resolved_haloses = [0,0,0,0,0,0,0,0,]
# min_resolved_haloses = [5]*num_jobs

# # min_frac_resolved_haloses = [0,0.16,0,0.16,0,0.16,0,0.16,
# #                              0,0.16,0,0.16,0,0.16,0,0.16,]
# # min_frac_resolved_haloses = [0,0,0,0,0,0,0,0,]
# min_frac_resolved_haloses = [0.16]*num_jobs

# # calc_thermo_propses = [True,True,True,True,False,False,False,False,
# #                        True,True,True,True,False,False,False,False,]
# # calc_thermo_propses = [True,True,False,False,
# #                        True,True,False,False,]
# # calc_thermo_propses = [True,False,True,False]
# calc_thermo_propses = [True]*num_jobs
    
# calc_metal_propses = [True]*num_jobs

# # calc_log_propses = [True,True,True,True,False,False,False,False,
# #                     True,True,True,True,False,False,False,False,]
# # calc_log_propses = [True,True,False,False,
# #                     True,True,False,False,]
# calc_log_propses = [True,False,True,False,]
    
# classify_haloses = [True]*num_jobs



with open('job_calc_profile_properties_and_statistics_template.sh', 'r') as f:
    job_data = f.read()

job_dir = './job_scripts_for_profiles'
if not os.path.isdir(job_dir):
    os.makedirs(job_dir, mode = 0o755, exist_ok = True)


## Create job file, substitute in parameters, submit, (delete job file?)
new_jobs = []
for ii in range(num_jobs):
    job_file = 'job_calc_profiles_properties_and_statistics_%s.sh' % ii
    param_dict = {
#         'JOB_NAME':job_names[ii],
        'WALL_TIME':wall_times[ii],
        'SNAPSHOT':snapshots[ii],
        'DIR':dirs[ii],
        'INPUT_NAME':names[ii],
#         'SUFFIX':'physical_xscale-' + suffixes[ii],
#         'SUFFIX':suffixes[ii],
        'SUFFIX':suffixes[ii],# + '-min_part_%s' % min_resolved_counts_per_halos[ii] + '-min_res_halos_%s' % min_resolved_counts_per_bins[ii],
#         'SUFFIX':'physical_xscale-' + suffixes[ii] + '-min_part_%s' % min_resolved_counts_per_halos[ii],
        'CODE':codes[ii],
        'SIZE':sizes[ii],
        'NDIM':ndims[ii],
        'CHANGE_XSCALE':change_xscales[ii],
        'NEW_XSCALE':new_xscales[ii],
        'AVG_TYPE':avg_types[ii],
        'USE_ERRORS':use_errorses[ii],
        'RES_FILE':res_files[ii],
        'PROFILES_ARE_2D_BUT_RES_PROFILE_IS_3D':profiles_are_2d_but_res_file_is_3ds[ii],
        'USE_3D_FOR_2D_RES':use_3d_for_2d_reses[ii],
        'THREE_D_FILE':three_d_files[ii],
        'TWO_D_FILE':two_d_files[ii],
        'INCLUDE_RESOLUTION_CUTS':include_resolution_cuts[ii],
        'MIN_RESOLVED_COUNTS_PER_HALO_BIN':min_resolved_counts_per_halos[ii],
        'MIN_RESOLVED_COUNTS_PER_HALO_CORE':min_resolved_counts_per_halo_cores[ii],
        'MIN_RESOLVED_COUNTS_PER_BIN':min_resolved_counts_per_bins[ii],
        'MIN_TOTAL_HALOS':min_total_haloses[ii],
        'MIN_RESOLVED_HALOS':min_resolved_haloses[ii],
        'MIN_FRAC_RESOLVED_HALOS':min_frac_resolved_haloses[ii],
        'CALC_THERMO_PROPS':calc_thermo_propses[ii],
        'CALC_METAL_PROPS':calc_metal_propses[ii],
        'CLASSIFY_HALOS':classify_haloses[ii],
        'CALC_LOG_PROPS':calc_log_propses[ii],
    }
    
    ## Create job_name
    job_name = '%s-%s-%s-%sd-%s-profiles' % (param_dict['SNAPSHOT'], param_dict['CODE'], param_dict['SIZE'], 
                                          param_dict['NDIM'], param_dict['SUFFIX'])
#     job_name = param_dict['INPUT_NAME']
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