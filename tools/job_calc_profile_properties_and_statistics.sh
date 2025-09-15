#!/bin/bash -l
#########################################################
##job_name=calcs_Simba-C_Dust_512/50_Nlgal_3_profiles_5e5K

##SBATCH -J calcs_Simba_512/50_5e5K_Nlgal3_2d_profiles
##SBATCH -J calcs_Simba-C_Dust_512/50_5e5K_Nlgal3_2d_profiles
##SBATCH -J calcs_Simba-C_NoDust_512/50_Nlgal_3_profiles_5e5K
##SBATCH -J calcs_Simba_1024/100_5e5K_Nlgal3_2d_median_profiles
#SBATCH -J calcs_Simba-C_Dust_1024/100_5e5K_Nlgal3_2d_median_profiles
##SBATCH -J calcs_Simba_NoAGN_512/50_5e5K_Nlgal3_2d_profiles
##SBATCH -J calcs_Simba_NoFB_512/50_5e5K_Nlgal3_2d_profiles
##SBATCH -J calcs_Simba_NoJet_512/50_Nlgal_3_profiles_5e5K

##SBATCH -J $job_name

#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-babul-ad
##SBATCH --account=aspadawe
##SBATCH -o /scratch/b/babul/aspadawe/analysis/simba/slurm_files/slurm-%j-calcs_Simba_1024-100_5e5K_Nlgal3_2d_median_v2_R500_profiles.out
##SBATCH -o /scratch/b/babul/aspadawe/analysis/simba/slurm_files/slurm-%j-calcs_Simba_512-50_5e5K_Nlgal3_2d_physical_profiles.out
#SBATCH -o /scratch/b/babul/aspadawe/analysis/simba/slurm_files/slurm-%j-calcs_Simba-C_Dust_1024-100_5e5K_Nlgal3_2d_median_v2_profiles.out
##SBATCH -o /scratch/b/babul/aspadawe/analysis/simba/slurm_files/slurm-%j-calcs_Simba-C_Dust_512-50_5e5K_Nlgal3_2d_physical_profiles.out
##SBATCH -o /scratch/b/babul/aspadawe/analysis/simba/slurm_files/slurm-%j-calcs_Simba_NoAGN_512-50_5e5K_Nlgal3_2d_profiles.out
#########################################################
##SBATCH --time=24:00:0
#SBATCH --time=2:59:0
##SBATCH --time=0:40:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=80
#SBATCH --mem=0
#########################################################

# conda activate XIGrM
conda activate calc_profiles


## Simba-C N256L25
# python calc_profile_properties_and_statistics_good.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c/new' --name='Simba-C_Kobayashi_Dust_N256L25_igrm_2d_cylindrical_profiles_2_Nlgal_3_windsremoved_r500_intermediatebins_testing__temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='all' --code='Simba-C' --size='N256L25' --ndim=2

## Simba N512L50
# python calc_profile_properties_and_statistics_good.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba/final' --name='Simba_N512L50_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_physical_units_intermediatebins_og_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='all_no_classification' --code='Simba' --size='N512L50' --ndim=2 #good
# python calc_profile_properties_and_statistics_new.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba/final' --name='Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='weighted_median_with_errors_all_no_classification' --code='Simba' --size='N1024L100' --ndim=3 --change_xscale=False --new_xscale=None --use_ref_bins=True --avg_type='median' --use_errors=True --min_resolved_counts_per_halo=64 --min_resolved_counts_per_bin=64 --min_total_halos=0 --min_resolved_halos=5 --min_frac_resolved_halos=0.16 #good

# python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_test_rbin_weighting_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba' --size='N512L50'

## Simba-C Dust N512L50
# python calc_profile_properties_and_statistics_good.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c/final' --name='Simba-C_Kobayashi_Dust_N512L50_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_physical_units_intermediatebins_og_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='all_no_classification' --code='Simba-C' --size='N512L50' --ndim=2 #good
# python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c' --name='Simba-C_Kobayashi_Dust_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_backup_halo_props_test_rbin_weighting_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba-C' --size='N512L50'

## Simba-C No Dust N512L50
#python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c/new' --name='Simba-C_Kobayashi_NoDust_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_from_reniers_project_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='all_noCC_classification' --code='Simba-C' --size='N512L50'

## Simba N1024L100
# python calc_profile_properties_and_statistics_good.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba/final' --name='Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='all_no_classification' --code='Simba' --size='N1024L100' --ndim=3  ## good
python calc_profile_properties_and_statistics_ref_bins.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba/final' --name='Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='mean_no_errors-all_props_no_classifications-v0' --code='Simba' --size='N1024L100' --ndim=3 --change_xscale=False --new_xscale=None --avg_type='mean' --use_errors=False --min_resolved_counts_per_halo=0 --min_resolved_counts_per_bin=0 --min_total_halos=0 --min_resolved_halos=0 --min_frac_resolved_halos=0 --calc_thermo_props=True --calc_metal_props=True --classify_halos=False --calc_log_props=True #real grood

#Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut\=5e5\ K_nh_cut\=0.13\ cm\*\*-3.npy 

## Simba-C Dust N1024L100
# python calc_profile_properties_and_statistics_good.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c/final' --name='Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='all_no_classification' --code='Simba-C' --size='N1024L100' --ndim=3  ## good
# python calc_profile_properties_and_statistics_ref_bins.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c/final' --name='Simba-C_Kobayashi_Dust_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='median-metal_props_no_log_no_classifications-v2' --code='Simba-C' --size='N1024L100' --ndim=2 --change_xscale=False --new_xscale=None --avg_type='median' --use_errors=True --min_resolved_counts_per_halo=0 --min_resolved_counts_per_bin=64 --min_total_halos=5 --min_resolved_halos=0 --min_frac_resolved_halos=0.16 --calc_thermo_props=False --calc_metal_props=True --classify_halos=False --calc_log_props=False #real good

# python calc_profile_properties_and_statistics_good.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c/new' --name='Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='metallicity_no_log_no_classification' --code='Simba-C' --size='N1024L100'
#python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c' --name='Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba-C' --size='N1024L100'


## Simba No AGN N512L50
# python calc_profile_properties_and_statistics_good.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba/final' --name='Simba_NoAGN_N512L50_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_from_reniers_scratch_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='metallicity_no_log_no_classification' --code='Simba' --size='N512L50' --ndim=2  ## good
#python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_NoAGN_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_from_reniers_scratch_test_rbin_weighting_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba' --size='N512L50'
#python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_NoAGN_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_from_reniers_project_test_rbin_weighting_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba' --size='N512L50'

## Simba No FB N512L50
# python calc_profile_properties_and_statistics_good.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba/final' --name='Simba_NoFB_N512L50_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good_from_reniers_scratch_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --suffix='metallicity_no_log_no_classification' --code='Simba' --size='N512L50' --ndim=2  ## good
#python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_NoFB_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_from_reniers_scratch_test_rbin_weighting_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba' --size='N512L50'
#python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_NoFB_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_from_reniers_project_test_rbin_weighting_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba' --size='N512L50'

## Simba No Jet N512L50
#python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_NoJet_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_from_reniers_project_test_rbin_weighting_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba' --size='N512L50'


########################################################












#python calc_profile_properties_and_statistics.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c' --name='Simba-C_Kobayashi_Dust_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_includeall' --code='Simba-C'python calc_profile_properties_and_statistics.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c' --name='Simba-C_Kobayashi_Dust_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba-C' --size='N512L50'
#python calc_profile_properties_and_statistics.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c' --name='Simba-C_Kobayashi_Dust_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba-C' --size='N512L50'
#python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c' --name='Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba-C' --size='N1024L100'
#python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba-c' --name='Simba-C_Kobayashi_Dust_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba-C' --size='N512L50'

#python calc_profile_properties_and_statistics.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins' --code='Simba' --size='N1024L1000'
#python calc_profile_properties_and_statistics_reduced.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba' --size='N1024L100'

#python calc_profile_properties_and_statistics.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_includeall' --code='Simba' --size='N512L50'
#python calc_profile_properties_and_statistics.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_temp_cut=5e5 K_nh_cut=0.13 cm**-3' --code='Simba' --size='N512L50'

#python calc_profile_properties_and_statistics.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_NoFB_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_includeall' --code='Simba'
#python calc_profile_properties_and_statistics.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_NoAGN_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_includeall' --code='Simba'
#python calc_profile_properties_and_statistics.py --snap='snapshot_151' --dir='/scratch/b/babul/aspadawe/analysis/simba/simba' --name='Simba_NoJet_N512L50_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_includeall' --code='Simba'

########################################################
