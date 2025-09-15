#!/bin/bash -l
#########################################################
#SBATCH -J calcs_JOB_NAME

#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-babul-ad
#SBATCH -o /scratch/b/babul/aspadawe/analysis/simba/slurm_files/slurm-%j-calcs_JOB_NAME.out
#########################################################
#SBATCH --time=WALL_TIME
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#########################################################

conda activate calc_profiles

python ../calc_profile_properties_and_statistics_ref_bins.py --snap='SNAPSHOT' --dir='DIR' --name='INPUT_NAME' --suffix='SUFFIX' --code='CODE' --size='SIZE' --ndim=NDIM --change_xscale=CHANGE_XSCALE --new_xscale='NEW_XSCALE' --avg_type='AVG_TYPE' --use_errors=USE_ERRORS --res_file='RES_FILE' --profiles_are_2d_but_res_file_is_3d=PROFILES_ARE_2D_BUT_RES_PROFILE_IS_3D --use_3d_for_2d_res=USE_3D_FOR_2D_RES --three_d_file='THREE_D_FILE' --two_d_file='TWO_D_FILE' --include_resolution_cuts=INCLUDE_RESOLUTION_CUTS --min_resolved_counts_per_halo=MIN_RESOLVED_COUNTS_PER_HALO_BIN --min_resolved_counts_per_halo_core=MIN_RESOLVED_COUNTS_PER_HALO_CORE --min_resolved_counts_per_bin=MIN_RESOLVED_COUNTS_PER_BIN --min_total_halos=MIN_TOTAL_HALOS --min_resolved_halos=MIN_RESOLVED_HALOS --min_frac_resolved_halos=MIN_FRAC_RESOLVED_HALOS --calc_thermo_props=CALC_THERMO_PROPS --calc_metal_props=CALC_METAL_PROPS --classify_halos=CLASSIFY_HALOS --calc_log_props=CALC_LOG_PROPS

########################################################
