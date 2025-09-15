#!/bin/bash -l
#########################################################
#SBATCH -J gen_JOB_NAME

#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-babul-ad
#SBATCH -o /scratch/b/babul/aspadawe/analysis/simba/slurm_files/slurm-%j-gen_JOB_NAME.out
#########################################################
#SBATCH --time=WALL_TIME
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#########################################################

conda activate XIGrM

python ../gen_nd_profiles_igrm.py --code='CODE' --redshift=REDSHIFT --snap='SNAPSHOT' --ahf_basename='AHF_BASENAME' --nthreads=NTHREADS --prop='PROP_FILE' --save='SAVE_FILE' --filter='FILTER' --profile_type='PROFILE_TYPE' --ndim=NDIM --weight_by='WEIGHT_BY' --xscale='XSCALE' --temp_cut='TEMP_CUT' --nh_cut='NH_CUT' --N_lgal=N_LGAL --log_M500_cut=LOG_M500_CUT

########################################################
