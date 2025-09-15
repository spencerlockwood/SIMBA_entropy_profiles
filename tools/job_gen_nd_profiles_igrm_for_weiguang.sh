#!/bin/bash -l
#########################################################
#SBATCH -J Simba_Simba-C_N1024L100_2d_3d_profiles_Nlgal_3

#SBATCH --mail-user=apadawer@uvic.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-babul-ad
#SBATCH -o /scratch/b/babul/aspadawe/analysis/simba/slurm_files/slurm-%j.out
#########################################################
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
#########################################################

conda activate XIGrM


## Simba N1024L100
# python gen_nd_profiles_igrm.py --code='simba' --snap='/scratch/b/babul/renierht/output_directory/Simba_N1024L100/snapshot_151.hdf5' --ahf_basename='/scratch/b/babul/wcui/cosmo_simba_ahf/AHF_MVIR_151.z0.000.AHF_' --nthreads=40 --prop='/scratch/b/babul/renierht/output_directory/Simba_N1024L100/halo_results/halo_z0_NW.hdf5' --save='/scratch/b/babul/aspadawe/analysis/simba/simba/new/snapshot_151_Simba_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good' --filter='Sphere' --profile_type='Profile' --ndim=3 --weight_by='mass' --xscale='R500' --temp_cut='5e5 K' --nh_cut='0.13 cm**-3' --N_lgal=3 --log_M500_cut=13.0

# python gen_nd_profiles_igrm.py --code='simba' --snap='/scratch/b/babul/renierht/output_directory/Simba_N1024L100/snapshot_151.hdf5' --ahf_basename='/scratch/b/babul/wcui/cosmo_simba_ahf/AHF_MVIR_151.z0.000.AHF_' --nthreads=40 --prop='/scratch/b/babul/renierht/output_directory/Simba_N1024L100/halo_results/halo_z0_NW.hdf5' --save='/scratch/b/babul/aspadawe/analysis/simba/simba/new/snapshot_151_Simba_N1024L100_igrm_2d_spherical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good' --filter='Sphere' --profile_type='Profile' --ndim=2 --weight_by='mass' --xscale='R500' --temp_cut='5e5 K' --nh_cut='0.13 cm**-3' --N_lgal=3 --log_M500_cut=13.0

python gen_nd_profiles_igrm.py --code='simba' --snap='/scratch/b/babul/renierht/output_directory/Simba_N1024L100/snapshot_151.hdf5' --ahf_basename='/scratch/b/babul/wcui/cosmo_simba_ahf/AHF_MVIR_151.z0.000.AHF_' --nthreads=40 --prop='/scratch/b/babul/renierht/output_directory/Simba_N1024L100/halo_results/halo_z0_NW.hdf5' --save='/scratch/b/babul/aspadawe/analysis/simba/simba/new/snapshot_151_Simba_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good' --filter='Cylinder' --profile_type='Profile' --ndim=2 --weight_by='mass' --xscale='R500' --temp_cut='5e5 K' --nh_cut='0.13 cm**-3' --N_lgal=3 --log_M500_cut=13.0


## Simba-C Dust N1024L100
# python gen_nd_profiles_igrm.py --code='simba-c' --snap='/scratch/b/babul/renierht/output_directory/Kobayashi_N1024L100/snapshot_151.hdf5' --ahf_basename='/scratch/b/babul/wcui/cosmo_simbac_ahf/AHF_MVIR_151.z0.000.AHF_' --nthreads=40 --prop='/scratch/b/babul/renierht/output_directory/Kobayashi_N1024L100/halo_results/halo_z0_NW.hdf5' --save='/scratch/b/babul/aspadawe/analysis/simba/simba-c/new/snapshot_151_Simba-C_Kobayashi_Dust_N1024L100_igrm_3d_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good' --filter='Sphere' --profile_type='Profile' --ndim=3 --weight_by='mass' --xscale='R500' --temp_cut='5e5 K' --nh_cut='0.13 cm**-3' --N_lgal=3 --log_M500_cut=13.0

# python gen_nd_profiles_igrm.py --code='simba-c' --snap='/scratch/b/babul/renierht/output_directory/Kobayashi_N1024L100/snapshot_151.hdf5' --ahf_basename='/scratch/b/babul/wcui/cosmo_simbac_ahf/AHF_MVIR_151.z0.000.AHF_' --nthreads=40 --prop='/scratch/b/babul/renierht/output_directory/Kobayashi_N1024L100/halo_results/halo_z0_NW.hdf5' --save='/scratch/b/babul/aspadawe/analysis/simba/simba-c/new/snapshot_151_Simba-C_Kobayashi_Dust_N1024L100_igrm_2d_spherical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good' --filter='Sphere' --profile_type='Profile' --ndim=2 --weight_by='mass' --xscale='R500' --temp_cut='5e5 K' --nh_cut='0.13 cm**-3' --N_lgal=3 --log_M500_cut=13.0

python gen_nd_profiles_igrm.py --code='simba-c' --snap='/scratch/b/babul/renierht/output_directory/Kobayashi_N1024L100/snapshot_151.hdf5' --ahf_basename='/scratch/b/babul/wcui/cosmo_simbac_ahf/AHF_MVIR_151.z0.000.AHF_' --nthreads=40 --prop='/scratch/b/babul/renierht/output_directory/Kobayashi_N1024L100/halo_results/halo_z0_NW.hdf5' --save='/scratch/b/babul/aspadawe/analysis/simba/simba-c/new/snapshot_151_Simba-C_Kobayashi_Dust_N1024L100_igrm_2d_cylindrical_profiles_Nlgal_3_windsremoved_r500_intermediatebins_good' --filter='Cylinder' --profile_type='Profile' --ndim=2 --weight_by='mass' --xscale='R500' --temp_cut='5e5 K' --nh_cut='0.13 cm**-3' --N_lgal=3 --log_M500_cut=13.0


########################################################
