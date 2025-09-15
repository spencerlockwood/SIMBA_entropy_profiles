'''
SLOW mpi/ array job? - DO NOT try access the same hdf5 file at the same time across ranks/array jobs!!! also run with srun not mpirun when doing interactive job



turb_venv
export PYTHONPATH=$PYTHONPATH:/home/b/babul/fjenning/scratch/xray_bubbles/analysis_scripts/
export PYTHONPATH=$PYTHONPATH:/home/b/babul/fjenning/scratch/xray_bubbles/PAPER2_scripts/
module load intel/2019u3  gsl  hdf5 intelmpi


GO_STRAIGHT_TO_PLOTTING_PRESAVED_ANALYSIS requires the analysis to have been performed once already
Will have filenotfound errors if you are trying to find data before you ran the analysis - e.g. your PADDING before the snapnums you want to analyse isnt big enough (e.g. ANALYSE_BUBBLE_PARTS_FROM_PREVIOUS_N_MYR)

'''

'''
"recent" - e.g. bubble particles or condensed particles which became bubbles/cold some time ago and we still keep track of them at current snapnum
"instantaneous" - e.g. bubble particles or condensed particles which became bubbles/cold between now and the previous snapshot



- check T vs time plots for particles to check they are genuinely condensed and not just dipping above
- Fix DUPLICATE CONDENSED PIDs PRESENT
- what should ANALYSE_SATELLITE_WAKE_FROM_PREVIOUS_N_MYR be?

'''




'''
What else needs to go in grackle_pars ? UV background?
- Note, have added a radius cut for condensed particles. Currently is 10kpc
- Should we use median or mean (current) radial profile when looking for perturbations?
- Particles tagged as bubbles even within min_dist_in_kpc must be recoupled wind particles, which we should include
- But we should remove stellar winds using the ID list from the bhwinds file - DONE, as best as an. Still will include heating by non-prevouslywind particles
- What should ANALYSE_SATELLITE_WAKE_FROM_PREVIOUS_N_MYR be?
- Must check that my method for determining particles in satellite wake vs in front of satellite is valid. Have the "in-front" particles from snapnum_this_rank - 1 ALWAYS moved behind the satellite before we reach snapnum_this_rank?
- Check that after my fix (sorting the set in fin_condensing_times_of_cold_gas.py before splitting to MPI rank) means that condensed_gas_unique_PIDs is now always the same when I re-reun fin_condensing_times_of_cold_gas.py
- Currently I am excluding satellites which have no particles identified as "in front of them" at the present snapshot. These are likely to be spurious detections by the halo-finder. Also exluding those with fewer than 3 particles
- Should check that bubbleparts are still "hot"
- Should only plot particles within a slice
- Think about difference in sphere radii between analysis scripts
- ATM central galaxy is being included in the satellite wake calcs
- check that 2dhist is SUMMING in each bin
- KDE for 2dhist?
- Results seem to get fucked up for the analysis if this script stops then is restarted to fill in the gaps... why??
- For 'recently' condensed gas for the radius vs time vs condensed_mass histplot, need to think carefully if we should include gas that may have heated up between condensing and being cold in the current snapshot - need to include a next_heating_snapshot field being added from the find_condensing_times_of_cold_gas.py scripy
- Make versions of the radius vs time vs condensed histograms for wake-condensed, bubble-condensed, from-medium condensed etc
- entropy radial plots, try plotting mean and also variance too
- Think about how to distinguish stripped gas from satellites
- Check everythin in np.in1d is a numpy array, otherwise will fail silently
- For histograms, will have condensed mass INCREASE with time due to not having a buffer before starting to record which particles condense
- Measure pre-existing gas that is just drawn out by feedback, rather than condense https://www.youtube.com/watch?v=37i3-3hPtmM
- add cold mass within say 30kpc, 100kpc, 300kpc to KE plot below histogram. Also can do median temperature and percentile bands, etc
- why does the cold gas entropy evolution look so fucked up? - have weirdly large values and artefact-looking regions where the radial profile massively spikes. (galaxies?)
- Pass condensed PIDs to yT to make a projection plot 
- Mass-weighted entropy vs radius divided by median over time at that radius

- make sure no ordering errors like this;


for condensed_i in range(len(condensed_gas_unique_PIDs)):
    if condensed_gas_snapshots_when_cold[condensed_i] <= snapnum_this_rank and condensed_gas_snapshots_when_cold[condensed_i] > snapnum_this_rank - ANALYSE_CONDENSED_PARTS_FROM_PREVIOUS_N_SNAPS:
        condensed_gas_PIDs_to_plot.append(condensed_gas_unique_PIDs[condensed_i])
        condensed_snaps_to_plot.append(condensed_gas_snapshots_when_cold[condensed_i])
    
# sp_condensedpart_mask_for_plot = np.where(np.in1d(sp["PartType0", "unique_PID"].d, condensed_gas_PIDs_to_plot))[0]

sp_condensedpart_mask_for_plot = np.in1d(sp["PartType0", "unique_PID"].d, condensed_gas_PIDs_to_plot)
condensed_gas_x_for_plot = sp["PartType0", "particle_position_x"][sp_condensedpart_mask_for_plot].to("kpc").value - halo_center_post[0].to("kpc").value
condensed_gas_y_for_plot = sp["PartType0", "particle_position_y"][sp_condensedpart_mask_for_plot].to("kpc").value - halo_center_post[1].to("kpc").value
condensed_gas_z_for_plot = sp["PartType0", "particle_position_z"][sp_condensedpart_mask_for_plot].to("kpc").value - halo_center_post[2].to("kpc").value

where e.g. condensed_gas_x_for_plot is used alongside condensed_snaps_to_plot but their particle-ordering is different!!



rescale y-axis of histogram to linearise in ballistic trajectories, to make ballistic particles lie ona straight line, with analytic expression integrated from cluster profile
'''



import PARAMS

DO_BUBBLE_PARTICLE_PLOTS                                    = PARAMS.DO_BUBBLE_PARTICLE_PLOTS                                    
DO_CONDENSED_PARTICLE_PLOTS                                 = PARAMS.DO_CONDENSED_PARTICLE_PLOTS                                 
DO_SATELLITES                                               = PARAMS.DO_SATELLITES                                               
DO_SATELLITE_PLOTS                                          = PARAMS.DO_SATELLITE_PLOTS                                          
DO_MULTIPLOT                                                = PARAMS.DO_MULTIPLOT                                                
PLOT_CONDENSATION_FROM_DIFFERENT_SOURCES_ON_AXIS_1          = PARAMS.PLOT_CONDENSATION_FROM_DIFFERENT_SOURCES_ON_AXIS_1          
USE_EXISTING_FITSFILE_FOR_PLOTS                             = PARAMS.USE_EXISTING_FITSFILE_FOR_PLOTS                             
PLOT_SATELLITES_WITHOUT_WAKES_FOR_TESTING                   = PARAMS.PLOT_SATELLITES_WITHOUT_WAKES_FOR_TESTING                   
PLOT_TO_SHOW_CONDENSING_CLOSE_TO_BUBBLES_AND_WAKE           = PARAMS.PLOT_TO_SHOW_CONDENSING_CLOSE_TO_BUBBLES_AND_WAKE           
PLOT_RADIUS_VS_TEMPERATURE_TESTING_PLOTS                    = PARAMS.PLOT_RADIUS_VS_TEMPERATURE_TESTING_PLOTS                    
DO_ANALYSIS                                                 = PARAMS.DO_ANALYSIS                                                 
PLOT_ANALYSIS                                               = PARAMS.PLOT_ANALYSIS                                               
DO_ENTROPY_PROFILES                                         = PARAMS.DO_ENTROPY_PROFILES                                         
FILL_IN_ANALYSIS_GAPS                                       = PARAMS.FILL_IN_ANALYSIS_GAPS                                       
DO_CLOSE_TO_WAKE_AND_BUBBLES_PLOTS                          = PARAMS.DO_CLOSE_TO_WAKE_AND_BUBBLES_PLOTS                          
DO_INSTANTANEOUSLY_CONDENSED_PLOTS                          = PARAMS.DO_INSTANTANEOUSLY_CONDENSED_PLOTS                          
GO_STRAIGHT_TO_PLOTTING_PRESAVED_ANALYSIS                   = PARAMS.GO_STRAIGHT_TO_PLOTTING_PRESAVED_ANALYSIS                   
RECORD_FULL_HISTORIES_OF_CONDENSING_PARTICLES               = PARAMS.RECORD_FULL_HISTORIES_OF_CONDENSING_PARTICLES               
JUST_DO_RECORD_FULL_HISTORIES_OF_CONDENSING_PARTICLES       = PARAMS.JUST_DO_RECORD_FULL_HISTORIES_OF_CONDENSING_PARTICLES       
FILL_IN_FULL_HISTORIES_OF_CONDENSING_PARTICLES              = PARAMS.FILL_IN_FULL_HISTORIES_OF_CONDENSING_PARTICLES              
kpc_distance_threshold_for_bubble_to_condense_parts         = PARAMS.kpc_distance_threshold_for_bubble_to_condense_parts         
kpc_distance_threshold_for_condensed_to_wake                = PARAMS.kpc_distance_threshold_for_condensed_to_wake                
kpc_distance_threshold_for_condensed_to_satellite_infront   = PARAMS.kpc_distance_threshold_for_condensed_to_satellite_infront   
PLOT_BHWIND_ENERGIES                                        = PARAMS.PLOT_BHWIND_ENERGIES                                        
KE_x_bins_coarsening                                        = PARAMS.KE_x_bins_coarsening                                        
JET_DELTA_V_FROM_QUASARS_inKm_per_s                         = PARAMS.JET_DELTA_V_FROM_QUASARS_inKm_per_s                         
onplot_text_color                                           = PARAMS.onplot_text_color                                           
PLOT_TYPE                                                   = PARAMS.PLOT_TYPE                                                   
PLOT_VMIN                                                   = PARAMS.PLOT_VMIN                                                   
PLOT_VMAX                                                   = PARAMS.PLOT_VMAX                                                   
PLOT_S                                                      = PARAMS.PLOT_S                                                      
PLOT_ALPHA                                                  = PARAMS.PLOT_ALPHA                                                  
MULTI_PLOT_CMAP_IF_CLEAN                                    = PARAMS.MULTI_PLOT_CMAP_IF_CLEAN                                    
MULTIPLOT_DICT                                              = PARAMS.MULTIPLOT_DICT                                              
OVERWRITE_PARENT_PLOT                                       = PARAMS.OVERWRITE_PARENT_PLOT                                      
MIN_BUBBLE_T_K                                              = PARAMS.MIN_BUBBLE_T_K                                             
GET_BH_DETAILS                                              = PARAMS.GET_BH_DETAILS                                             
ANALYSE_BUBBLE_PARTS_FROM_PREVIOUS_N_MYR                    = PARAMS.ANALYSE_BUBBLE_PARTS_FROM_PREVIOUS_N_MYR                   
PLOT_CONDENSED_PARTS_FROM_PREVIOUS_N_MYR                    = PARAMS.PLOT_CONDENSED_PARTS_FROM_PREVIOUS_N_MYR                   
ANALYSE_SATELLITE_WAKE_FROM_PREVIOUS_N_MYR                  = PARAMS.ANALYSE_SATELLITE_WAKE_FROM_PREVIOUS_N_MYR                 
csr_template                                                = PARAMS.csr_template                                               
MIN_DISTANCE_FOR_VALID_SATELLITES_FRAC_OF_R500              = PARAMS.MIN_DISTANCE_FOR_VALID_SATELLITES_FRAC_OF_R500



print("MUST ADD IN DENOISING AND REMOVAL OF CENTRAL PARTS FOR SATELLITE PARTICLES FROM SCRIPT - C !!!!")

import yt 
import numpy as np
import h5py
import caesar
import sys
import socket
import os
from functools import partial
import unyt
import cmasher as cmr
import matplotlib.pyplot as plt
import smplotlib
from MAIN_ANALYSIS_FUNCTIONS_PAPER2 import *
from extract_sim_data import *
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpi4py.util import pkl5
from unyt import parsec as UNYT_PARSEC 
UNYT_KPC = 1000 * UNYT_PARSEC
import matplotlib
from matplotlib.colors import Normalize, SymLogNorm, CSS4_COLORS
from pprint import pprint
# import pygrackle
# from pygrackle import add_grackle_fields
from readgadget import readsnap
# from pygrackle.utilities.data_path import grackle_data_dir
# from pygrackle.utilities.model_tests import model_test_format_version
from datetime import datetime
from astropy.io import fits
from functools import partial
from scipy.spatial import cKDTree
import time
from scipy.stats import binned_statistic

DATE_STR = datetime.today().strftime('%Y-%m-%d')
# DATE_STR = "2024-11-08"

from mpi4py import MPI
''' https://stackoverflow.com/questions/77989495/how-to-mpi-send-and-receive-large-arrays-with-dtype-object-in-python '''
COMM = MPI.COMM_WORLD #pkl5.Intracomm(MPI.COMM_WORLD)  #   MPI.COMM_WORLD
SIZE = COMM.Get_size()
RANK = COMM.Get_rank()





if FILL_IN_FULL_HISTORIES_OF_CONDENSING_PARTICLES:
    RECORD_FULL_HISTORIES_OF_CONDENSING_PARTICLES = True
    JUST_DO_RECORD_FULL_HISTORIES_OF_CONDENSING_PARTICLES = True 









if ALLOW_FILES_MISSING_FOR_TESTING:
    if RANK == 0: print("WARNING! ALLOW_FILES_MISSING_FOR_TESTING is switched on !!")
    time.sleep(10)






ORIGINAL_SNAPNUMS_WHOLE_RANGE = SNAPNUMS_WHOLE_RANGE

if FILL_IN_FULL_HISTORIES_OF_CONDENSING_PARTICLES:
    SNAPNUMS_WHOLE_RANGE = [x for x in SNAPNUMS_WHOLE_RANGE if not os.path.exists(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PRODUCTS/zm{zoomnum}_Level{zoomlevel}_sn{str(x).zfill(3)}_all_condensing_particles_history_data.h5")]
    
elif FILL_IN_ANALYSIS_GAPS:
    SNAPNUMS_WHOLE_RANGE = [x for x in SNAPNUMS_WHOLE_RANGE if not os.path.exists(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PRODUCTS/zm{zoomnum}_Level{zoomlevel}_sn{str(x).zfill(3)}_data.h5")]

# print("len(SNAPNUMS_WHOLE_RANGE)",len(SNAPNUMS_WHOLE_RANGE))

# exit()

NUM_R500_FOR_PROJ = 3
proj_thick_factor = 0.1
CAESAR_HALO_IDX = 0




if DO_MULTIPLOT:
    multi_fig, ((multi_ax1, multi_ax2), (multi_ax3, multi_ax4)) = plt.subplots(2, 2)
    multi_fig.set_size_inches((8,8))
    plt.subplots_adjust(wspace=0, hspace=0)
    for axis_num in [1,2,3,4]:
        if PLOT_TYPE not in MULTIPLOT_DICT[f"AXIS_{axis_num}"]:
            raise NotImplementedError("Multiplot plots must be same as PLOT_TYPE currently")
        



ANALYSE_BUBBLE_PARTS_FROM_PREVIOUS_N_SNAPS = int(int(ANALYSE_BUBBLE_PARTS_FROM_PREVIOUS_N_MYR) / DT)
ANALYSE_SATELLITE_WAKE_FROM_PREVIOUS_N_SNAPS  = int(int(ANALYSE_SATELLITE_WAKE_FROM_PREVIOUS_N_MYR) / DT)
ANALYSE_CONDENSED_PARTS_FROM_PREVIOUS_N_SNAPS = int(int(PLOT_CONDENSED_PARTS_FROM_PREVIOUS_N_MYR) / DT)



if DO_MULTIPLOT:

    DO_BUBBLE_PARTICLE_PLOTS  = True
    DO_CONDENSED_PARTICLE_PLOTS = True
    DO_SATELLITES             = True
    DO_SATELLITE_PLOTS        = True


if JUST_DO_RECORD_FULL_HISTORIES_OF_CONDENSING_PARTICLES:
    DO_SATELLITES = False
    GET_BH_DETAILS = False

SNAPNUMS_WHOLE_RANGE = list(SNAPNUMS_WHOLE_RANGE)
SNAPNUMS_WHOLE_RANGE.sort() ### MUST be sorted
snapnums_this_rank      = np.array_split(SNAPNUMS_WHOLE_RANGE, SIZE)[RANK]
csr_template            = f"/home/b/babul/fjenning/scratch/xray_bubbles/PAPER2_Level{zoomlevel}_zoom_{zoomnum}/Groups/Caesar_FT2Myr_halo_{zoomnum}" + "_{}.hdf5"
csr_file                = f"/home/b/babul/fjenning/scratch/xray_bubbles/PAPER2_Level{zoomlevel}_zoom_{zoomnum}/Caesar_halo_{zoomnum}_{str(RESTART_SNAPNUM).zfill(3)}.hdf5"
RESTART_SNAPFILE        = f"/home/b/babul/fjenning/scratch/xray_bubbles/PAPER2_Level{zoomlevel}_zoom_{zoomnum}/snaps/snap_halo_{zoomnum}_{str(RESTART_SNAPNUM).zfill(3)}.hdf5"



COOLING_TABLE = "CloudyData_UVB=HM2012_shielded.h5"

def _unique_PIDs(field, data):
    # print("data['PartType0','ParticleIDs']", data['PartType0','ParticleIDs'])
    # print("unique PID function", [f"{int(data['PartType0','ParticleIDs'][i].value)}_{int(data['PartType0','ID_Generations'][i].value)}" for i in range(len(data['PartType0','ParticleIDs']))])
    return ( [f"{int(data['PartType0','ParticleIDs'][i].value)}_{int(data['PartType0','ID_Generations'][i].value)}" for i in range(len(data['PartType0','ParticleIDs']))])
'''https://stackoverflow.com/questions/5655574/converting-string-to-ascii-using-ord'''
def convert_to_ascii(text):
    return " ".join(str(ord(char)) for char in text)


def _add_specific_thermal_energy(field, data):
    # print("data['PartType0','ParticleIDs']", data['PartType0','ParticleIDs'])
    # print("unique PID function", [f"{int(data['PartType0','ParticleIDs'][i].value)}_{int(data['PartType0','ID_Generations'][i].value)}" for i in range(len(data['PartType0','ParticleIDs']))])
    return  data['PartType0','specific_thermal_energy']
yt.add_field(
    name=("gas", "specific_thermal_energy"),
    function=_add_specific_thermal_energy,
    sampling_type="local",
    units="J/kg",
)


def _xray_entropy(field, data):
    ''' https://ned.ipac.caltech.edu/level5/Sept09/Bohringer/Bohringer3.html'''
    return unyt.physical_constants.kb * data["gas", "temperature"] / data["gas", "El_number_density"]**(2/3)


yt.add_field(
    name=("gas", "xray_entropy"),
    function=_xray_entropy,
    sampling_type="local",
    units="keV* cm**2",
)


def _pressure(field, data):
    gamma = 5/3
    return (
        (gamma - 1.0)
        * data["gas", "density"]
        * data["gas", "specific_thermal_energy"]
    )
    
yt.add_field(
    name=("gas", "pressure"),
    function=_pressure,
    sampling_type="local",
    units="dyne/cm**2",
)


######## MAIN ##################


os.makedirs(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/", exist_ok=True)
os.makedirs(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/testing_PNGs_for_combined_data/", exist_ok=True)


if USE_CAESAR:
    if RANK == 0:
        csr_obj = caesar.load(csr_file)
        csr_massive_halo = csr_obj.halos[CAESAR_HALO_IDX]
        halo_r500   = csr_massive_halo.virial_quantities["r500c"]
        halo_m500   = csr_massive_halo.virial_quantities["m500c"]
        
        
        # print(csr_massive_halo.central_galaxy.GroupID)
        # print(csr_obj.galaxies[0].GroupID)
        # exit()
        
        # central_gal = [x for x in csr_obj.galaxies if x.GroupID == csr_massive_halo.central_galaxy][0]
        central_gal_stellar_mass = csr_massive_halo.central_galaxy.masses['stellar']
        
        print('%.2f' % float(np.log10(central_gal_stellar_mass.to("Msun").value)))
        print('%.2f' % float(np.log10(halo_m500.to("Msun").value)))
        # exit()
        halo_center = csr_massive_halo.minpotpos
        halo_center_original = csr_massive_halo.minpotpos

        ds_fid = yt.load(RESTART_SNAPFILE)
        start_time = ds_fid.current_time.to("Myr").value
        start_scale_factor = ds_fid.scale_factor
        sp_fid = ds_fid.sphere(center = halo_center , radius = halo_r500)
        SMBH_pid = int(sp_fid["PartType5", "particle_index"][np.argmax(sp_fid["PartType5", "BH_Mass"])])
    else:
        halo_r500       = None
        start_time      = None
        SMBH_pid        = None
        halo_center_original = None
        start_scale_factor = None
    halo_r500       = COMM.bcast(halo_r500       , root=0)
    start_time      = COMM.bcast(start_time      , root=0)
    SMBH_pid        = COMM.bcast(SMBH_pid        , root=0)
    halo_center_original = COMM.bcast(halo_center_original, root=0)
    start_scale_factor = COMM.bcast(start_scale_factor, root = 0)
    print(RANK, SMBH_pid)

else:
    raise NotImplementedError()
    # ds_fid = yt.load(RESTART_SNAPFILE)
    # # ad_fid = ds_fid.all_data()
    # start_time = ds_fid.current_time.to("Myr").value
    # # SMBH_pid = int(ad_fid["PartType5", "particle_index"][np.argmax(ad_fid["PartType5", "BH_Mass"])])
    # SMBH_pid = 314452702
    # halo_r500 = 1000 * 1000 * UNYT_PARSEC

print("halo R500", halo_r500)
bubble_gas_unique_PIDS_dict_list = []
condensed_between_snaps_cold_gas_unique_PIDS_dict_list = []
cold_gas_unique_PIDS_dict_list = []


grackle_data_file = os.path.join(RESTART_SNAPFILE.split("/snaps/")[0], COOLING_TABLE)

# print(pygrackle.yt_fields._field_map)


# exit()
grackle_pars = {'grackle_data_file': grackle_data_file,
                'UVbackground': 1,
                'h2_on_dust': 1}
cuts_dict    = [{'field':('gas','density'), '<':(0.1*unyt.physical_constants.proton_mass)*unyt.cm**-3},
                    {'field':("gas", "temperature"), '>': 2e5*unyt.K},
                    {'field':("PartType0", 'Density'), '<':(0.1*unyt.physical_constants.proton_mass)*unyt.cm**-3},
                    {'field':("PartType0", "temperature"), '>': 2e5*unyt.K},
                    {'field':("PartType0", "DelayTime"), '=':0},
                    {'field':("PartType0", "StarFormationRate"), '=':0}
                    ]
def _filtered_gas(pfilter, data, dataset_cuts):
    pfilter = True
    for i, cut in enumerate(dataset_cuts):
        if cut.get("=",None) != None:
            pfilter &= data[cut["field"]] == cut["="]
        else:
            if cut.get(">",None) != None:
                pfilter &= data[cut["field"]] > cut[">"]
            if cut.get("<",None) != None:
                pfilter &= data[cut["field"]] < cut["<"]  
                        
    return pfilter

yt.add_particle_filter(
    "filtered_gas",
    function=partial(_filtered_gas, dataset_cuts = cuts_dict),
    filtered_type="PartType0",
    requires=['density', 'temperature', 'DelayTime', 'StarFormationRate'],
)


if PLOT_BHWIND_ENERGIES and len(snapnums_this_rank) == 0:
    ### just load in one snapfile so we can define cosmology
    snapnums_this_rank = [RESTART_SNAPNUM+3,]
for snapnum_this_rank in snapnums_this_rank:
    
    
    print(RANK, " Doing Snapnum ", snapnum_this_rank)
    if GO_STRAIGHT_TO_PLOTTING_PRESAVED_ANALYSIS and not PLOT_BHWIND_ENERGIES:
        break
    

    
    snapfile = None

    ''' Load in all needed snapshots. The current one (post), the one we want to look at pre-condensing, and the one we want to look at pre-bubble''' 
    if int(snapnum_this_rank) > RESTART_SNAPNUM:
        snapfile            = f"/home/b/babul/fjenning/scratch/xray_bubbles/PAPER2_Level{zoomlevel}_zoom_{zoomnum}/snaps/FT2Myr_snap_halo_{zoomnum}_{str(snapnum_this_rank).zfill(3)}.hdf5"
    else:
        snapfile            = f"/home/b/babul/fjenning/scratch/xray_bubbles/PAPER2_Level{zoomlevel}_zoom_{zoomnum}/snaps/snap_halo_{zoomnum}_{str(snapnum_this_rank).zfill(3)}.hdf5"
        
        
    ''' We need to load the previous snapfile in order to remove spurious condensed particles which just enter the sphere at snapnum_this_rank but were outside the sphere (and hence appeared as 'hot' in the condense-snap script) in snapnum_this_rank - 1'''
    if int(snapnum_this_rank-DELTA_SNAP_FOR_CONDENSATION) > RESTART_SNAPNUM:
        snapfile_pre_for_infall_checking            = f"/home/b/babul/fjenning/scratch/xray_bubbles/PAPER2_Level{zoomlevel}_zoom_{zoomnum}/snaps/FT2Myr_snap_halo_{zoomnum}_{str(snapnum_this_rank-DELTA_SNAP_FOR_CONDENSATION).zfill(3)}.hdf5"
    else:
        snapfile_pre_for_infall_checking            = f"/home/b/babul/fjenning/scratch/xray_bubbles/PAPER2_Level{zoomlevel}_zoom_{zoomnum}/snaps/snap_halo_{zoomnum}_{str(snapnum_this_rank-DELTA_SNAP_FOR_CONDENSATION).zfill(3)}.hdf5"
        
        


    
    ds = yt.load(snapfile)
    
    
    ''' Get the PIDs of wind particles from the BH, so we can filter out stellar winds for the bubble-tagging'''
    H0      =   100*ds.hubble_constant
    omlam0  =   ds.omega_lambda
    om0     =   ds.omega_matter
    Tcmb0   =   2.7 
    cosmo   =   FlatLambdaCDM(H0=H0, Om0 = om0, Ob0 = 0.048, Tcmb0=Tcmb0, )
    define_units(cosmo)
    
    if GO_STRAIGHT_TO_PLOTTING_PRESAVED_ANALYSIS:
        break
    
    redshift_post = ds.current_redshift
    age_post = cosmo.age(redshift_post)
    dminus = 20*ASTROPY_Myr
    age_post_minus_20Myr = age_post - (20*ASTROPY_Myr)
    redshift_post_minus_20Myr = z_at_value(cosmo.age, age_post_minus_20Myr)  
    
    scale_fac_post             = 1./(1.+redshift_post)
    scale_fac_post_minus_20Myr = 1./(1.+redshift_post_minus_20Myr)



    ds.add_field(
        name=("PartType0", "unique_PID"),
        function=_unique_PIDs,
        sampling_type="local",
        units=None,
        dimensions = None,
    )  
    
    # ds_pre_for_infall_checking = yt.load(snapfile_pre_for_infall_checking)
    # ds_pre_for_infall_checking.add_field(
    #     name=("PartType0", "unique_PID"),
    #     function=_unique_PIDs,
    #     sampling_type="local",
    #     units=None,
    #     dimensions = None,
    # )  



    ds.add_particle_filter("filtered_gas",)
    if RANK == 0: print("Added filtered gas filter", flush = True)
    ad = ds.all_data()
    halo_center_post = ad["PartType5", "particle_position"][ad["PartType5", "particle_index"] == SMBH_pid][0]
    
    
    sp  = ds.sphere(center = halo_center_post , radius = (NUM_R500_FOR_PROJ*halo_r500.to("kpc").value,'kpc'))
    sp_bulk_vel = sp.quantities.bulk_velocity()
    sp.set_field_parameter("bulk_velocity", sp_bulk_vel)  ### MUST do this before accessing any radial velocity fields
    
    
    # print(sp["gas","entropy"][0])
    # exit()
    
    # ad_pre = ds_pre_for_infall_checking.all_data()
    # halo_center_pre = ad_pre["PartType5", "particle_position"][ad_pre["PartType5", "particle_index"] == SMBH_pid][0]
    # sp_pre_for_infall_checking  = ds_pre_for_infall_checking.sphere(center = halo_center_pre , radius = (NUM_R500_FOR_PROJ*halo_r500.to("kpc").value,'kpc'))

    
    
    
    # print(ds.ActiveDimensions)
    # exit()
    
    
    # ''' get segmentation faults if try to add to both then use them later on'''
    # add_grackle_fields(ds_pre_for_condensegas, parameters=grackle_pars)
    
    
    # add_grackle_fields(ds, parameters=grackle_pars)
    
    # ad_pre = ds_pre_for_condensegas.all_data()
    
    
    
    
    
    '''########################################'''
    '''                                        '''
    '''                                        '''
    '''  Read Black Hole Feedback Time Series  '''
    '''                                        '''
    '''                                        '''
    '''########################################'''
    if GET_BH_DETAILS:
        ''' Load the BH wind and properties histories'''
        bhdetails_dir = snapfile.split("/snaps/")[0] + "/snaps/blackhole_details/"


        bh_wind_files = [bhdetails_dir + f"/{f}" for f in os.listdir(bhdetails_dir) if "bhwind" in f]
        bh_details_files = [bhdetails_dir + f"/{f}" for f in os.listdir(bhdetails_dir) if "blackhole_details" in f]  
        bhwind_history = get_bhwind_history(SMBH_pid,  files= bh_wind_files, min_a = scale_fac_post_minus_20Myr, max_a = scale_fac_post)
        
        
        if bhwind_history == []:
            print("No BH winds this snap")
            bhwind_PIDs = []
        else:
            test_units=True
            bhinfo_history  = get_bh_info(SMBH_pid,  files=bh_details_files, min_a = scale_fac_post_minus_20Myr, max_a = scale_fac_post, a_list_winds = bhwind_history["a_list_winds"], )
        
            if RUN_TESTS: check_dict_units(bhwind_history, "bhwind_history")
            if RUN_TESTS: check_dict_units(bhinfo_history, "bhinfo_history")
            if scale_fac_post < max(bhwind_history["a_list_winds"]) or scale_fac_post_minus_20Myr > min(bhwind_history["a_list_winds"]):
                raise ValueError("Scale factor range for winds is bad")
            
            bhwind_PIDs = [int(id) for id in bhwind_history["bhwind_id_history"]]



    '''#############################'''
    '''                             '''
    '''                             '''
    '''   Read Satellite Particles  '''
    '''  and generate list of wake  '''
    '''and in-front PIDs from snaps '''
    '''    previous to current      '''
    '''                             '''
    '''                             '''
    '''#############################'''
    
    

        
    if DO_SATELLITES:
        if PLOT_SATELLITES_WITHOUT_WAKES_FOR_TESTING:
            
            
            sat_test_fig, sat_test_axis  = plt.subplots()
            sat_test_fig.set_size_inches((10,10))
            
            buff_size = (800,800)
            raw_fields  = [("filtered_gas", "temperature"),]
            sat_test_extent = [(1200,'kpc'),(1200,'kpc'),(20,'kpc')]
            prj_sat_test         = yt.ProjectionPlot(ds, normal = "z",  fields=raw_fields, weight_field=("filtered_gas", "density"),center = halo_center_post.to("kpc"), width = sat_test_extent, buff_size = buff_size)
            # prj_sat_test.save(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/testing_PNGs_for_combined_data/TEST_SATELLITE_POS:zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}.png")
            # # plt.close()
            
            # exit()
            sat_test_cmap = "Greys_r"
            sat_test_img_data = np.asarray(prj_sat_test.frb[raw_fields[0]])
            half_width_x =    sat_test_extent[0][0]    /2            
            half_width_y =    sat_test_extent[1][0]    /2
            
            # print("sat_test_img_data", sat_test_img_data, np.nanmax(sat_test_img_data))
           
            im = sat_test_axis.imshow(sat_test_img_data, cmap = sat_test_cmap, #cmap="cmr.apple",
            extent = [-half_width_x, half_width_x, -half_width_y, half_width_y],
            # norm=LogNorm(vmin=1e5, vmax=1e7,))
            # norm=LogNorm())
            )
        t_start_sat = time.time()
        wake_PIDS                   = []
        infront_of_satellite_PIDs   = []
        
        for sat_snap in range(snapnum_this_rank - ANALYSE_SATELLITE_WAKE_FROM_PREVIOUS_N_SNAPS, snapnum_this_rank+1):
            
            if sat_snap < RESTART_SNAPNUM:
                continue
            if sat_snap > snapnum_this_rank:
                continue
            
            try:
                '''File is Written in find_and_save_cold_gas_parts.py'''
                with h5py.File(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_DATA/zm{zoomnum}_Level{zoomlevel}_satellite_data_sn{sat_snap}.h5",'r') as satellite_f:
                    
                    ''' These are the satellite caesar indexes'''
                    sat_keys = [x.split("/")[-1] for x in list(satellite_f[f"SATELLITE/{sat_snap}/"].keys())]
                    sat_keys.sort()

                    if '0' in sat_keys:
                        raise ValueError("Central should not be in list")

                    for galaxy_i in sat_keys:
                        this_sat_NN_infront_gas_IDs                 = satellite_f[f"SATELLITE/{sat_snap}/{galaxy_i}/in_front_COLD_IDs"]
                        this_sat_NN_infront_gas_unique_identifiers  = satellite_f[f"SATELLITE/{sat_snap}/{galaxy_i}/in_front_COLD_Unique_Identifier"]
                        this_sat_NN_behind_gas_IDs                  = satellite_f[f"SATELLITE/{sat_snap}/{galaxy_i}/behind_COLD_IDs"]
                        this_sat_NN_behind_gas_unique_identifiers   = satellite_f[f"SATELLITE/{sat_snap}/{galaxy_i}/behind_COLD_Unique_Identifier" ]
                        this_sat_vel_rel_to_group_com               = satellite_f[f"SATELLITE/{sat_snap}/{galaxy_i}/gal_velocity_rel_to_group_com"]
                        this_sat_pos                                = satellite_f[f"SATELLITE/{sat_snap}/{galaxy_i}/position"]


                        ''' Don;t include this satellite particles if a significant number of particles lie within the center. Cannot do later when we combine all satellite PIDs, because strong outflows can escape the center and these particles will still be classified, erroneously, as satellite parts. Also, we must be careful here not to discount satellite which just have some noise in the center region.'''
                        this_sat_mask    = (np.in1d(sp["PartType0", "unique_PID"].d, this_sat_NN_infront_gas_unique_PIDs))|(np.in1d(sp["PartType0", "unique_PID"].d, this_sat_NN_behind_gas_unique_PIDs))
                        print('sum(sp_post["PartType0", "unique_PID"][this_sat_mask].to("kpc").value < MIN_RADIUS_FOR_VALID_SATELLITES_FRAC_OF_R500*halo_r500.to("kpc").value )', sum(sp["PartType0", "particle_radius"][this_sat_mask].to("kpc").value < MIN_RADIUS_FOR_VALID_SATELLITES_FRAC_OF_R500*halo_r500.to("kpc").value ))
                        if sum(sp_post["PartType0", "particle_radius"][this_sat_mask].to("kpc").value < MIN_RADIUS_FOR_VALID_SATELLITES_FRAC_OF_R500*halo_r500.to("kpc").value ) > 5:
                            continue 
                        # print("this_sat_pos[0] - sp.center[0].to('kpc').value", this_sat_pos[0],this_sat_pos[0] - sp.center[0].to('kpc').value)
                        # print("this_sat_pos[1] - sp.center[1].to('kpc').value", this_sat_pos[1], this_sat_pos[1] - sp.center[1].to('kpc').value)
                        if PLOT_SATELLITES_WITHOUT_WAKES_FOR_TESTING and sat_snap == snapnum_this_rank:
                            sat_test_axis.scatter(x=(this_sat_pos[0] - halo_center_post[0].to('kpc').value), 
                                              y=(this_sat_pos[1] - halo_center_post[1].to('kpc').value),
                                              transform = sat_test_axis.transData,
                                              color = "red")
                        
                        
                        # CAESARFILE = csr_template.format(str(sat_snap).zfill(3))
                        # this_csr_obj = caesar.load(CAESARFILE)
                        # this_sat = this_csr_obj.galaxies[int(galaxy_i)]
                        # this_sat_pos_from_csr = this_sat.pos.to("kpc").value
                        
                        # print("this_sat_pos_from_csr/this_sat_pos_from_file", this_sat_pos_from_csr/this_sat_pos)
                        # print(snapnum_this_rank, "this_sat_pos rel to center", '%.2f' % float(np.linalg.norm(this_sat_pos - halo_center_post.to("kpc").value)))
                            
                        this_sat_NN_infront_gas_unique_PIDs         =   ["{}_{}".format(a_, b_) for a_, b_ in zip(this_sat_NN_infront_gas_IDs, this_sat_NN_infront_gas_unique_identifiers)]
                        this_sat_NN_behind_gas_unique_PIDs          =   ["{}_{}".format(a_, b_) for a_, b_ in zip(this_sat_NN_behind_gas_IDs, this_sat_NN_behind_gas_unique_identifiers)]
            
            
                        if len(this_sat_NN_infront_gas_unique_PIDs) == 0 or len(this_sat_NN_behind_gas_unique_PIDs) == 0:
                            if PRINT_WARNING_MESSAGES: print("Either no wake or no in-front of satellite gas, probably a spurious halo finder issue so skipping!")
                            continue
                        if len(this_sat_NN_infront_gas_unique_PIDs) + len(this_sat_NN_behind_gas_unique_PIDs) < 3:
                            if PRINT_WARNING_MESSAGES: print("< 3 particles in wake + in-front satellite parts, so skipping...!")
                            continue       

                        wake_PIDS.extend(this_sat_NN_behind_gas_unique_PIDs)
                        
                        if sat_snap == snapnum_this_rank:
                            infront_of_satellite_PIDs.extend(this_sat_NN_infront_gas_unique_PIDs)
                        else:
                            wake_PIDS.extend(this_sat_NN_infront_gas_unique_PIDs) ### These particles SHOULD all have moved into the wake by the present snapshot
            except Exception as e:
                if not ALLOW_FILES_MISSING_FOR_TESTING:
                    raise RuntimeError(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_DATA/zm{zoomnum}_Level{zoomlevel}_satellite_data_sn{sat_snap}.h5 \n",e)   
                else:
                    if PRINT_WARNING_MESSAGES: print(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_DATA/zm{zoomnum}_Level{zoomlevel}_satellite_data_sn{sat_snap}.h5 \n",e, "\n continuing since ALLOW_FILES_MISSING_FOR_TESTING is switched on")
        
        
        if PLOT_SATELLITES_WITHOUT_WAKES_FOR_TESTING:
            plt.figure(sat_test_fig)
            plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/testing_PNGs_for_combined_data/TEST_SATELLITE_POS:zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}.png")
            plt.close()
            
            
            
        ''' Create satellite masks'''
        in_front_MASK    = (np.in1d(sp["PartType0", "unique_PID"].d, infront_of_satellite_PIDs))&(~np.in1d(sp["PartType0", "unique_PID"].d, wake_PIDS))
        wake_MASK        = (~np.in1d(sp["PartType0", "unique_PID"].d, infront_of_satellite_PIDs))&(np.in1d(sp["PartType0", "unique_PID"].d, wake_PIDS))
        
        if PRINT_WARNING_MESSAGES: print(" - if not equal then we have duplicate id's, might explain far-away particles in wake", len(sp["PartType0", "unique_PID"].d[wake_MASK]), len(set(sp["PartType0", "unique_PID"].d[wake_MASK])))
        
        
        test_satellite_actually_unique_PIDs, test_satellite_unique_indexing,test_satellite_unique_inverse, test_satellite_counts = np.unique(sp["PartType0","unique_PID"].d,return_index =True, return_inverse=True,
                        return_counts=True)
        ''' Get an array same length as all the UniquePIDs with entries corresponding to the number of times that PID at that index appears'''
        test_satellite_num_counts_size_of_orig_arr = test_satellite_counts[test_satellite_unique_inverse]
        test_satellite_unique_PIDs_mask = np.asarray(test_satellite_num_counts_size_of_orig_arr) == 1
        in_front_MASK &= (test_satellite_unique_PIDs_mask)
        wake_MASK &= (test_satellite_unique_PIDs_mask)
        
        
        # exit()
    
        in_front_gas_x = sp["PartType0", "particle_position_x"][in_front_MASK].to("kpc").value - halo_center_post[0].to("kpc").value
        in_front_gas_y = sp["PartType0", "particle_position_y"][in_front_MASK].to("kpc").value - halo_center_post[1].to("kpc").value
        in_front_gas_z = sp["PartType0", "particle_position_z"][in_front_MASK].to("kpc").value - halo_center_post[2].to("kpc").value
    
        wake_gas_x = sp["PartType0", "particle_position_x"][wake_MASK].to("kpc").value - halo_center_post[0].to("kpc").value
        wake_gas_y = sp["PartType0", "particle_position_y"][wake_MASK].to("kpc").value - halo_center_post[1].to("kpc").value
        wake_gas_z = sp["PartType0", "particle_position_z"][wake_MASK].to("kpc").value - halo_center_post[2].to("kpc").value
            
            
            
        satellite_wake_tree                 = cKDTree(np.vstack((wake_gas_x,wake_gas_y,wake_gas_z)).T)
        satellite_infront_tree              = cKDTree(np.vstack((in_front_gas_x,in_front_gas_y,in_front_gas_z)).T)
        ''' Find and root out "noisy" particles which are loners and may have been ejected etc'''
        list_of_neighbours = [len(x)-1 for x in satellite_wake_tree.query_ball_tree(satellite_wake_tree, r = MIN_DIST_FOR_ISOLATED_SATELLITE_PARTS)]  ### minus 1 to remove self-counting
        print("min list_of_neighbours", min(list_of_neighbours), list_of_neighbours[:10])
        print("list_of_neighbours", np.shape(list_of_neighbours))
        wake_noise_MASK = (np.asarray(list_of_neighbours)>MIN_NUM_NEIGHBOURS_WITHIN_MIN_DIST_FOR_ISOLATED_SATELLITE_PARTS)
        
        list_of_neighbours = [len(x)-1 for x in satellite_infront_tree.query_ball_tree(satellite_infront_tree, r = MIN_DIST_FOR_ISOLATED_SATELLITE_PARTS)]
        print("min list_of_neighbours", min(list_of_neighbours))
        in_front_noise_MASK = (np.asarray(list_of_neighbours)>MIN_NUM_NEIGHBOURS_WITHIN_MIN_DIST_FOR_ISOLATED_SATELLITE_PARTS)
        
        ''' Get curent positions of all satellite wake and in-front particles, then build trees'''
        in_front_gas_x = sp["PartType0", "particle_position_x"][in_front_MASK][in_front_noise_MASK].to("kpc").value - halo_center_post[0].to("kpc").value
        in_front_gas_y = sp["PartType0", "particle_position_y"][in_front_MASK][in_front_noise_MASK].to("kpc").value - halo_center_post[1].to("kpc").value
        in_front_gas_z = sp["PartType0", "particle_position_z"][in_front_MASK][in_front_noise_MASK].to("kpc").value - halo_center_post[2].to("kpc").value
    
        wake_gas_x = sp["PartType0", "particle_position_x"][wake_MASK][wake_noise_MASK].to("kpc").value - halo_center_post[0].to("kpc").value
        wake_gas_y = sp["PartType0", "particle_position_y"][wake_MASK][wake_noise_MASK].to("kpc").value - halo_center_post[1].to("kpc").value
        wake_gas_z = sp["PartType0", "particle_position_z"][wake_MASK][wake_noise_MASK].to("kpc").value - halo_center_post[2].to("kpc").value
        
            
            
            
            
            
            
            
            
        print(RANK, "Time to do satellites:", time.time() - t_start_sat)        
    
    
    
    
    
    '''##################################'''
    '''                                  '''
    '''                                  '''
    '''Load All Condensing particles data'''
    '''                                  '''
    '''                                  '''
    '''##################################'''
    t_start_condense_block = time.time()
    ''' ANALYISIS WITH DATA AT THIS SNAPNUM'''
    with h5py.File(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_CONDENSING_HISTORY/data/CONDENSING_DATA.h5",'r') as f:
        condensed_gas_IDs                            =   list(f[f"CONDENSING_GAS/COLD_IDs"] )
        condensed_gas_unique_identifiers             =   list(f[f"CONDENSING_GAS/COLD_Unique_Identifier"] )
        condensed_gas_snapshots_when_cold            =   list(f[f"CONDENSING_GAS/condensing_snaps"])
        condensed_gas_snapshots_when_hot_precondense =   [sn-1 for sn in condensed_gas_snapshots_when_cold]
        condensed_gas_unique_PIDs                   =   ["{}_{}".format(a_, b_) for a_, b_ in zip(condensed_gas_IDs, condensed_gas_unique_identifiers)]


    if RECORD_FULL_HISTORIES_OF_CONDENSING_PARTICLES: 
        PIDs_gas_that_ever_condensed  = [condensed_gas_unique_PIDs[i] for i in range(len(condensed_gas_unique_PIDs))]
        
        ''' Some particles may condense more than once'''
        PIDs_gas_that_ever_condensed = list(set(PIDs_gas_that_ever_condensed))
        ''' Sort here so that each data file will be correctly/easily ordered'''
        
        print("PIDs_gas_that_ever_condensed len1", len(PIDs_gas_that_ever_condensed))
        PIDs_gas_that_ever_condensed.sort()
        # print("PIDs_gas_that_ever_condensed len2", len(PIDs_gas_that_ever_condensed))
        
    ''' THESE TWO PIDS ARE FOR TWO SEPARTE CONDENSING INTERVALS'''
    PIDs_gas_immediately_pre_condensing  = [condensed_gas_unique_PIDs[i] for i in range(len(condensed_gas_unique_PIDs)) if condensed_gas_snapshots_when_hot_precondense[i] == snapnum_this_rank ]
    PIDs_gas_immediately_post_condensing = [condensed_gas_unique_PIDs[i] for i in range(len(condensed_gas_unique_PIDs)) if condensed_gas_snapshots_when_cold[i] == snapnum_this_rank]
    
    
    
    
    
    
    
    # if RECORD_FULL_HISTORIES_OF_CONDENSING_PARTICLES: 
        
        
    #     ''' Make a mask for the sphere which will select out the particles in PIDs_gas_that_ever_condensed that have unique PIDs in sp'''
    #     gas_that_ever_condensed_mask_this_snap = np.in1d(sp["PartType0", "unique_PID"].d, PIDs_gas_that_ever_condensed)
    #     # print("PIDs_gas_that_ever_condensed len3", len(PIDs_gas_that_ever_condensed))
    #     condensed_ever_actually_unique_PIDs, condensed_ever_unique_indexing,condensed_ever_unique_inverse, condensed_ever_counts = np.unique(sp["PartType0","unique_PID"].d,return_index =True, return_inverse=True,
    #                     return_counts=True)
    #     ''' Get an array same length as all the UniquePIDs with entries corresponding to the number of times that PID at that index appears'''
    #     condensed_ever_num_counts_size_of_orig_arr = condensed_ever_counts[condensed_ever_unique_inverse]
    #     condensed_ever_unique_PIDs_mask = np.asarray(condensed_ever_num_counts_size_of_orig_arr) == 1
    #     gas_that_ever_condensed_mask_this_snap &= condensed_ever_unique_PIDs_mask  ### remove duplicated condensed gas PIDs in this snap
    #     # sp_condensedpart_mask_for_plot &= condensed_unique_PIDs_mask
    #     parts_condensed_ever_unique_PIDs_this_snap      = sp["PartType0", "unique_PID"].value[gas_that_ever_condensed_mask_this_snap]
    #     parts_condensed_ever_mass_this_snap             = sp["PartType0", "particle_mass"].to("Msun").value[gas_that_ever_condensed_mask_this_snap]
    #     parts_condensed_ever_radial_velocity_this_snap  = sp["PartType0", "radial_velocity"].to("km/s").value[gas_that_ever_condensed_mask_this_snap]
    #     parts_condensed_ever_radius_this_snap           = sp["PartType0", "particle_radius"].to("kpc").value[gas_that_ever_condensed_mask_this_snap]
    #     parts_condensed_ever_xray_entropy_this_snap     = sp["gas", "xray_entropy"].to("keV* cm**2").value[gas_that_ever_condensed_mask_this_snap]
    #     parts_condensed_ever_pressure_this_snap         = sp["gas", "pressure"].to("dyne/cm**2").value[gas_that_ever_condensed_mask_this_snap]
    #     parts_condensed_ever_temperature_this_snap      = sp["PartType0", "temperature"].to("K").value[gas_that_ever_condensed_mask_this_snap]
        
    
    #     print("len parts_condensed_ever_unique_PIDs_this_snap 1", len(parts_condensed_ever_unique_PIDs_this_snap))
    
    #     ''' Some IDs that ever condense may not be in the sphere yet. Alternatively they might have duplicate PIDs in this snap in sp. Include these as nans in the above array to be saved, to preserve the uniformity of the PIDs list across all saved files'''
    #     ids_missing_from_sp = np.setdiff1d(PIDs_gas_that_ever_condensed,parts_condensed_ever_unique_PIDs_this_snap)    # yields the elements in 1st arg that are NOT in 2nd arg
    #     len_ids_missing_from_sp = np.size(ids_missing_from_sp)

    #     array_of_nans_len_ids_missing_from_sp = np.empty(len_ids_missing_from_sp)
    #     array_of_nans_len_ids_missing_from_sp[:] = np.nan
    #     parts_condensed_ever_unique_PIDs_this_snap      = np.concatenate((parts_condensed_ever_unique_PIDs_this_snap      , ids_missing_from_sp))
    #     parts_condensed_ever_mass_this_snap             = np.concatenate((parts_condensed_ever_mass_this_snap,array_of_nans_len_ids_missing_from_sp ))
    #     parts_condensed_ever_radial_velocity_this_snap  = np.concatenate((parts_condensed_ever_radial_velocity_this_snap  , array_of_nans_len_ids_missing_from_sp))
    #     parts_condensed_ever_radius_this_snap           = np.concatenate((parts_condensed_ever_radius_this_snap           , array_of_nans_len_ids_missing_from_sp))
    #     parts_condensed_ever_xray_entropy_this_snap     = np.concatenate((parts_condensed_ever_xray_entropy_this_snap     , array_of_nans_len_ids_missing_from_sp))
    #     parts_condensed_ever_pressure_this_snap         = np.concatenate((parts_condensed_ever_pressure_this_snap    , array_of_nans_len_ids_missing_from_sp))
    #     parts_condensed_ever_temperature_this_snap      = np.concatenate((parts_condensed_ever_temperature_this_snap, array_of_nans_len_ids_missing_from_sp))
        
    #     sort_array_for_condensed_ever_unique_PIDs = np.argsort(parts_condensed_ever_unique_PIDs_this_snap)

    #     parts_condensed_ever_unique_PIDs_this_snap      =   parts_condensed_ever_unique_PIDs_this_snap[sort_array_for_condensed_ever_unique_PIDs]
    #     parts_condensed_ever_mass_this_snap             =   parts_condensed_ever_mass_this_snap[sort_array_for_condensed_ever_unique_PIDs]
    #     parts_condensed_ever_radial_velocity_this_snap  =   parts_condensed_ever_radial_velocity_this_snap[sort_array_for_condensed_ever_unique_PIDs]
    #     parts_condensed_ever_radius_this_snap           =   parts_condensed_ever_radius_this_snap[sort_array_for_condensed_ever_unique_PIDs]
    #     parts_condensed_ever_xray_entropy_this_snap     =   parts_condensed_ever_xray_entropy_this_snap[sort_array_for_condensed_ever_unique_PIDs]
    #     parts_condensed_ever_pressure_this_snap         =   parts_condensed_ever_pressure_this_snap[sort_array_for_condensed_ever_unique_PIDs]
    #     parts_condensed_ever_temperature_this_snap      =   parts_condensed_ever_temperature_this_snap[sort_array_for_condensed_ever_unique_PIDs]
        
    #     # print("len(parts_condensed_ever_unique_PIDs_this_snap), len(set(parts_condensed_ever_unique_PIDs_this_snap))", len(parts_condensed_ever_unique_PIDs_this_snap), len(set(parts_condensed_ever_unique_PIDs_this_snap)))
        
        
    #     ''' Make sure that for each snapshot, the arrays are saved in the identical order, so we can efficiently read them at the end.'''
    #     assert(np.array_equal(parts_condensed_ever_unique_PIDs_this_snap, PIDs_gas_that_ever_condensed))
    #     # print("np.array_equal(parts_condensed_ever_unique_PIDs_this_snap, PIDs_gas_that_ever_condensed)", np.array_equal(parts_condensed_ever_unique_PIDs_this_snap, PIDs_gas_that_ever_condensed))
    #     # exit()
        
    #     ''' Save data for all PIDs which ever condense, for this snapshot'''
    #     os.makedirs(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PRODUCTS/", exist_ok=True)
    #     f = h5py.File(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PRODUCTS/zm{zoomnum}_Level{zoomlevel}_sn{str(snapnum_this_rank).zfill(3)}_all_condensing_particles_history_data.h5",'w')
        
        
    #     f.create_dataset(f"PARTICLES_THAT_EVER_CONDENSE/{snapnum_this_rank}/masses", data=parts_condensed_ever_mass_this_snap)
    #     f.create_dataset(f"PARTICLES_THAT_EVER_CONDENSE/{snapnum_this_rank}/radial_velocity",  data=parts_condensed_ever_radial_velocity_this_snap)
    #     f.create_dataset(f"PARTICLES_THAT_EVER_CONDENSE/{snapnum_this_rank}/radius", data = parts_condensed_ever_radius_this_snap)
    #     f.create_dataset(f"PARTICLES_THAT_EVER_CONDENSE/{snapnum_this_rank}/xray_entropy", data = parts_condensed_ever_xray_entropy_this_snap)
    #     f.create_dataset(f"PARTICLES_THAT_EVER_CONDENSE/{snapnum_this_rank}/pressure", data = parts_condensed_ever_pressure_this_snap )
    #     f.create_dataset(f"PARTICLES_THAT_EVER_CONDENSE/{snapnum_this_rank}/temperature", data = parts_condensed_ever_temperature_this_snap )
        
        
    #     f.create_dataset(f"PARTICLES_THAT_EVER_CONDENSE/{snapnum_this_rank}/IDs", data =[int(x.split("_")[0]) for x in parts_condensed_ever_unique_PIDs_this_snap])
    #     f.create_dataset(f"PARTICLES_THAT_EVER_CONDENSE/{snapnum_this_rank}/Unique_Identifier", data =[int(x.split("_")[1]) for x in parts_condensed_ever_unique_PIDs_this_snap])
    #     f.close()
    
    
    
    
    
    
    
    
    
    
    
    if JUST_DO_RECORD_FULL_HISTORIES_OF_CONDENSING_PARTICLES:
        continue
    if RANK == 0: print("Number of condensing particles", len(condensed_gas_IDs), flush = True)
    # exit()













    
    if len(PIDs_gas_immediately_pre_condensing) > 0:
        
        sp_pre_condensing_mask_for_analysis  = np.in1d(sp["PartType0", "unique_PID"].d, PIDs_gas_immediately_pre_condensing)
        sp_post_condensing_mask_for_analysis = np.in1d(sp["PartType0", "unique_PID"].d, PIDs_gas_immediately_post_condensing)
        
        radii_pre_condensing = sp["PartType0", "particle_radius"][sp_pre_condensing_mask_for_analysis].to("kpc").value
        temperature_pre_condensing = sp["PartType0", "temperature"][sp_pre_condensing_mask_for_analysis].to("K").value
        
        
        radii_post_condensing = sp["PartType0", "particle_radius"][sp_post_condensing_mask_for_analysis].to("kpc").value
        temperature_post_condensing = sp["PartType0", "temperature"][sp_post_condensing_mask_for_analysis].to("K").value
        
        
        
        if PLOT_RADIUS_VS_TEMPERATURE_TESTING_PLOTS:
            fig = plt.figure(figsize = (6,6))
            plt.scatter(np.log10(temperature_post_condensing), np.log10(radii_post_condensing), color = "blue")
            plt.scatter(np.log10(temperature_pre_condensing), np.log10(radii_pre_condensing), color = "red")

            plt.xlabel("log temperature")
            plt.ylabel("log radius")
            plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/testing_PNGs_for_combined_data/TEST:zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_radius_vs_temp_of_condensing_parts.png", dpi=300)
            plt.close
            print(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/testing_PNGs_for_combined_data/TEST:zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_radius_vs_temp_of_condensing_parts.png")
            
            
            fig = plt.figure(figsize = (6,6))
            plt.scatter(np.log10(temperature_pre_condensing), np.log10(radii_pre_condensing), color = "red")
            # plt.scatter(np.log10(temperature_post_condensing), np.log10(radii_post_condensing), color = "blue")
            plt.xlabel("log temperature")
            plt.ylabel("log radius")
            plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/testing_PNGs_for_combined_data/TEST_hot_only:zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_radius_vs_temp_of_condensing_parts.png", dpi=300)
            plt.close
        
        
        
        '''####################################################'''
        '''                                                    '''
        '''                                                    '''
        '''Load All Cold Particle Data for this current snapnum'''
        '''   (not needed in this script, commented out)       '''
        '''                                                    '''
        '''####################################################'''
        # h5_filename = f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_DATA/" + "zm{}_Level{}_cold_data_sn{}.h5"
        # with h5py.File(h5_filename.format(zoomnum, zoomlevel, snapnum_this_rank), "r",) as f:
        #     cold_gas_IDs_this_snap                  =   list(f[f"COLD_GAS/{snapnum_this_rank}/COLD_IDs"] )
        #     cold_gas_unique_identifier_this_snap    =   list(f[f"COLD_GAS/{snapnum_this_rank}/COLD_Unique_Identifier"] )
        #     cold_gas_radii_this_snap                =   np.asarray(f[f"COLD_GAS/{snapnum_this_rank}/COLD_GAS_RADIUS"])
        #     cold_gas_temperature_this_snap          =   np.asarray(f[f"COLD_GAS/{snapnum_this_rank}/COLD_GAS_TEMP"])
            
        #     cold_gas_unique_PIDs_this_snap          =   ["{}_{}".format(a_, b_) for a_, b_ in zip(cold_gas_IDs_this_snap, cold_gas_unique_identifier_this_snap)]
            
            
        # # print("h5_filename, len(cold_gas_IDs_this_snap)",h5_filename, len(cold_gas_IDs_this_snap))
        # # exit()
        
        
        # ''' Make sure not particles are initially cold pre-condensing '''    
        # for bad_part_i, bad_part_PID in enumerate(sp["PartType0", "unique_PID"][sp_pre_condensing_mask_for_analysis].d[temperature_pre_condensing < 10**5]):
            
        #     if bad_part_PID in cold_gas_unique_PIDs_this_snap:
                
        #         h5_file_index = np.where(np.asarray(cold_gas_unique_PIDs_this_snap) == bad_part_PID)[0]
        #         this_script_bad_temp    = temperature_pre_condensing[temperature_pre_condensing > 10**5]
        #         bad_part_h5_temp        =  cold_gas_temperature_this_snap[h5_file_index]
                
        #         print("Bad particle IS in the h5 file pre-condensing!! Temperatures are", '%.2f' % float(np.log10(this_script_bad_temp)), '%.2f' % float(np.log10(bad_part_h5_temp)))
        
        #     else:
        #         print("Bad particle is not in the h5 file pre-condensing, which is WEIRD!!")
        
        # if len(sp["PartType0", "unique_PID"][sp_pre_condensing_mask_for_analysis].d[temperature_pre_condensing < 10**5]) > 0:
        #     raise RuntimeError()
        
        
        # cold_pids_this_snap = sp["PartType0", "unique_PID"].d[sp["PartType0", "temperature"].to("K").value < 10**5]
        
        # # print(h5_filename.format(zoomnum, zoomlevel, snapnum_this_rank), "\n", snapfile, " \n lens cold pids", len(cold_pids_this_snap) , len(cold_gas_unique_PIDs_this_snap), np.sum(np.in1d(np.asarray(cold_pids_this_snap), np.asarray(cold_gas_unique_PIDs_this_snap))) , "\n -------------- \n")
        # # print("cold_pids_this_snap[0], cold_gas_IDs_this_snap[0]", cold_pids_this_snap[0], cold_gas_unique_PIDs_this_snap[0] )
        # if not len(cold_pids_this_snap) == np.sum(np.in1d(np.asarray(cold_pids_this_snap), np.asarray(cold_gas_unique_PIDs_this_snap))):
        #     if PRINT_WARNING_MESSAGES: print("not len(cold_pids_this_snap) != np.sum(np.in1d(cold_pids_this_snap, cold_gas_IDs_this_snap)) !!")

    if RANK == 0: print(RANK, " Time to do condensepart loading block: ", time.time() - t_start_condense_block, flush = True)
    
    




    '''##############################'''
    '''                              '''
    '''                              '''
    '''Load All Bubble particles data'''
    '''                              '''
    '''                              '''
    '''##############################'''
    bubble_parts_dict = {}
    t_read_bubble_particles = time.time()
    print(RANK, "Reading Bubble Particle Data", flush = True)
    for bubble_snap in  np.arange( snapnum_this_rank - ANALYSE_BUBBLE_PARTS_FROM_PREVIOUS_N_SNAPS, snapnum_this_rank + 1, step = 1):
        
        # print("Considering bubble snap", bubble_snap)
        
        if bubble_snap > snapnum_this_rank: ### defend against bad coding skills
            continue
        
        if bubble_snap < RESTART_SNAPNUM:
            continue
        
        # print("Doing bubble snap", bubble_snap)
        
        
        if ALLOW_FILES_MISSING_FOR_TESTING and not os.path.exists(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_DATA/zm{zoomnum}_Level{zoomlevel}_cold_data_sn{bubble_snap}.h5"):
            if PRINT_WARNING_MESSAGES: print(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_DATA/zm{zoomnum}_Level{zoomlevel}_cold_data_sn{bubble_snap}.h5  does not exist but allowing this since we have ALLOW_FILES_MISSING_FOR_TESTING switched on...")
            continue
        
        with h5py.File(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_DATA/zm{zoomnum}_Level{zoomlevel}_cold_data_sn{bubble_snap}.h5",'r') as f:
            bubble_gas_IDs_this_snap                =   f[f"BUBBLE_PARTS/{bubble_snap}/IDs"] 
            bubble_gas_unique_identifier_this_snap  =   f[f"BUBBLE_PARTS/{bubble_snap}/Unique_Identifier"]
            # print("bubble_gas_IDs_this_snap", np.asarray(bubble_gas_IDs_this_snap))
            bubble_gas_unique_PIDs_this_snap        =   ["{}_{}".format(a_, b_) for a_, b_ in zip(bubble_gas_IDs_this_snap, bubble_gas_unique_identifier_this_snap)]
            bubble_parts_dict[bubble_snap]          = bubble_gas_unique_PIDs_this_snap
            
            
            
            
            
            
    
    print(RANK, "time to read bubble parts", time.time() - t_read_bubble_particles, flush = True)
    
    
    
    '''' NOW DO PLOTS WHERE WE MUST ALSO CONSIDER PREVIOUS SNAPSHOTS'''
    ''' Load the condensing data '''

    if DO_CONDENSED_PARTICLE_PLOTS or DO_BUBBLE_PARTICLE_PLOTS or DO_SATELLITE_PLOTS:
        
        
        t_start_projection_plot = time.time()
        
        
        os.makedirs(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/FITSFILES/", exist_ok=True)
        
        if PLOT_TYPE == "FILTERED_GAS_TEMPERATURE":
            
            fits_name = f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/FITSFILES/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_z.fits"
            
            
            buff_size = (800,800)
            raw_fields  = [("filtered_gas", "temperature"),]
            # parent_plot_width = [0.3*halo_r500.to("kpc"),0.3*halo_r500.to("kpc"), 0.3*halo_r500.to("kpc")]
            parent_plot_width = [2*NUM_R500_FOR_PROJ*halo_r500.to("kpc"),2*NUM_R500_FOR_PROJ*halo_r500.to("kpc"), 2*NUM_R500_FOR_PROJ*halo_r500.to("kpc")]
            parent_plot_extent = [-parent_plot_width[0].to('kpc').value/2,parent_plot_width[0].to('kpc').value/2,-parent_plot_width[1].to('kpc').value/2,parent_plot_width[1].to('kpc').value/2]
            
            if os.path.exists(fits_name) and not OVERWRITE_PARENT_PLOT:
                print(fits_name, " exists and OVERWRITE_PARENT_PLOT = ", OVERWRITE_PARENT_PLOT)
                pass
            else:
                prj_parent         = yt.ProjectionPlot(ds, normal = "z",  fields=raw_fields, weight_field=("filtered_gas", "density"),center = halo_center_post.to("kpc"), width = parent_plot_width, data_source = sp, buff_size = buff_size)
                fid_pw      = prj_parent.to_fits_data(fields=raw_fields, length_unit="kpc", image_res = buff_size[0], origin="image")
                fid_pw.writeto(fits_name, overwrite=True)   

            ''' reload data and check that the extent is correct'''
            hdulist = fits.open(fits_name)
            # print(repr(hdulist[0].header), flush = True)
            print(hdulist[0].header["NAXIS1"] * hdulist[0].header["CDELT1"], NUM_R500_FOR_PROJ*halo_r500.to('kpc').value, flush = True)
            dx = hdulist[0].header["CDELT1"]
            dy = hdulist[0].header["CDELT2"]
            half_width_x = hdulist[0].header["CDELT1"] * hdulist[0].header["NAXIS1"] / 2
            half_width_y = hdulist[0].header["CDELT2"] * hdulist[0].header["NAXIS2"] / 2
            extent = [-half_width_x, half_width_x, -half_width_y, half_width_y]
            # print(extent, parent_plot_extent)
            
            if not np.array_equal(np.around(np.asarray(extent), 2) , np.around(np.asarray(parent_plot_extent), 2)):
                raise RuntimeError()
            
            parent_cmap = "Greys_r"
            parent_img_data = np.asarray(hdulist[0].data)
            
         
        # if PLOT_TYPE == "TEMPERATURE":   
        #     if USE_EXISTING_FITSFILE_FOR_PLOTS:
                
        #         tmp_fits_name = f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_gas_proj_densweighted_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_kT.fits"
        #         if os.path.exists(tmp_fits_name):
        #             hdulist = fits.open(tmp_fits_name)
        #         else:
        #             fits_name = f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/FITSFILES/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_z.fits"
                    
                    
        #             buff_size = (800,800)
        #             raw_fields  = [("gas", "temperature"),]
        #             # parent_plot_width = [0.3*halo_r500.to("kpc"),0.3*halo_r500.to("kpc"), 0.3*halo_r500.to("kpc")]
        #             parent_plot_width = [2*NUM_R500_FOR_PROJ*halo_r500.to("kpc"),2*NUM_R500_FOR_PROJ*halo_r500.to("kpc"), 2*NUM_R500_FOR_PROJ*halo_r500.to("kpc")]
        #             parent_plot_extent = [-parent_plot_width[0].to('kpc').value/2,parent_plot_width[0].to('kpc').value/2,-parent_plot_width[1].to('kpc').value/2,parent_plot_width[1].to('kpc').value/2]
                    
        #             if os.path.exists(fits_name) and not OVERWRITE_PARENT_PLOT:
        #                 print(fits_name, " exists and OVERWRITE_PARENT_PLOT = ", OVERWRITE_PARENT_PLOT)
        #                 pass
        #             else:
        #                 prj_parent         = yt.ProjectionPlot(ds, normal = "z",  fields=raw_fields, weight_field=("gas", "density"),center = halo_center_post.to("kpc"), width = parent_plot_width, data_source = sp, buff_size = buff_size)
        #                 fid_pw      = prj_parent.to_fits_data(fields=raw_fields, length_unit="kpc", image_res = buff_size[0], origin="image")
        #                 fid_pw.writeto(fits_name, overwrite=True)   

        #             ''' reload data and check that the extent is correct'''
        #             hdulist = fits.open(fits_name)
                    
                    
                    
                    
                    
                    
            
            
        if PLOT_TYPE == "FILTERED_GAS_DENSITY":
            raise NotImplementedError()

    
        if DO_MULTIPLOT:
            print(RANK, "Initialising Multiplot", flush = True)
            for ax_i, axis in enumerate([multi_ax1, multi_ax2, multi_ax3, multi_ax4]):
                
                
                if len(MULTIPLOT_DICT[f"AXIS_{int(ax_i+1)}"]) == 1 and not PLOT_CONDENSATION_FROM_DIFFERENT_SOURCES_ON_AXIS_1:
                    im = axis.imshow(parent_img_data,  cmap = MULTI_PLOT_CMAP_IF_CLEAN, #cmap="cmr.apple",
                    extent = parent_plot_extent,
                    norm=LogNorm(vmin=PLOT_VMIN, vmax=PLOT_VMAX,)) 
                else:
                    im = axis.imshow(parent_img_data,  cmap = parent_cmap, #cmap="cmr.apple",
                    extent = parent_plot_extent,
                    norm=LogNorm(vmin=PLOT_VMIN, vmax=PLOT_VMAX,))
                    
                
                if ax_i <2: axis.set_xticklabels([])
                if ax_i == 1 or ax_i == 3: axis.set_yticklabels([])
                plt.figure(multi_fig)
                plt.subplots_adjust(wspace=0, hspace=0)
    
        if RANK == 0: print(RANK, "Time to do projection plot:", time.time()- t_start_projection_plot)
    if DO_SATELLITE_PLOTS:
        
        t_start_satellite_plots = time.time()
        
        fig, ax = plt.subplots(layout='constrained')
        fig, ax = plt.subplots()
        fig.set_size_inches((10,10))
        im = ax.imshow(parent_img_data,  cmap = parent_cmap, #cmap="cmr.apple",
                extent = parent_plot_extent,
                norm=LogNorm(vmin=PLOT_VMIN, vmax=PLOT_VMAX,))
         
        plt.scatter(x=wake_gas_x,y=wake_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "orange",  transform = ax.transData)
        plt.scatter(x=in_front_gas_x,y=in_front_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "blue",  transform = ax.transData)
        
        
        
        if DO_MULTIPLOT:
            if "SATELLITE_PARTS" in MULTIPLOT_DICT["AXIS_1"]:
                multi_ax1.scatter(x=in_front_gas_x,y=in_front_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "blue",  transform = multi_ax1.transData, zorder = 20_000)
                multi_ax1.scatter(x=wake_gas_x,y=wake_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "orange",  transform = multi_ax1.transData, zorder = 10_000)
            if "SATELLITE_PARTS" in MULTIPLOT_DICT["AXIS_2"]:
                multi_ax2.scatter(x=in_front_gas_x,y=in_front_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "blue",  transform = multi_ax2.transData, zorder = 20_000)
                multi_ax2.scatter(x=wake_gas_x,y=wake_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "orange",  transform = multi_ax2.transData, zorder = 10_000)
            if "SATELLITE_PARTS" in MULTIPLOT_DICT["AXIS_3"]:
                multi_ax3.scatter(x=in_front_gas_x,y=in_front_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "blue",  transform = multi_ax3.transData, zorder = 20_000)
                multi_ax3.scatter(x=wake_gas_x,y=wake_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "orange",  transform = multi_ax3.transData, zorder = 10_000)
            if "SATELLITE_PARTS" in MULTIPLOT_DICT["AXIS_4"]:
                multi_ax4.scatter(x=in_front_gas_x,y=in_front_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "blue",  transform = multi_ax4.transData, zorder = 20_000)
                multi_ax4.scatter(x=wake_gas_x,y=wake_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "orange",  transform = multi_ax4.transData, zorder = 10_000)    
            if "SATELLITE_PARTS" in MULTIPLOT_DICT["AXIS_1"]:
                multi_ax1.text(0.1,0.8,s = "Satellite-wake Particles", transform=multi_ax1.transAxes, color = onplot_text_color)
            if "SATELLITE_PARTS" in MULTIPLOT_DICT["AXIS_2"]:
                multi_ax2.text(0.1,0.8,s = "Satellite-wake Particles", transform=multi_ax2.transAxes, color = onplot_text_color)
            if "SATELLITE_PARTS" in MULTIPLOT_DICT["AXIS_3"]:
                multi_ax3.text(0.1,0.8,s = "Satellite-wake Particles", transform=multi_ax3.transAxes, color = onplot_text_color)
            if "SATELLITE_PARTS" in MULTIPLOT_DICT["AXIS_4"]:
                multi_ax4.text(0.1,0.8,s = "Satellite-wake Particles", transform=multi_ax4.transAxes, color = onplot_text_color)
        plt.figure(fig) 
        
        plt.text(0.1,0.8,s = "Satellite - wake", transform=ax.transAxes, color = onplot_text_color)
        plt.figure(fig)       
        plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_satellite_z.png", dpi=300)
        print(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_satellite_z.png", flush = True)
        plt.close()
        
        
        # exit()
        
        
        if RANK == 0: print(RANK, "Time to do satellite plots:", time.time() - t_start_satellite_plots)
    if DO_CONDENSED_PARTICLE_PLOTS:
        
        t_start_condense_plots = time.time()
        
        
        fig, ax = plt.subplots(layout='constrained')
        fig, ax = plt.subplots()
        fig.set_size_inches((10,10))
        im = ax.imshow(parent_img_data,  cmap = "Greys_r", # "cmr.lilac", #parent_cmap, #cmap="cmr.apple",
                extent = parent_plot_extent,
                norm=LogNorm(vmin=PLOT_VMIN, vmax=PLOT_VMAX,))
        ent_cax = ax.inset_axes([1.01, 0, 0.05, 1])
        cbar = fig.colorbar(im, label = "T (K)", cax = ent_cax)

        # Take colors at regular intervals spanning the colormap.
        # cmap = matplotlib.colormaps['cmr.cosmic_r']
        
        cmap = matplotlib.colormaps['cmr.neon_r']

        keep_condensed_parts_for_N_snaps_for_plot = 10000
        max_condensed_plot_lookback_color = 50
        # condensed_plot_colors = cmap(np.linspace(0, 1, ANALYSE_CONDENSED_PARTS_FROM_PREVIOUS_N_SNAPS + 1))
        condensed_plot_colors = cmap(np.linspace(0, 1, max_condensed_plot_lookback_color + 1))
        print("condensed_plot_colors", condensed_plot_colors)
        condensed_snaps_to_plot = []
        condensed_gas_PIDs_to_plot = []
        
        
        for condensed_i in range(len(condensed_gas_unique_PIDs)):
            if condensed_gas_snapshots_when_cold[condensed_i] <= snapnum_this_rank and condensed_gas_snapshots_when_cold[condensed_i] > snapnum_this_rank - ANALYSE_CONDENSED_PARTS_FROM_PREVIOUS_N_SNAPS:
                condensed_gas_PIDs_to_plot.append(condensed_gas_unique_PIDs[condensed_i])
                condensed_snaps_to_plot.append(condensed_gas_snapshots_when_cold[condensed_i])
            
        # sp_condensedpart_mask_for_plot = np.where(np.in1d(sp["PartType0", "unique_PID"].d, condensed_gas_PIDs_to_plot))[0]
        
        
        sp_condensedpart_mask_for_plot = np.in1d(sp["PartType0", "unique_PID"].d, condensed_gas_PIDs_to_plot)
        
        
        
        condensed_actually_unique_PIDs, condensed_unique_indexing,condensed_unique_inverse, condensed_counts = np.unique(sp["PartType0","unique_PID"].d,return_index =True, return_inverse=True,
                        return_counts=True)
        ''' Get an array same length as all the UniquePIDs with entries corresponding to the number of times that PID at that index appears'''
        condensed_num_counts_size_of_orig_arr = condensed_counts[condensed_unique_inverse]
        condensed_unique_PIDs_mask = np.asarray(condensed_num_counts_size_of_orig_arr) == 1
        sp_condensedpart_mask_for_plot &= condensed_unique_PIDs_mask
        
        
        condensed_gas_x_for_plot = sp["PartType0", "particle_position_x"][sp_condensedpart_mask_for_plot].to("kpc").value - halo_center_post[0].to("kpc").value
        condensed_gas_y_for_plot = sp["PartType0", "particle_position_y"][sp_condensedpart_mask_for_plot].to("kpc").value - halo_center_post[1].to("kpc").value
        condensed_gas_z_for_plot = sp["PartType0", "particle_position_z"][sp_condensedpart_mask_for_plot].to("kpc").value - halo_center_post[2].to("kpc").value



        if len(condensed_gas_x_for_plot) > 0:
            for i_condense in range(len(condensed_gas_x_for_plot)):
                # print("snapnum_this_rank - condensed_snaps_to_plot[i_condense]", snapnum_this_rank - condensed_snaps_to_plot[i_condense], len(condensed_plot_colors) )
                
                if snapnum_this_rank - condensed_snaps_to_plot[i_condense] > keep_condensed_parts_for_N_snaps_for_plot:
                    continue
                
                ''' Must do 1-by-1 otherwise the following sp-derived lists will be mixed up compared to condensed_snaps_to_plot'''
                sp_condensedpart_mask_for_plot = np.in1d(sp["PartType0", "unique_PID"].d, condensed_gas_PIDs_to_plot[i_condense])
                condensed_gas_x_for_plot = sp["PartType0", "particle_position_x"][sp_condensedpart_mask_for_plot].to("kpc").value - halo_center_post[0].to("kpc").value
                condensed_gas_y_for_plot = sp["PartType0", "particle_position_y"][sp_condensedpart_mask_for_plot].to("kpc").value - halo_center_post[1].to("kpc").value
                condensed_gas_z_for_plot = sp["PartType0", "particle_position_z"][sp_condensedpart_mask_for_plot].to("kpc").value - halo_center_post[2].to("kpc").value
                
                if len(condensed_gas_x_for_plot) > 1:
                    print("WARNING!! DUPLICATE CONDENSED PIDs PRESENT!")
                    continue
                    
                if condensed_gas_x_for_plot > parent_plot_extent[0] and condensed_gas_x_for_plot < parent_plot_extent[1] and condensed_gas_y_for_plot > parent_plot_extent[2] and condensed_gas_y_for_plot < parent_plot_extent[3]:
                    
                    # print("snapnum_this_rank - condensed_snaps_to_plot[i_condense]", snapnum_this_rank - condensed_snaps_to_plot[i_condense])
                    

                    tmp_color = condensed_plot_colors[min(snapnum_this_rank - condensed_snaps_to_plot[i_condense], max_condensed_plot_lookback_color)]
                    bubble_plot_alpha_fac = 15/(2*min(snapnum_this_rank - condensed_snaps_to_plot[i_condense], max_condensed_plot_lookback_color) + 10)
                    bubble_plot_s_fac = 30/(min(snapnum_this_rank - condensed_snaps_to_plot[i_condense], max_condensed_plot_lookback_color) + 10)
                    condense_zorder = 1/min(snapnum_this_rank - condensed_snaps_to_plot[i_condense], max_condensed_plot_lookback_color) 
                    
                    plt.scatter(x=condensed_gas_x_for_plot,y=condensed_gas_y_for_plot, s=bubble_plot_s_fac*PLOT_S,  alpha = bubble_plot_alpha_fac*PLOT_ALPHA, color = tmp_color, zorder = condense_zorder,  transform = ax.transData)
                    if DO_MULTIPLOT:
                        if "CONDENSED_PARTS" in MULTIPLOT_DICT["AXIS_1"]:
                            multi_ax1.scatter(x=condensed_gas_x_for_plot,y=condensed_gas_y_for_plot, s=bubble_plot_s_fac*PLOT_S,  alpha = bubble_plot_alpha_fac*PLOT_ALPHA, color = tmp_color, zorder = condense_zorder,  transform = multi_ax1.transData)
                        if "CONDENSED_PARTS" in MULTIPLOT_DICT["AXIS_2"]:
                            multi_ax2.scatter(x=condensed_gas_x_for_plot,y=condensed_gas_y_for_plot, s=bubble_plot_s_fac*PLOT_S,  alpha = bubble_plot_alpha_fac*PLOT_ALPHA, color = tmp_color, zorder = condense_zorder,  transform = multi_ax2.transData)
                        if "CONDENSED_PARTS" in MULTIPLOT_DICT["AXIS_3"]:
                            multi_ax3.scatter(x=condensed_gas_x_for_plot,y=condensed_gas_y_for_plot, s=bubble_plot_s_fac*PLOT_S,  alpha = bubble_plot_alpha_fac*PLOT_ALPHA, color = tmp_color, zorder = condense_zorder,  transform = multi_ax3.transData)
                        if "CONDENSED_PARTS" in MULTIPLOT_DICT["AXIS_4"]:
                            multi_ax4.scatter(x=condensed_gas_x_for_plot,y=condensed_gas_y_for_plot, s=bubble_plot_s_fac*PLOT_S,  alpha = bubble_plot_alpha_fac*PLOT_ALPHA, color = tmp_color, zorder = condense_zorder,  transform = multi_ax4.transData)
                            
                    plt.figure(fig)     
                    
                    
                       
        
        if DO_MULTIPLOT:
            if "CONDENSED_PARTS" in MULTIPLOT_DICT["AXIS_1"]:
                multi_ax1.text(0.1,0.8,s = "Condensed Particles", transform=multi_ax1.transAxes, color = onplot_text_color)
            if "CONDENSED_PARTS" in MULTIPLOT_DICT["AXIS_2"]:
                multi_ax2.text(0.1,0.8,s = "Condensed Particles", transform=multi_ax2.transAxes, color = onplot_text_color)
            if "CONDENSED_PARTS" in MULTIPLOT_DICT["AXIS_3"]:
                multi_ax3.text(0.1,0.8,s = "Condensed Particles", transform=multi_ax3.transAxes, color = onplot_text_color)
            if "CONDENSED_PARTS" in MULTIPLOT_DICT["AXIS_4"]:
                multi_ax4.text(0.1,0.8,s = "Condensed Particles", transform=multi_ax4.transAxes, color = onplot_text_color)
        plt.figure(fig) 
        
        
        
        
        # plt.text(0.1,0.8,s = "Condensed Particles", transform=ax.transAxes, color = onplot_text_color)
        
        
        plt.text(0.07,0.07,s = f"{'%.2f' % float(ds.current_time.to('Myr').value - start_time)} Myr", transform=ax.transAxes, color = "white", size = 30, bbox=dict(boxstyle="round", facecolor='black', edgecolor='yellow'))
        
        plt.figure(fig)     
        # try:
        #     xticklocs = ax.get_xticks() 
        #     xtick_labels = ax.get_xticklabels()
        #     ax2 = ax.twiny()
        #     ax2.set_xlim(ax.get_xlim())
        #     ax2.set_xticks(xticklocs)
        #     ax2.set_xticklabels(['%.2f' % float(float(x)/halo_r500) for x in  xtick_labels])
        #     ax2.set_xlabel("R/R500")
        # except Exception as e:
        #     print(e)
            
        ax.set_xlabel("x (kpc)")
        ax.set_ylabel("y (kpc)")
        plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_with_condenseparts_z.png", dpi=300)
        print(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_with_condenseparts_z.png")
        plt.close()


        print(RANK, "time to do condensed particle plots:", time.time() - t_start_condense_plots, flush = True)
    if DO_BUBBLE_PARTICLE_PLOTS:
        print(RANK, "Doing Bubble Particle Plots", flush = True)
        
        cmap = matplotlib.colormaps['cmr.neon']
        max_bubble_plot_lookback_color = 30
        # bubble_plot_colors = cmap(np.linspace(0, 1, ANALYSE_BUBBLE_PARTS_FROM_PREVIOUS_N_SNAPS + 1))
        bubble_plot_colors = cmap(np.linspace(0, 1, max_bubble_plot_lookback_color + 1))
        
        t_bubble_part_plots = time.time()
        
        fig, ax = plt.subplots(layout='constrained')
        fig, ax = plt.subplots()
        fig.set_size_inches((10,10))
        im = ax.imshow(parent_img_data,  cmap = parent_cmap, #cmap="cmr.apple",
                extent = parent_plot_extent,
                norm=LogNorm(vmin=PLOT_VMIN, vmax=PLOT_VMAX,))
        
        ent_cax = ax.inset_axes([1.01, 0, 0.05, 1])
        cbar = fig.colorbar(im, label = "T (K)", cax = ent_cax)
        
        

        plt.figure(fig) 
        plt.figure(fig)      
        plt.text(0.07,0.07,s = f"{'%.2f' % float(ds.current_time.to('Myr').value - start_time)} Myr", transform=ax.transAxes, color = "white", size = 30, bbox=dict(boxstyle="round", facecolor='black', edgecolor='yellow'))
        
        ax.set_xlabel("x (kpc)")
        ax.set_ylabel("y (kpc)")
        plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_with_bubble_parts_z.png", dpi=300)
        print(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_with_bubble_parts_z.png")
        plt.close()
            
            
        # if DO_BUBBLE_PARTICLE_PLOTS:
        #     plt.text(0.1,0.8,s = "Bubble Particles", transform=ax.transAxes, color = onplot_text_color)
        if DO_MULTIPLOT:
            if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_1"]:
                multi_ax1.text(0.1,0.8,s = "Bubble Particles", transform=multi_ax1.transAxes, color = onplot_text_color)
            if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_2"]:
                multi_ax2.text(0.1,0.8,s = "Bubble Particles", transform=multi_ax2.transAxes, color = onplot_text_color)
            if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_3"]:
                multi_ax3.text(0.1,0.8,s = "Bubble Particles", transform=multi_ax3.transAxes, color = onplot_text_color)
            if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_4"]:
                multi_ax4.text(0.1,0.8,s = "Bubble Particles", transform=multi_ax4.transAxes, color = onplot_text_color)
        





    '''#########################################################'''
    '''                                                         '''
    '''                                                         '''
    '''  Now we must see which bubble particles have cooled by  '''
    '''the current snapnum. Also record still-hot particles and '''
    '''           bubble particles launched this snap           '''
    '''                                                         '''
    '''                                                         '''
    '''#########################################################'''

    if DO_BUBBLE_PARTICLE_PLOTS or DO_ANALYSIS:
        t_bubble_part_analysis = time.time()
        
        recent_still_hot_bubble_particle_masses_arr          = []
        recent_still_hot_bubble_particle_radii_arr           = []
        recent_still_hot_bubble_particle_x_arr           = []
        recent_still_hot_bubble_particle_y_arr           = []
        recent_still_hot_bubble_particle_z_arr           = []
        recent_still_hot_bubble_particle_ages_arr            = []
        recent_still_hot_bubble_particle_temperature_arr     = []
        
        recent_cooled_bubble_unique_PIDs_arr = []
        recent_cooled_bubble_particle_masses_arr             = []
        recent_cooled_bubble_particle_radii_arr              = []
        recent_cooled_bubble_particle_x_arr           = []
        recent_cooled_bubble_particle_y_arr           = []
        recent_cooled_bubble_particle_z_arr           = []
        recent_cooled_bubble_particle_ages_arr               = []
        recent_cooled_bubble_particle_temperature_arr        = []
        launched_this_snap_bubble_gas_radius                 = None
        launched_this_snap_bubble_gas_mass                   = None
        launched_this_snap_bubble_gas_temperature            = None
        

        for bubble_i, bubble_snap in enumerate(np.arange( snapnum_this_rank - ANALYSE_BUBBLE_PARTS_FROM_PREVIOUS_N_SNAPS, snapnum_this_rank + 1, step = 1)):
            
            if bubble_snap not in list(bubble_parts_dict.keys()):
                # print(f"bubble_snap {bubble_snap} not in list(bubble_parts_dict.keys())    -- maybe no bubble particles this snap")
                continue
            this_bubble_unique_PIDS = bubble_parts_dict[bubble_snap]
            # sp_bubblepart_mask = np.where(np.in1d(sp["PartType0", "unique_PID"].d, this_bubble_unique_PIDS))[0]
            sp_hot_bubblepart_mask = np.in1d(sp["PartType0", "unique_PID"].d, this_bubble_unique_PIDS)&(sp["PartType0", "temperature"].to("K").value > MIN_BUBBLE_T_K)
            hot_bubble_gas_x = sp["PartType0", "particle_position_x"][sp_hot_bubblepart_mask].to("kpc").value - halo_center_post[0].to("kpc").value
            hot_bubble_gas_y = sp["PartType0", "particle_position_y"][sp_hot_bubblepart_mask].to("kpc").value - halo_center_post[1].to("kpc").value
            hot_bubble_gas_z = sp["PartType0", "particle_position_z"][sp_hot_bubblepart_mask].to("kpc").value - halo_center_post[2].to("kpc").value
            hot_bubble_gas_radius = sp["PartType0", "particle_radius"][sp_hot_bubblepart_mask].to("kpc").value
            hot_bubble_gas_mass = sp["PartType0", "particle_mass"][sp_hot_bubblepart_mask].to("Msun").value
            hot_bubble_gas_temperature = sp["PartType0", "temperature"][sp_hot_bubblepart_mask].to("K").value
            hot_bubble_gas_ages = np.ones_like(hot_bubble_gas_mass) * int(snapnum_this_rank-bubble_snap) * DT
            
            sp_cooled_bubblepart_mask = np.in1d(sp["PartType0", "unique_PID"].d, this_bubble_unique_PIDS)&(sp["PartType0", "temperature"].to("K").value <= MIN_BUBBLE_T_K)
            cooled_bubble_gas_unique_PID = sp["PartType0", "unique_PID"][sp_cooled_bubblepart_mask].value
            cooled_bubble_gas_x = sp["PartType0", "particle_position_x"][sp_cooled_bubblepart_mask].to("kpc").value - halo_center_post[0].to("kpc").value
            cooled_bubble_gas_y = sp["PartType0", "particle_position_y"][sp_cooled_bubblepart_mask].to("kpc").value - halo_center_post[1].to("kpc").value
            cooled_bubble_gas_z = sp["PartType0", "particle_position_z"][sp_cooled_bubblepart_mask].to("kpc").value - halo_center_post[2].to("kpc").value
            cooled_bubble_gas_radius = sp["PartType0", "particle_radius"][sp_cooled_bubblepart_mask].to("kpc").value
            cooled_bubble_gas_mass = sp["PartType0", "particle_mass"][sp_cooled_bubblepart_mask].to("Msun").value
            cooled_bubble_gas_temperature = sp["PartType0", "temperature"][sp_cooled_bubblepart_mask].to("K").value
            cooled_bubble_gas_ages = np.ones_like(cooled_bubble_gas_mass) * int(snapnum_this_rank-bubble_snap) * DT

            recent_still_hot_bubble_particle_masses_arr.extend(list(hot_bubble_gas_mass))   
            recent_still_hot_bubble_particle_radii_arr.extend(list(hot_bubble_gas_radius))    
            recent_still_hot_bubble_particle_ages_arr.extend(list(hot_bubble_gas_ages))    
            recent_still_hot_bubble_particle_temperature_arr.extend(list(hot_bubble_gas_temperature)) 
            recent_still_hot_bubble_particle_x_arr.extend(list(hot_bubble_gas_x))
            recent_still_hot_bubble_particle_y_arr.extend(list(hot_bubble_gas_y))
            recent_still_hot_bubble_particle_z_arr.extend(list(hot_bubble_gas_z))
            recent_cooled_bubble_unique_PIDs_arr.extend(list(cooled_bubble_gas_unique_PID))
            recent_cooled_bubble_particle_masses_arr.extend(list(cooled_bubble_gas_mass))    
            recent_cooled_bubble_particle_radii_arr.extend(list(cooled_bubble_gas_radius))  
            recent_cooled_bubble_particle_x_arr.extend(list(cooled_bubble_gas_x))
            recent_cooled_bubble_particle_y_arr.extend(list(cooled_bubble_gas_y))
            recent_cooled_bubble_particle_z_arr.extend(list(cooled_bubble_gas_z))  
            recent_cooled_bubble_particle_ages_arr.extend(list(cooled_bubble_gas_ages))    
            recent_cooled_bubble_particle_temperature_arr.extend(list(cooled_bubble_gas_temperature))

            if bubble_snap == snapnum_this_rank:
                launched_this_snap_bubble_gas_radius        = sp["PartType0", "particle_radius"][sp_hot_bubblepart_mask].to("kpc").value
                launched_this_snap_bubble_gas_mass          = sp["PartType0", "particle_mass"][sp_hot_bubblepart_mask].to("Msun").value
                launched_this_snap_bubble_gas_temperature   = sp["PartType0", "temperature"][sp_hot_bubblepart_mask].to("K").value    
  
            if not np.allclose(np.asarray(hot_bubble_gas_radius), np.sqrt(hot_bubble_gas_x**2 + hot_bubble_gas_y**2 + hot_bubble_gas_z**2)): ### Make sure yT is doing all the centering stuff properly for the radius field

                print(np.shape(hot_bubble_gas_radius), np.shape(np.sqrt(hot_bubble_gas_x**2 + hot_bubble_gas_y**2 + hot_bubble_gas_z**2)))
                print(hot_bubble_gas_radius - np.sqrt(hot_bubble_gas_x**2 + hot_bubble_gas_y**2 + hot_bubble_gas_z**2))
                print(hot_bubble_gas_radius , np.sqrt(hot_bubble_gas_x**2 + hot_bubble_gas_y**2 + hot_bubble_gas_z**2))
                raise ValueError()

            if DO_BUBBLE_PARTICLE_PLOTS:
                # print("bubble_i", bubble_i)
                # bubble_color = bubble_plot_colors[bubble_i]
                bubble_color = bubble_plot_colors[min(snapnum_this_rank-bubble_snap, max_bubble_plot_lookback_color)]
                
                plt.scatter(x=cooled_bubble_gas_x,y=cooled_bubble_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "grey",  transform = ax.transData)
                plt.scatter(x=hot_bubble_gas_x,y=hot_bubble_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = bubble_color,  transform = ax.transData)
                
                if DO_MULTIPLOT:
                    
                    if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_1"]:
                        multi_ax1.scatter(x=cooled_bubble_gas_x,y=cooled_bubble_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "grey",  transform = multi_ax1.transData)
                    if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_2"]:
                        multi_ax2.scatter(x=cooled_bubble_gas_x,y=cooled_bubble_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "grey",  transform = multi_ax2.transData)
                    if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_3"]:
                        multi_ax3.scatter(x=cooled_bubble_gas_x,y=cooled_bubble_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "grey",  transform = multi_ax3.transData)
                    if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_4"]:
                        multi_ax4.scatter(x=cooled_bubble_gas_x,y=cooled_bubble_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = "grey",  transform = multi_ax4.transData)
                        
                                
                    if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_1"]:
                        multi_ax1.scatter(x=hot_bubble_gas_x,y=hot_bubble_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = bubble_color,  transform = multi_ax1.transData)
                    if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_2"]:
                        multi_ax2.scatter(x=hot_bubble_gas_x,y=hot_bubble_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = bubble_color,  transform = multi_ax2.transData)
                    if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_3"]:
                        multi_ax3.scatter(x=hot_bubble_gas_x,y=hot_bubble_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = bubble_color,  transform = multi_ax3.transData)
                    if "BUBBLE_PARTS" in MULTIPLOT_DICT["AXIS_4"]:
                        multi_ax4.scatter(x=hot_bubble_gas_x,y=hot_bubble_gas_y, s=PLOT_S,  alpha = PLOT_ALPHA, color = bubble_color,  transform = multi_ax4.transData)
                              
                plt.figure(fig) 
            

            # if RANK == 0: print(RANK, "Time to do bubble particle plots: ", time.time() - t_bubble_part_plots)
            if RANK == 0: print(RANK, "Time to do bubble part analysis", time.time() - t_bubble_part_analysis, flush = True)

    if DO_MULTIPLOT:
        plt.figure(multi_fig)   
        
        os.makedirs(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/MULTIFIGS/",    exist_ok=True)    
        plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/MULTIFIGS/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_multifig_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_z.png", dpi=300)
        print(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/MULTIFIGS/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_multifig_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_z.png")
        plt.close()
    
    
    
    
    '''##############################################'''
    '''                                              '''
    '''                                              '''
    '''Build trees for satellite and bubble particles'''
    ''' and find particles which condense near them  '''
    '''                                              '''
    '''                                              '''
    '''##############################################'''
    if RANK == 0: print(RANK, "Now doing the Analysis Part", flush=True)
    
    
    if DO_ANALYSIS:
    
        instantaneously_condensed_gas_PIDs_for_analysis = []
        
        for condensed_i in range(len(condensed_gas_unique_PIDs)):
            if condensed_gas_snapshots_when_cold[condensed_i] == snapnum_this_rank:
                instantaneously_condensed_gas_PIDs_for_analysis.append(condensed_gas_unique_PIDs[condensed_i])

            
        sp_condensedpart_mask_for_analysis = np.in1d(sp["PartType0", "unique_PID"].d, instantaneously_condensed_gas_PIDs_for_analysis)
        instantaneously_condensed_gas_PIDs_for_analysis = None   #### the order of this array is different to the order of those derived from tghe above mask. kill this variable here for safety
        instantaneously_condensed_gas_unique_PID_for_analysis = sp["PartType0", "unique_PID"][sp_condensedpart_mask_for_analysis].value
        instantaneously_condensed_gas_x_for_analysis = sp["PartType0", "particle_position_x"][sp_condensedpart_mask_for_analysis].to("kpc").value - halo_center_post[0].to("kpc").value
        instantaneously_condensed_gas_y_for_analysis = sp["PartType0", "particle_position_y"][sp_condensedpart_mask_for_analysis].to("kpc").value - halo_center_post[1].to("kpc").value
        instantaneously_condensed_gas_z_for_analysis = sp["PartType0", "particle_position_z"][sp_condensedpart_mask_for_analysis].to("kpc").value - halo_center_post[2].to("kpc").value
        instantaneously_condensed_gas_mass_for_analysis = sp["PartType0", "particle_mass"][sp_condensedpart_mask_for_analysis].to("Msun").value
        instantaneously_condensed_gas_radii_for_analysis = sp["PartType0", "particle_radius"][sp_condensedpart_mask_for_analysis].to("kpc").value
        
        
        
        instantaneously_condensed_gas_tree  = cKDTree(np.vstack((instantaneously_condensed_gas_x_for_analysis,instantaneously_condensed_gas_y_for_analysis,instantaneously_condensed_gas_z_for_analysis)).T)
        recent_still_hot_bubble_tree        = cKDTree(np.vstack((recent_still_hot_bubble_particle_x_arr,recent_still_hot_bubble_particle_y_arr,recent_still_hot_bubble_particle_z_arr)).T)
        satellite_wake_tree                 = cKDTree(np.vstack((wake_gas_x,wake_gas_y,wake_gas_z)).T)
        satellite_infront_tree              = cKDTree(np.vstack((in_front_gas_x,in_front_gas_y,in_front_gas_z)).T)

        
        
        ''' Nearness of condensation to bubbles'''
        close_to_bubbleparts_indexes = instantaneously_condensed_gas_tree.query_ball_tree(recent_still_hot_bubble_tree, r = kpc_distance_threshold_for_bubble_to_condense_parts)
        close_to_bubbleparts_MASK = np.asarray([len(x)>0 for x in close_to_bubbleparts_indexes])
        close_to_bubbleparts_MASK = (close_to_bubbleparts_MASK==1)
        # print("close_to_bubbleparts_MASK", close_to_bubbleparts_MASK)
        # print("close_to_bubbleparts_MASK",close_to_bubbleparts_MASK)
        assert len(close_to_bubbleparts_MASK) == len(instantaneously_condensed_gas_x_for_analysis)
        close_to_bubbleparts_PIDs = instantaneously_condensed_gas_unique_PID_for_analysis[close_to_bubbleparts_MASK]
        close_to_bubbleparts_masses = instantaneously_condensed_gas_mass_for_analysis[close_to_bubbleparts_MASK]
        close_to_bubbleparts_radii = instantaneously_condensed_gas_radii_for_analysis[close_to_bubbleparts_MASK]
        
        assert(len(set(close_to_bubbleparts_PIDs)) == len(close_to_bubbleparts_PIDs))   ### Make sure all unique PIDs and the mask has operated correctly and not indexed instead of masked
        
        ''' Condensation from cooled bubble material itself'''
        condensed_bubble_material_MASK = np.in1d(np.asarray(instantaneously_condensed_gas_unique_PID_for_analysis), np.asarray(recent_cooled_bubble_unique_PIDs_arr))
        assert len(condensed_bubble_material_MASK) == len(instantaneously_condensed_gas_x_for_analysis)
        condensed_bubble_material_PIDs = instantaneously_condensed_gas_unique_PID_for_analysis[condensed_bubble_material_MASK]
        condensed_bubble_material_masses = instantaneously_condensed_gas_mass_for_analysis[condensed_bubble_material_MASK]
        condensed_bubble_material_radii = instantaneously_condensed_gas_radii_for_analysis[condensed_bubble_material_MASK]
        assert(len(set(condensed_bubble_material_PIDs)) == len(condensed_bubble_material_PIDs))   ### Make sure all unique PIDs and the mask has operated correctly and not indexed instead of masked
        
        
        ''' Condensation near/in satellite wake'''
        close_to_satellite_wake_indexes = instantaneously_condensed_gas_tree.query_ball_tree(satellite_wake_tree, r = kpc_distance_threshold_for_condensed_to_wake)
        close_to_satellite_wake_MASK = np.asarray([len(x)>0 for x in close_to_satellite_wake_indexes])   
        close_to_satellite_wake_MASK = (close_to_satellite_wake_MASK == True)|(np.in1d(np.asarray(instantaneously_condensed_gas_unique_PID_for_analysis), np.asarray(wake_PIDS)))
        assert len(close_to_satellite_wake_MASK) == len(instantaneously_condensed_gas_x_for_analysis)
        close_to_satellite_wake_PIDs = instantaneously_condensed_gas_unique_PID_for_analysis[close_to_satellite_wake_MASK]
        close_to_satellite_wake_masses = instantaneously_condensed_gas_mass_for_analysis[close_to_satellite_wake_MASK]
        close_to_satellite_wake_radii = instantaneously_condensed_gas_radii_for_analysis[close_to_satellite_wake_MASK]
        assert(len(set(close_to_satellite_wake_PIDs)) == len(close_to_satellite_wake_PIDs))   ### Make sure all unique PIDs and the mask has operated correctly and not indexed instead of masked
        
        
        ''' Condensation in front of satellite'''
        close_to_satellite_infront_indexes = instantaneously_condensed_gas_tree.query_ball_tree(satellite_infront_tree, r = kpc_distance_threshold_for_condensed_to_satellite_infront)
        close_to_satellite_infront_MASK = np.asarray([len(x)>0 for x in close_to_satellite_infront_indexes])
        close_to_satellite_infront_MASK = (close_to_satellite_infront_MASK == True)|(np.in1d(np.asarray(instantaneously_condensed_gas_unique_PID_for_analysis), np.asarray(infront_of_satellite_PIDs)))
        assert len(close_to_satellite_infront_MASK) == len(instantaneously_condensed_gas_x_for_analysis)
        close_to_satellite_infront_PIDs = instantaneously_condensed_gas_unique_PID_for_analysis[close_to_satellite_infront_MASK]
        close_to_satellite_infront_masses = instantaneously_condensed_gas_mass_for_analysis[close_to_satellite_infront_MASK]
        close_to_satellite_infront_radii = instantaneously_condensed_gas_radii_for_analysis[close_to_satellite_infront_MASK]

        assert(len(set(close_to_satellite_infront_PIDs)) == len(close_to_satellite_infront_PIDs))   ### Make sure all unique PIDs and the mask has operated correctly and not indexed instead of masked
        

                            
        if DO_MULTIPLOT and PLOT_CONDENSATION_FROM_DIFFERENT_SOURCES_ON_AXIS_1:
            multi_ax1.scatter(instantaneously_condensed_gas_x_for_analysis[close_to_bubbleparts_MASK] ,instantaneously_condensed_gas_y_for_analysis[close_to_bubbleparts_MASK],label = "From Cavities", alpha = PLOT_ALPHA, color = "red",  transform = multi_ax1.transData)
            multi_ax1.scatter(instantaneously_condensed_gas_x_for_analysis[close_to_satellite_wake_MASK] ,instantaneously_condensed_gas_y_for_analysis[close_to_satellite_wake_MASK],label = "From Wake", alpha = PLOT_ALPHA, color = "orange",  transform = multi_ax1.transData)
            multi_ax1.scatter(instantaneously_condensed_gas_x_for_analysis[close_to_satellite_infront_MASK] ,instantaneously_condensed_gas_y_for_analysis[close_to_satellite_infront_MASK],label = "From sat in-front", alpha = PLOT_ALPHA, color = "blue",  transform = multi_ax1.transData)
            # print(snapnum_this_rank, "PLOT_CONDENSATION_FROM_DIFFERENT_SOURCES_ON_AXIS_1 instantaneously_condensed_gas_x_for_analysis[close_to_bubbleparts_MASK]", instantaneously_condensed_gas_x_for_analysis[close_to_bubbleparts_MASK])
            # print(snapnum_this_rank, "PLOT_CONDENSATION_FROM_DIFFERENT_SOURCES_ON_AXIS_1 instantaneously_condensed_gas_x_for_analysis[close_to_satellite_wake_MASK]", instantaneously_condensed_gas_x_for_analysis[close_to_satellite_wake_MASK])

            multi_ax1.legend()
        
        
        if PLOT_TO_SHOW_CONDENSING_CLOSE_TO_BUBBLES_AND_WAKE:
            print("Plotting show_condensing_near_bubbles_and_wake")
            fig, ax = plt.subplots(layout='constrained')
            fig, ax = plt.subplots()
            fig.set_size_inches((10,10))

            buff_size_tmp = (800,800)
            raw_fields_tmp  = [("filtered_gas", "temperature"),]
            # parent_plot_width = [0.3*halo_r500.to("kpc"),0.3*halo_r500.to("kpc"), 0.3*halo_r500.to("kpc")]
            tmp_plot_width = [2*NUM_R500_FOR_PROJ*halo_r500.to("kpc"),2*NUM_R500_FOR_PROJ*halo_r500.to("kpc"), 2*NUM_R500_FOR_PROJ*halo_r500.to("kpc")]
            tmp_plot_extent = [-tmp_plot_width[0].to('kpc').value/2,tmp_plot_width[0].to('kpc').value/2,-tmp_plot_width[1].to('kpc').value/2,tmp_plot_width[1].to('kpc').value/2]
            prj_tmp         = yt.ProjectionPlot(ds, normal = "z",  fields=raw_fields_tmp, weight_field=("filtered_gas", "density"),center = halo_center_post.to("kpc"), width = tmp_plot_width, data_source = sp, buff_size = buff_size_tmp)
            # fid_tmp      = prj_parent.to_fits_data(fields=raw_fields_tmp, length_unit="kpc", image_res = buff_size_tmp[0], origin="image")
            im = ax.imshow(prj_tmp.frb[raw_fields_tmp[0]],  cmap = "Greys_r", #cmap="cmr.apple",
                    extent = tmp_plot_extent,
                    norm=LogNorm(vmin=PLOT_VMIN, vmax=PLOT_VMAX,))
            
            
            
            ''' These can be overlaid on each other'''
            ax.scatter(instantaneously_condensed_gas_x_for_analysis[close_to_satellite_wake_MASK] ,instantaneously_condensed_gas_y_for_analysis[close_to_satellite_wake_MASK], s = PLOT_S, label = "From Wake", alpha = 1, color = "green",  transform = ax.transData)
            ax.scatter(instantaneously_condensed_gas_x_for_analysis[close_to_satellite_infront_MASK] ,instantaneously_condensed_gas_y_for_analysis[close_to_satellite_infront_MASK], s = 0.7*PLOT_S, label = "From sat in-front", alpha = 1, color = "blue",  transform = ax.transData)
            ax.scatter(instantaneously_condensed_gas_x_for_analysis[close_to_bubbleparts_MASK] ,instantaneously_condensed_gas_y_for_analysis[close_to_bubbleparts_MASK], s = 0.5*PLOT_S, label = "From Cavities", alpha = 1, color = "red",  transform = ax.transData)
            ax.legend()
            
            

            plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/show_condensing_near_bubbles_and_wake_zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_z.png", dpi=300)
            print(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_AND_BUBBLES_COMBINED_DATA/PNGs_for_combined_data/show_condensing_near_bubbles_and_wake_zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_{PLOT_TYPE}_projection_no_wcs_{'%.2f' % float(NUM_R500_FOR_PROJ*halo_r500.to('kpc').value)}kpc_z.png")
            plt.close()
        
        # exit()
        
        
        
        
        
        
        '''#######################################################'''
        '''                                                       '''
        '''                                                       '''
        '''      Find which particles have condensed within       '''
        '''the past ANALYSE_CONDENSED_PARTS_FROM_PREVIOUS_N_SNAPS '''
        '''                  to include for the                   '''
        '''       analysis of recently-condensed particles.       '''
        '''                                                       '''
        '''                                                       '''
        '''#######################################################'''
        recently_condensed_gas_PIDs_for_analysis = []
        for condensed_i in range(len(condensed_gas_unique_PIDs)):
            if condensed_gas_snapshots_when_cold[condensed_i] <= snapnum_this_rank and condensed_gas_snapshots_when_cold[condensed_i] > snapnum_this_rank - ANALYSE_CONDENSED_PARTS_FROM_PREVIOUS_N_SNAPS:
                recently_condensed_gas_PIDs_for_analysis.append(condensed_gas_unique_PIDs[condensed_i])
                
        ''' Mask into current datasphere of recently-condensed particles that exist in the current datasphere'''
        sp_condensedpart_mask_for_analysis = np.in1d(sp["PartType0", "unique_PID"].d, recently_condensed_gas_PIDs_for_analysis)
        
        if len(sp["PartType0", "unique_PID"][sp_condensedpart_mask_for_analysis].value) != recently_condensed_gas_PIDs_for_analysis:
                if PRINT_WARNING_MESSAGES :print("Warning! len(recently_condensed_gas_PIDs_for_analysis_found_in_sphere without temperature mask) != recently_condensed_gas_PIDs_for_analysis did you use different sphere sizes between scripts?")
                if not ALLOW_FILES_MISSING_FOR_TESTING:
                    raise RuntimeError()
                
                
        ''' ONLY include particles which are still cold in this current snapshot'''
        sp_condensedpart_mask_for_analysis &= (sp["PartType0", "temperature"].to("K").value < 1e5)
        recently_condensed_gas_x_for_analysis       = sp["PartType0", "particle_position_x"][sp_condensedpart_mask_for_analysis].to("kpc").value - halo_center_post[0].to("kpc").value
        recently_condensed_gas_y_for_analysis       = sp["PartType0", "particle_position_y"][sp_condensedpart_mask_for_analysis].to("kpc").value - halo_center_post[1].to("kpc").value
        recently_condensed_gas_z_for_analysis       = sp["PartType0", "particle_position_z"][sp_condensedpart_mask_for_analysis].to("kpc").value - halo_center_post[2].to("kpc").value
        recently_condensed_gas_mass_for_analysis    = sp["PartType0", "particle_mass"][sp_condensedpart_mask_for_analysis].to("Msun").value
        recently_condensed_gas_PIDs_for_analysis_found_in_sphere = sp["PartType0", "unique_PID"][sp_condensedpart_mask_for_analysis].value  ### This can be fewer than recently_condensed_gas_PIDs_for_analysis because of the aditional temperature
        recently_condensed_gas_radii_for_analysis   = sp["PartType0", "particle_radius"][sp_condensedpart_mask_for_analysis].to("kpc").value
        
        # if len(recently_condensed_gas_PIDs_for_analysis_found_in_sphere) != recently_condensed_gas_PIDs_for_analysis:
        #     print("Warning! len(recently_condensed_gas_PIDs_for_analysis_found_in_sphere) != recently_condensed_gas_PIDs_for_analysis did you use different sphere sizes between scripts?")
            
        #     if not ALLOW_FILES_MISSING_FOR_TESTING:
        #         raise RuntimeError()
        

        
        

        
        
        
        # n_bins = 100
        # ptype = "filtered_gas"
        # radius_field = (ptype, "particle_radius")
        # filtered_gas_entropy_profile = yt.create_profile(sp,radius_field, extrema = {radius_field:(1e-5*halo_r500.to('kpc'), 1.05*NUM_R500_FOR_PROJ*halo_r500.to('kpc'))},
        #                         fields=[(ptype, 'xray_entropy')],
        #                         units={radius_field: "kpc", (ptype, 'xray_entropy'):"cm**2*keV"},logs={radius_field: True},
        #                         weight_field = (ptype, "particle_mass"),  ### need this, or else is an accumulation
        #                         accumulation = False,
        #                         n_bins = n_bins)  
        

        
        
        # ptype = "gas"
        # radius_field = (ptype, "particle_radius")
        # all_gas_entropy_profile = yt.create_profile(sp,radius_field, extrema = {radius_field:(1e-5*halo_r500.to('kpc'), 1.05*NUM_R500_FOR_PROJ*halo_r500.to('kpc'))},
        #                         fields=[(ptype, 'xray_entropy')],
        #                         units={radius_field: "kpc", (ptype, 'xray_entropy'):"cm**2*keV"},logs={radius_field: True},
        #                         weight_field = (ptype, "particle_mass"),  ### need this, or else is an accumulation
        #                         accumulation = False,
        #                         n_bins = n_bins)  
        

        
        
# def _xray_entropy(field, data):
#     ''' https://ned.ipac.caltech.edu/level5/Sept09/Bohringer/Bohringer3.html'''
#     return unyt.physical_constants.kb * data["gas", "temperature"] / data["gas", "El_number_density"]**(2/3)
# yt.add_field(
#     name=("gas", "xray_entropy"),
#     function=_xray_entropy,
#     sampling_type="local",
#     units="keV* cm**2",
# )



    '''#############################'''
    '''                             '''
    '''                             '''
    '''Make entropy radial profiles '''
    '''                             '''
    '''                             '''
    '''                             '''
    '''#############################'''

    if DO_ENTROPY_PROFILES:   
        ''' For some reason entropy won't be added by yT to filtered_gas... lets just do it manually '''
        
        print(RANK, "doing binned entropy")
        entropy_radius_bins = np.linspace(-3,np.log10(NUM_R500_FOR_PROJ), 500)
        binned_entropy_all_gas_medians, entropy_bin_edges, _ = binned_statistic(np.log10(sp["gas","radius"].to("kpc").value/halo_r500),
                                                        (unyt.physical_constants.kb * sp["gas", "temperature"] / sp["gas", "El_number_density"]**(2/3)).to("keV*cm**2").value, bins = entropy_radius_bins)
        
        binned_entropy_filtered_gas_medians, _, _ = binned_statistic(np.log10(sp["filtered_gas","radius"].to("kpc").value/halo_r500),
                                                        (unyt.physical_constants.kb * sp["filtered_gas", "temperature"] / sp["filtered_gas", "El_number_density"]**(2/3)).to("keV*cm**2").value, bins = entropy_radius_bins)
        os.makedirs(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/ENTROPY_PROFILES/", exist_ok=True)
        entropy_fig = plt.figure(figsize=(6,4))
        plt.plot( np.asarray(entropy_bin_edges[:-1]) + 0.5*np.diff(entropy_bin_edges)[0], np.log10(binned_entropy_filtered_gas_medians), )
        plt.xlabel("Log R/R500")
        plt.ylabel("Log Sx")
        plt.xlim(-1.8, 0.5)
        plt.ylim(0.2, 6)
        plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/ENTROPY_PROFILES/zm{zoomnum}_Level{zoomlevel}_{str(snapnum_this_rank).zfill(3)}_xray_entropy_filtered_gas.png")
        plt.close(entropy_fig)
        
        
        
    '''#################################'''
    '''                                 '''
    '''                                 '''
    '''Save all required data to h5 file'''
    '''                                 '''
    '''                                 '''
    '''                                 '''
    '''#################################'''
    if DO_ANALYSIS:
        os.makedirs(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PRODUCTS/", exist_ok=True)
        f = h5py.File(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PRODUCTS/zm{zoomnum}_Level{zoomlevel}_sn{str(snapnum_this_rank).zfill(3)}_data.h5",'w')

        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/masses", data=instantaneously_condensed_gas_mass_for_analysis)
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/radii",  data=instantaneously_condensed_gas_radii_for_analysis)
        
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_bubbleparts_MASK", data = list(close_to_bubbleparts_MASK))
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/condensed_bubble_material_MASK", data = list(condensed_bubble_material_MASK))
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_infront_MASK ", data = list(close_to_satellite_infront_MASK) )
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_wake_MASK", data = list(close_to_satellite_wake_MASK))
        
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_bubbleparts_IDs", data =[int(x.split("_")[0]) for x in close_to_bubbleparts_PIDs])
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_bubbleparts_Unique_Identifier", data =[int(x.split("_")[1]) for x in close_to_bubbleparts_PIDs])
        if len(close_to_bubbleparts_PIDs) > 0:
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_bubbleparts_masses", data =close_to_bubbleparts_masses)
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_bubbleparts_radii", data =close_to_bubbleparts_radii)
        else:
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_bubbleparts_masses", data = 0)
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_bubbleparts_radii", data = 0 )
        
        
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/condensed_bubble_material_IDs", data =[int(x.split("_")[0]) for x in condensed_bubble_material_PIDs])
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/condensed_bubble_material_Unique_Identifier", data =[int(x.split("_")[1]) for x in condensed_bubble_material_PIDs])
        if len(condensed_bubble_material_PIDs) > 0:
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/condensed_bubble_material_masses", data =condensed_bubble_material_masses)
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/condensed_bubble_material_radii", data =condensed_bubble_material_radii)
        else:
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/condensed_bubble_material_masses", data = 0)
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/condensed_bubble_material_radii", data = 0 )
        
        
        
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_wake_IDs", data =[int(x.split("_")[0]) for x in close_to_satellite_wake_PIDs])
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_wake_Unique_Identifier", data =[int(x.split("_")[1]) for x in close_to_satellite_wake_PIDs])
        if len(close_to_satellite_wake_PIDs) > 0:
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_wake_masses", data =close_to_satellite_wake_masses)
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_wake_radii", data =close_to_satellite_wake_radii)
        else:
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_wake_masses", data = 0)
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_wake_radii", data = 0 )
            
        
        
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_infront_IDs", data =[int(x.split("_")[0]) for x in close_to_satellite_infront_PIDs])
        f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_infront_Unique_Identifier", data =[int(x.split("_")[1]) for x in close_to_satellite_infront_PIDs])
        if len(close_to_satellite_infront_PIDs) > 0:
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_infront_masses", data =close_to_satellite_infront_masses)
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_infront_radii", data =close_to_satellite_infront_radii)
        else:
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_infront_masses", data = 0)
            f.create_dataset(f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum_this_rank}/close_to_satellite_infront_radii", data = 0 )
        
        f.create_dataset(f"CONDENSED_PARTICLES/{snapnum_this_rank}/current_time",  data=ds.current_time.to("Myr").value)
        
        f.create_dataset(f"RECENTLY_CONDENSED_PARTICLES/{snapnum_this_rank}/masses", data=recently_condensed_gas_mass_for_analysis)
        f.create_dataset(f"RECENTLY_CONDENSED_PARTICLES/{snapnum_this_rank}/radii",  data=recently_condensed_gas_radii_for_analysis)
        f.create_dataset(f"RECENTLY_CONDENSED_PARTICLES/{snapnum_this_rank}/IDs",  data=[int(x.split("_")[0]) for x in recently_condensed_gas_PIDs_for_analysis_found_in_sphere])
        f.create_dataset(f"RECENTLY_CONDENSED_PARTICLES/{snapnum_this_rank}/Unique_Identifier",  data=[int(x.split("_")[1]) for x in recently_condensed_gas_PIDs_for_analysis])
        f.create_dataset(f"RECENTLY_CONDENSED_PARTICLES/{snapnum_this_rank}/ANALYSE_CONDENSED_PARTS_FROM_PREVIOUS_N_SNAPS",  data=ANALYSE_CONDENSED_PARTS_FROM_PREVIOUS_N_SNAPS)
        
        
        if DO_ENTROPY_PROFILES:
            f.create_dataset(f"RADIAL_PROFILES/{snapnum_this_rank}/FILTERED_GAS_ENTROPY",  data=binned_entropy_filtered_gas_medians)
            f.create_dataset(f"RADIAL_PROFILES/{snapnum_this_rank}/FILTERED_GAS_ENTROPY_x",  data=np.asarray(entropy_bin_edges[:-1]) + 0.5*np.diff(entropy_bin_edges)[0]) ### save bin midpoints
            
            f.create_dataset(f"RADIAL_PROFILES/{snapnum_this_rank}/ALL_GAS_ENTROPY",  data= binned_entropy_all_gas_medians)
            f.create_dataset(f"RADIAL_PROFILES/{snapnum_this_rank}/ALL_GAS_ENTROPY_x",  data=np.asarray(entropy_bin_edges[:-1]) + 0.5*np.diff(entropy_bin_edges)[0]) ### save bin midpoints
            
        
        
        ''' Bubbles'''
        
        f.create_dataset(f"BUBBLE_PARTS/{snapnum_this_rank}/recent_still_hot_bubble_particle_masses_arr", data = recent_still_hot_bubble_particle_masses_arr)      
        f.create_dataset(f"BUBBLE_PARTS/{snapnum_this_rank}/recent_still_hot_bubble_particle_radii_arr", data = recent_still_hot_bubble_particle_radii_arr)       
        f.create_dataset(f"BUBBLE_PARTS/{snapnum_this_rank}/recent_still_hot_bubble_particle_ages_arr", data = recent_still_hot_bubble_particle_ages_arr)        
        f.create_dataset(f"BUBBLE_PARTS/{snapnum_this_rank}/recent_still_hot_bubble_particle_temperature_arr", data = recent_still_hot_bubble_particle_temperature_arr) 
        
        f.create_dataset(f"BUBBLE_PARTS/{snapnum_this_rank}/recent_cooled_bubble_particle_masses_arr", data = recent_cooled_bubble_particle_masses_arr)         
        f.create_dataset(f"BUBBLE_PARTS/{snapnum_this_rank}/recent_cooled_bubble_particle_radii_arr", data = recent_cooled_bubble_particle_radii_arr)          
        f.create_dataset(f"BUBBLE_PARTS/{snapnum_this_rank}/recent_cooled_bubble_particle_ages_arr", data = recent_cooled_bubble_particle_ages_arr)           
        f.create_dataset(f"BUBBLE_PARTS/{snapnum_this_rank}/recent_cooled_bubble_particle_temperature_arr", data = recent_cooled_bubble_particle_temperature_arr)    
        
        f.create_dataset(f"BUBBLE_PARTS/{snapnum_this_rank}/launched_this_snap_bubble_gas_radius", data = launched_this_snap_bubble_gas_radius)             
        f.create_dataset(f"BUBBLE_PARTS/{snapnum_this_rank}/launched_this_snap_bubble_gas_mass", data = launched_this_snap_bubble_gas_mass)               
        f.create_dataset(f"BUBBLE_PARTS/{snapnum_this_rank}/launched_this_snap_bubble_gas_temperature", data = launched_this_snap_bubble_gas_temperature)        
        
        f.close()
        
        
        
        
        
        # os.makedirs(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_TESTING/", exist_ok=True)
        # fig = plt.figure(figsize = (6,4))
        
        # test_recently_condensed_masses_present_in_plot = []
        # test_snapnums_arr = []
        
        # for snapnum in ORIGINAL_SNAPNUMS_WHOLE_RANGE:
        #     condensed_snaps_to_plot = []
        #     condensed_gas_PIDs_to_plot = []
        #     for condensed_i in range(len(condensed_gas_unique_PIDs)):
        #         if condensed_gas_snapshots_when_cold[condensed_i] <= snapnum and condensed_gas_snapshots_when_cold[condensed_i] > snapnum - ANALYSE_CONDENSED_PARTS_FROM_PREVIOUS_N_SNAPS:
        #             condensed_gas_PIDs_to_plot.append(condensed_gas_unique_PIDs[condensed_i])
        #             condensed_snaps_to_plot.append(condensed_gas_snapshots_when_cold[condensed_i])
                
        #     sp_condensedpart_mask_for_plot = np.in1d(sp["PartType0", "unique_PID"].d, condensed_gas_PIDs_to_plot)
            
        #     if np.sum(sp_condensedpart_mask_for_plot) > 0:
        #         condensed_gas_mass_for_plot = float(np.sum(sp["PartType0", "particle_mass"][sp_condensedpart_mask_for_plot].to("Msun").value))
        #     else:
        #         condensed_gas_mass_for_plot = 1e-1
            
        #     test_recently_condensed_masses_present_in_plot.append(np.log10(condensed_gas_mass_for_plot))
            
        #     test_snapnums_arr.append(snapnum)
        # plt.plot(test_snapnums_arr, test_recently_condensed_masses_present_in_plot)
        # plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_TESTING/test_loading_condensed_parts_same_across_all_ranks_{snapnum_this_rank}.png")
        
        # plt.close()   
COMM.Barrier()




'''#####################################################'''
'''                                                     '''
'''                                                     '''
'''All ranks have finished and saved their data for each'''
''' required snapnum. On RANK 0 combine this data and   '''
'''                  perform analysis                   '''
'''                                                     '''
'''                                                     '''
'''                                                     '''
'''#####################################################'''

if (GO_STRAIGHT_TO_PLOTTING_PRESAVED_ANALYSIS or PLOT_ANALYSIS): # and RANK == 0:
    
    instantaneous_all_masses        = np.asarray([])
    instantaneous_all_radii         = np.asarray([])
    instantaneously_condensed_gas_all_times_Myr                  = np.asarray([])
    single_times_arr                = np.asarray([])
    total_condensed_mass_in_snap    = np.asarray([])
    
    ''' Now for RECENTLY condensed (not just in current snapshot) particles'''
    recent_all_masses               = np.asarray([])
    recent_all_radii                = np.asarray([])
    recent_condensed_gas_all_times_Myr = np.asarray([])
    
    ''' entropy'''
    all_gas_entropy_list_to_stack           = []
    filtered_gas_entropy_list_to_stack      = []
    
    
    
    ''' PID arrays'''
    close_to_bubbleparts_Unique_PIDs_arr        = np.asarray([])
    condensed_bubble_material_Unique_PIDs_arr   = np.asarray([])
    close_to_satellite_wake_Unique_PIDs_arr     = np.asarray([])
    close_to_satellite_infront_Unique_PIDs_arr  = np.asarray([])
    recently_condensed_gas_Unique_PIDs_arr      = np.asarray([])
    
    close_to_bubbleparts_times                  = np.asarray([])
    close_to_satellite_wake_times               = np.asarray([])
    close_to_satellite_infront_times            = np.asarray([])
    
    
    ''' Bubbles'''
    recent_still_hot_bubble_particle_masses_arr       = np.asarray([])
    recent_still_hot_bubble_particle_radii_arr        = np.asarray([])
    recent_still_hot_bubble_particle_ages_arr         = np.asarray([])
    recent_still_hot_bubble_particle_temperature_arr  = np.asarray([])
    recent_cooled_bubble_particle_masses_arr          = np.asarray([])
    recent_cooled_bubble_particle_radii_arr           = np.asarray([])
    recent_cooled_bubble_particle_ages_arr            = np.asarray([])
    recent_cooled_bubble_particle_temperature_arr     = np.asarray([])
    launched_this_snap_bubble_gas_radius              = np.asarray([])
    launched_this_snap_bubble_gas_mass                = np.asarray([])
    launched_this_snap_bubble_gas_temperature         = np.asarray([])
    recent_still_hot_bubble_particle_times_arr        = np.asarray([])
    recent_cooled_bubble_particle_times_arr           = np.asarray([])
                
    
    
    snapnums_read = []
    
    
    
    ''' General '''
    instantaneous_close_to_bubbleparts_total_mass_arr        = []
    instantaneous_condensed_bubble_material_total_mass_arr   = []
    instantaneous_close_to_satellite_wake_total_mass_arr     = []
    instantaneous_close_to_satellite_infront_total_mass_arr  = []
    
    
    
    print("ORIGINAL_SNAPNUMS_WHOLE_RANGE", ORIGINAL_SNAPNUMS_WHOLE_RANGE, len(ORIGINAL_SNAPNUMS_WHOLE_RANGE))
    os.makedirs(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/INDIVIDUAL_SNAPNUMS/",exist_ok=True)
    for snapnum in ORIGINAL_SNAPNUMS_WHOLE_RANGE:
        file = f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PRODUCTS/zm{zoomnum}_Level{zoomlevel}_sn{str(snapnum).zfill(3)}_data.h5"
        # 
        print("Loading for snapnum ", snapnum, flush = True)
        if 1==1:
            
            
            if not os.path.exists(file) and ALLOW_FILES_MISSING_FOR_TESTING:
                print("Warning!", file," does not exist!")
                continue
            
            try:
                with h5py.File(file,'r') as f:
                    pass
            except Exception as e:
                print("Warning!", file," has a problem and cannot be read!")
                print(e)
                continue
                
            
            with h5py.File(file,'r') as f:
                
                # print(np.asarray(f[f"CONDENSED_PARTICLES/{snapnum}/current_time"]))
                # continue
                # print(list(f.keys()))
                # continue
                try:
                    print('f[f"CONDENSED_PARTICLES/', list(f[f"CONDENSED_PARTICLES/{snapnum}/"].keys()))
                    current_time_minus_start_time = float(float( np.asarray(f[f"CONDENSED_PARTICLES/{snapnum}/current_time"])) - start_time)
                
                    single_times_arr                                            = np.append(single_times_arr, np.asarray(f[f"CONDENSED_PARTICLES/{snapnum}/current_time"]))
                except Exception as e:
                    print("Warning 1!", file," has a problem and cannot be read!")
                    print(e)
                    # exit()
                    continue
                
                
                
                if DO_INSTANTANEOUSLY_CONDENSED_PLOTS:
                    instantaneous_condensed_gas_mass_for_analysis_this_snap     = np.asarray(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/masses"])
                    instantaneous_condensed_gas_radii_for_analysis_this_snap    = np.asarray(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/radii"])
                    instantaneous_all_masses                                    = np.concatenate((instantaneous_all_masses, instantaneous_condensed_gas_mass_for_analysis_this_snap))
                    instantaneous_all_radii                                     = np.concatenate((instantaneous_all_radii, instantaneous_condensed_gas_radii_for_analysis_this_snap))
                    
                    
                    
                    # exit()
                    if type(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_bubbleparts_masses"][()] ) == int :
                        instantaneous_close_to_bubbleparts_total_mass_arr.append(0)
                    else:
                        instantaneous_close_to_bubbleparts_total_mass_arr.append(          np.sum(np.asarray(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_bubbleparts_masses"])))      

                    if type(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/condensed_bubble_material_masses"][()]  ) == int :
                        instantaneous_condensed_bubble_material_total_mass_arr.append(0)
                    else:  
                        instantaneous_condensed_bubble_material_total_mass_arr.append(     np.sum(np.asarray(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/condensed_bubble_material_masses"]))) 

                    if type(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_satellite_wake_masses"][()]    ) == int :
                        instantaneous_close_to_satellite_wake_total_mass_arr.append(0)
                    else:  
                        instantaneous_close_to_satellite_wake_total_mass_arr.append(       np.sum(np.asarray(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_satellite_wake_masses"])))   

                    if type(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_satellite_infront_masses"][()] ) == int :
                        instantaneous_close_to_satellite_infront_total_mass_arr.append(0)
                    else:  
                        instantaneous_close_to_satellite_infront_total_mass_arr.append(    np.sum(np.asarray(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_satellite_infront_masses"])))

                    
                    
                    
                    total_condensed_mass_in_snap                                = np.concatenate((total_condensed_mass_in_snap, np.asarray([np.sum(instantaneous_condensed_gas_mass_for_analysis_this_snap),])))
                
                recent_condensed_gas_mass_for_analysis_this_snap            = np.asarray(f[f"RECENTLY_CONDENSED_PARTICLES/{snapnum}/masses"])
                recent_condensed_gas_radii_for_analysis_this_snap           = np.asarray(f[f"RECENTLY_CONDENSED_PARTICLES/{snapnum}/radii"])
                recent_all_masses                                           = np.concatenate((recent_all_masses, recent_condensed_gas_mass_for_analysis_this_snap))
                recent_all_radii                                            = np.concatenate((recent_all_radii, recent_condensed_gas_radii_for_analysis_this_snap))


                recently_condensed_gas_IDs                                  =     list(f[f"RECENTLY_CONDENSED_PARTICLES/{snapnum}/IDs"])
                recently_condensed_gas_Unique_Identifier                    =     list(f[f"RECENTLY_CONDENSED_PARTICLES/{snapnum}/Unique_Identifier"])
                recently_condensed_gas_Unique_PIDs          =   np.asarray(["{}_{}".format(a_, b_) for a_, b_ in zip(recently_condensed_gas_IDs,    recently_condensed_gas_Unique_Identifier)] )









                if DO_CLOSE_TO_WAKE_AND_BUBBLES_PLOTS:
                    # close_to_bubbleparts_PIDs_condensed_this_snap               = list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_bubbleparts_PIDs"])
                    # condensed_bubble_material_PIDs_condensed_this_snap          = list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/condensed_bubble_material_PIDs"])
                    # close_to_satellite_wake_PIDs_condensed_this_snap            = list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_satellite_wake_PIDs"])
                    # close_to_satellite_infront_PIDs_condensed_this_snap         = list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_satellite_infront_PIDs"])
                    close_to_bubbleparts_IDs                                    =     list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_bubbleparts_IDs"])
                    condensed_bubble_material_IDs                               =     list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/condensed_bubble_material_IDs"])
                    close_to_satellite_wake_IDs                                 =     list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_satellite_wake_IDs"])
                    close_to_satellite_infront_IDs                              =     list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_satellite_infront_IDs"])
                    close_to_bubbleparts_Unique_Identifier                      =     list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_bubbleparts_Unique_Identifier"])
                    condensed_bubble_material_Unique_Identifier                 =     list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/condensed_bubble_material_Unique_Identifier"])
                    close_to_satellite_wake_Unique_Identifier                   =     list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_satellite_wake_Unique_Identifier"])
                    close_to_satellite_infront_Unique_Identifier                =     list(f[f"INSTANTANEOUSLY_CONDENSED_PARTICLES/{snapnum}/close_to_satellite_infront_Unique_Identifier"])

    
                    close_to_bubbleparts_Unique_PIDs            =   np.asarray(["{}_{}".format(a_, b_) for a_, b_ in zip(close_to_bubbleparts_IDs,    close_to_bubbleparts_Unique_Identifier)] )
                    condensed_bubble_material_Unique_PIDs       =   np.asarray(["{}_{}".format(a_, b_) for a_, b_ in zip(condensed_bubble_material_IDs,    condensed_bubble_material_Unique_Identifier)] )
                    close_to_satellite_wake_Unique_PIDs         =   np.asarray(["{}_{}".format(a_, b_) for a_, b_ in zip(close_to_satellite_wake_IDs,    close_to_satellite_wake_Unique_Identifier)] )
                    close_to_satellite_infront_Unique_PIDs      =   np.asarray(["{}_{}".format(a_, b_) for a_, b_ in zip(close_to_satellite_infront_IDs,    close_to_satellite_infront_Unique_Identifier)] )
                
                    close_to_bubbleparts_Unique_PIDs_arr                = np.concatenate((close_to_bubbleparts_Unique_PIDs_arr, close_to_bubbleparts_Unique_PIDs       )       )
                    condensed_bubble_material_Unique_PIDs_arr           = np.concatenate((condensed_bubble_material_Unique_PIDs_arr, condensed_bubble_material_Unique_PIDs  ) ) 
                    close_to_satellite_wake_Unique_PIDs_arr             = np.concatenate((close_to_satellite_wake_Unique_PIDs_arr, close_to_satellite_wake_Unique_PIDs    )   ) 
                    close_to_satellite_infront_Unique_PIDs_arr          = np.concatenate((close_to_satellite_infront_Unique_PIDs_arr, close_to_satellite_infront_Unique_PIDs )) 
                
        
                # print(len(recent_condensed_gas_mass_for_analysis_this_snap), len(recently_condensed_gas_Unique_PIDs))
                assert len(recent_condensed_gas_mass_for_analysis_this_snap) == len(recently_condensed_gas_Unique_PIDs)
   


                    
                # print("recently_condensed_gas_Unique_PIDs", recently_condensed_gas_Unique_PIDs)
                recently_condensed_gas_Unique_PIDs_arr              = np.concatenate((recently_condensed_gas_Unique_PIDs_arr, recently_condensed_gas_Unique_PIDs     )    ) 
                # print("recently_condensed_gas_Unique_PIDs_arr",recently_condensed_gas_Unique_PIDs_arr)
                # print("shape recently_condensed_gas_Unique_PIDs_arr", np.shape(recently_condensed_gas_Unique_PIDs_arr))
                
                
                ''' BUBBLE_DATA'''
                recent_still_hot_bubble_particle_masses_arr         = np.concatenate(( recent_still_hot_bubble_particle_masses_arr     ,np.asarray(f[f"BUBBLE_PARTS/{snapnum}/recent_still_hot_bubble_particle_masses_arr"])            ))
                recent_still_hot_bubble_particle_radii_arr          = np.concatenate(( recent_still_hot_bubble_particle_radii_arr      ,np.asarray(f[f"BUBBLE_PARTS/{snapnum}/recent_still_hot_bubble_particle_radii_arr"])            ))
                recent_still_hot_bubble_particle_ages_arr           = np.concatenate(( recent_still_hot_bubble_particle_ages_arr       ,np.asarray(f[f"BUBBLE_PARTS/{snapnum}/recent_still_hot_bubble_particle_ages_arr"])            ))
                recent_still_hot_bubble_particle_temperature_arr    = np.concatenate(( recent_still_hot_bubble_particle_temperature_arr,np.asarray(f[f"BUBBLE_PARTS/{snapnum}/recent_still_hot_bubble_particle_temperature_arr"])            ))
                recent_cooled_bubble_particle_masses_arr            = np.concatenate(( recent_cooled_bubble_particle_masses_arr        ,np.asarray(f[f"BUBBLE_PARTS/{snapnum}/recent_cooled_bubble_particle_masses_arr"])        ))
                recent_cooled_bubble_particle_radii_arr             = np.concatenate(( recent_cooled_bubble_particle_radii_arr         ,np.asarray(f[f"BUBBLE_PARTS/{snapnum}/recent_cooled_bubble_particle_radii_arr"])        ))
                recent_cooled_bubble_particle_ages_arr              = np.concatenate(( recent_cooled_bubble_particle_ages_arr          ,np.asarray(f[f"BUBBLE_PARTS/{snapnum}/recent_cooled_bubble_particle_ages_arr"])        ))
                # recent_cooled_bubble_particle_temperature_arr       = np.concatenate(( recent_cooled_bubble_particle_temperature_arr   ,np.asarray(f[f"BUBBLE_PARTS/{snapnum}/recent_cooled_bubble_particle_temperature_arr"])))
                launched_this_snap_bubble_gas_radius                = np.concatenate(( launched_this_snap_bubble_gas_radius            ,np.asarray(f[f"BUBBLE_PARTS/{snapnum}/launched_this_snap_bubble_gas_radius"])        ))
                launched_this_snap_bubble_gas_mass                  = np.concatenate(( launched_this_snap_bubble_gas_mass              ,np.asarray(f[f"BUBBLE_PARTS/{snapnum}/launched_this_snap_bubble_gas_mass"])        ))
                launched_this_snap_bubble_gas_temperature           = np.concatenate(( launched_this_snap_bubble_gas_temperature       ,np.asarray(f[f"BUBBLE_PARTS/{snapnum}/launched_this_snap_bubble_gas_temperature"])        ))
                
                recent_still_hot_bubble_particle_times_arr          = np.concatenate((recent_still_hot_bubble_particle_times_arr, current_time_minus_start_time * np.ones_like(np.asarray(f[f"BUBBLE_PARTS/{snapnum}/recent_still_hot_bubble_particle_masses_arr"])            ) ))
                recent_cooled_bubble_particle_times_arr             = np.concatenate((recent_cooled_bubble_particle_times_arr, current_time_minus_start_time * np.ones_like(np.asarray(f[f"BUBBLE_PARTS/{snapnum}/recent_cooled_bubble_particle_masses_arr"])) ))
                
                
                if DO_CLOSE_TO_WAKE_AND_BUBBLES_PLOTS:
                    ''' record the snapshot where particles have condensed near to bubbles'''
                    close_to_bubbleparts_time                       = current_time_minus_start_time * np.ones(len(close_to_bubbleparts_Unique_PIDs))
                    close_to_bubbleparts_times  = np.concatenate((close_to_bubbleparts_times, close_to_bubbleparts_time ))
                    close_to_wake_time                       = current_time_minus_start_time * np.ones(len(close_to_satellite_wake_Unique_PIDs))
                    close_to_infront_time                       = current_time_minus_start_time * np.ones(len(close_to_satellite_infront_Unique_PIDs))
                    close_to_satellite_wake_times  = np.concatenate((close_to_satellite_wake_times, close_to_wake_time ))
                    close_to_satellite_infront_times  = np.concatenate((close_to_satellite_infront_times, close_to_infront_time ))
                
                ''' record the snapshot where particles have condensed near to wake'''
                
                
                
                
                
                
                
                if DO_INSTANTANEOUSLY_CONDENSED_PLOTS:
                    instantaneously_condensed_gas_times_in_Myr              = current_time_minus_start_time * np.ones_like(instantaneous_condensed_gas_mass_for_analysis_this_snap)
                    instantaneously_condensed_gas_all_times_Myr             = np.concatenate((instantaneously_condensed_gas_all_times_Myr, instantaneously_condensed_gas_times_in_Myr))
                recent_condensed_gas_times_in_Myr                       = current_time_minus_start_time * np.ones_like(recent_condensed_gas_mass_for_analysis_this_snap)
                recent_condensed_gas_all_times_Myr                      = np.concatenate((recent_condensed_gas_all_times_Myr, recent_condensed_gas_times_in_Myr))




                
   
   
   
                ''' ENTROPY '''
                if DO_ENTROPY_PROFILES:
                    all_gas_entropy_list_to_stack.append(list(f[f"RADIAL_PROFILES/{snapnum}/ALL_GAS_ENTROPY"]))    
                    filtered_gas_entropy_list_to_stack.append(list(f[f"RADIAL_PROFILES/{snapnum}/FILTERED_GAS_ENTROPY"]))
                    if snapnum == ORIGINAL_SNAPNUMS_WHOLE_RANGE[0]:
                        entropy_radial_profiles_xvals = np.asarray(f[f"RADIAL_PROFILES/{snapnum}/ALL_GAS_ENTROPY_x"])
                        
                    else:
                        ''' check that all entropy radial x values are homogeneous, across snaps'''
                        if not np.array_equal(entropy_radial_profiles_xvals,np.asarray(f[f"RADIAL_PROFILES/{snapnum}/ALL_GAS_ENTROPY_x"])):
                            raise NotImplementedError()
                        if not np.array_equal(entropy_radial_profiles_xvals,np.asarray(f[f"RADIAL_PROFILES/{snapnum}/FILTERED_GAS_ENTROPY_x"])):
                            raise NotImplementedError()
                    
        # except Exception as e:
        #     if not ALLOW_FILES_MISSING_FOR_TESTING:
        #         raise RuntimeError(file,e)   
        #     else:
        #         print(f"{file}, \n",e, "\n continuing since ALLOW_FILES_MISSING_FOR_TESTING is switched on")
        #         # os.remove(file)

                snapnums_read.append(snapnum)

    COMM.Barrier()
    
    
            
    define_units(cosmo)
    ''' Load the BH wind and properties histories'''
    bhdetails_dir = f"/home/b/babul/fjenning/scratch/xray_bubbles/PAPER2_Level{zoomlevel}_zoom_{zoomnum}/" + "/snaps/blackhole_details/"


    bh_wind_files = [bhdetails_dir + f"/{f}" for f in os.listdir(bhdetails_dir) if "bhwind" in f]
    bh_details_files = [bhdetails_dir + f"/{f}" for f in os.listdir(bhdetails_dir) if "blackhole_details" in f]  
    bhwind_history = get_bhwind_history(SMBH_pid,  files= bh_wind_files, min_a = start_scale_factor-0.1, max_a = 1)
    if bhwind_history == []:
        print("No BH winds this snap")
        bhwind_PIDs = []
    else:
        test_units=True
        bhinfo_history  = get_bh_info(SMBH_pid,  files=bh_details_files, min_a = start_scale_factor-0.1, max_a = 1, a_list_winds = bhwind_history["a_list_winds"], )
    
        if RUN_TESTS: check_dict_units(bhwind_history, "bhwind_history")
        if RUN_TESTS: check_dict_units(bhinfo_history, "bhinfo_history")
        
        
        if 1 < max(bhwind_history["a_list_winds"]) or start_scale_factor-0.1 > min(bhwind_history["a_list_winds"]):
            raise ValueError("Scale factor range for winds is bad")
        
        bhwind_PIDs = [int(id) for id in bhwind_history["bhwind_id_history"]]

    bhwind_history = add_relative_velocites(bhwind_history, bhinfo_history)
    bhwind_history["bhwind_vels"] = np.copy(bhwind_history["bhwind_bhrelative_vels"])
    min_vel = bhinfo_history["quasar_vel_from_mdot_history_at_bhwind_times"]
    
    print("min max quasar vel", min(min_vel), max(min_vel))
    jet_mask = [bhwind_history["bhwind_vels"][i].to("km/s").value > min_vel[i].to("km/s").value + JET_DELTA_V_FROM_QUASARS_inKm_per_s  for i in range(max(len(bhwind_history["bhwind_vels"]),len(min_vel)))]
    quasar_mask = [bhwind_history["bhwind_vels"][i].to("km/s").value < min_vel[i].to("km/s").value + JET_DELTA_V_FROM_QUASARS_inKm_per_s  for i in range(max(len(bhwind_history["bhwind_vels"]),len(min_vel)))]

    KE = (0.5 *  (bhwind_history["bhwind_vels"]**2)  * bhwind_history["bhwind_mass"]).to("ergs")
    bhwinds_Myr_times = cosmo.age([(1/a - 1) for a in bhwind_history['a_list_winds'] ]).to("Myr").value
    bhwinds_Myr_times = bhwinds_Myr_times #- start_time
    
    jet_times = bhwinds_Myr_times[jet_mask]
    jet_KE = KE[jet_mask]
    

    print("total_condensed_mass_in_snap",total_condensed_mass_in_snap)
    print("instantaneous_close_to_bubbleparts_total_mass_arr)",instantaneous_close_to_bubbleparts_total_mass_arr)
    print("instantaneous_condensed_bubble_material_total_mass_arr)",instantaneous_condensed_bubble_material_total_mass_arr)
    print("instantaneous_close_to_satellite_wake_total_mass_arr)",instantaneous_close_to_satellite_wake_total_mass_arr)
    print("instantaneous_close_to_satellite_infront_total_mass_arr",instantaneous_close_to_satellite_infront_total_mass_arr)
    fig, ax1 = plt.subplots() 
    fig.set_size_inches((12,6))
    # plt.plot(single_times_arr, )
    plt.plot(single_times_arr, np.log10(np.cumsum(np.asarray(total_condensed_mass_in_snap))), label = "Total")
    plt.plot(single_times_arr, np.log10(np.cumsum(np.asarray(instantaneous_close_to_bubbleparts_total_mass_arr))), label = "Close to Bubbles")
    plt.plot(single_times_arr, np.log10(np.cumsum(np.asarray(instantaneous_condensed_bubble_material_total_mass_arr))), label = "Condensed Bubble Material")
    plt.plot(single_times_arr, np.log10(np.cumsum(np.asarray(instantaneous_close_to_satellite_wake_total_mass_arr))), label = "Satellite Wake")
    plt.plot(single_times_arr, np.log10(np.cumsum(np.asarray(instantaneous_close_to_satellite_infront_total_mass_arr))), label = "Satellite In-Front")
    plt.xlabel("dt (Myr)")
    plt.ylabel("log M_condense")
    
    jets_mask_for_condensing_timescale = (jet_times >= np.min(single_times_arr))&(jet_times <= np.max(single_times_arr))
    ax2 = ax1.twinx()
    ax2.plot(jet_times[jets_mask_for_condensing_timescale], np.log10(np.cumsum(jet_KE[jets_mask_for_condensing_timescale].to("erg").value)), label = "JET")
    ax2.set_ylabel('Jet KE', color='b')
    
    ax1.legend()
    ax2.legend()
    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_time_vs_component_mass_condensing.png")
    plt.close()
    exit()
    
    
    
    
    
    


    if DO_CLOSE_TO_WAKE_AND_BUBBLES_PLOTS:
        recent_condensed_gas_all_times_Myr_just_from_near_bubbles   = []
        recent_condensed_gas_masses_just_from_near_bubbles          = []
        recent_condensed_gas_radii_just_from_near_bubbles           = []
        
        ranksplitted_recent_all_masses = np.array_split(recent_all_masses,SIZE)[RANK]
        ranksplitted_recent_all_radii = np.array_split(recent_all_radii,SIZE)[RANK]
        ranksplitted_recently_condensed_gas_Unique_PIDs_arr = np.array_split(recently_condensed_gas_Unique_PIDs_arr,SIZE)[RANK]
        rank_splitted_recent_condensed_gas_all_times_Myr = np.array_split(recent_condensed_gas_all_times_Myr, SIZE)[RANK]
        
        
        recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake   = []
        recent_condensed_gas_masses_just_from_near_satellite_wake          = []
        recent_condensed_gas_radii_just_from_near_satellite_wake           = []
        
        recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront   = []
        recent_condensed_gas_masses_just_from_near_satellite_infront          = []
        recent_condensed_gas_radii_just_from_near_satellite_infront           = []
        
        ranksplitted_recent_all_masses = np.array_split(recent_all_masses,SIZE)[RANK]
        ranksplitted_recent_all_radii = np.array_split(recent_all_radii,SIZE)[RANK]
        ranksplitted_recently_condensed_gas_Unique_PIDs_arr = np.array_split(recently_condensed_gas_Unique_PIDs_arr,SIZE)[RANK]
        rank_splitted_recent_condensed_gas_all_times_Myr = np.array_split(recent_condensed_gas_all_times_Myr, SIZE)[RANK]
        
        

        
        
        
        for i in range(len(ranksplitted_recent_all_masses)):
            if i % 1000 == 0 and RANK==0: print(i,"/",len(ranksplitted_recent_all_masses), flush = True)
            # print("recently_condensed_gas_Unique_PIDs_arr", recently_condensed_gas_Unique_PIDs_arr)
            gas_PID = ranksplitted_recently_condensed_gas_Unique_PIDs_arr[i]
            
            
            index_into_close_to_bubbleparts = np.where(close_to_bubbleparts_Unique_PIDs_arr == gas_PID)[0]
            if len(index_into_close_to_bubbleparts) == 0:  ### particle never condensed near to bubble
                continue
            elif close_to_bubbleparts_times[min(index_into_close_to_bubbleparts)] < rank_splitted_recent_condensed_gas_all_times_Myr[i]:   ### Gas particle condensed near to bubble before the current time
                recent_condensed_gas_masses_just_from_near_bubbles.append(ranksplitted_recent_all_masses[i])
                recent_condensed_gas_all_times_Myr_just_from_near_bubbles.append(rank_splitted_recent_condensed_gas_all_times_Myr[i])
                recent_condensed_gas_radii_just_from_near_bubbles.append(ranksplitted_recent_all_radii[i])
                
            index_into_close_to_satellite_wake = np.where(close_to_satellite_wake_Unique_PIDs_arr == gas_PID)[0]
            if len(index_into_close_to_satellite_wake) == 0:  ### particle never condensed near to wake
                continue
            elif close_to_satellite_wake_times[min(index_into_close_to_satellite_wake)] < rank_splitted_recent_condensed_gas_all_times_Myr[i]:   ### Gas particle condensed near to wake before the current time
                recent_condensed_gas_masses_just_from_near_satellite_wake.append(ranksplitted_recent_all_masses[i])
                recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake.append(rank_splitted_recent_condensed_gas_all_times_Myr[i])
                recent_condensed_gas_radii_just_from_near_satellite_wake.append(ranksplitted_recent_all_radii[i])
            
            
            index_into_close_to_satellite_infront = np.where(close_to_satellite_infront_Unique_PIDs_arr == gas_PID)[0]
            if len(index_into_close_to_satellite_infront) == 0:  ### particle never condensed near to infront
                continue
            elif close_to_satellite_infront_times[min(index_into_close_to_satellite_infront)] < rank_splitted_recent_condensed_gas_all_times_Myr[i]:   ### Gas particle condensed near to infront before the current time
                recent_condensed_gas_masses_just_from_near_satellite_infront.append(ranksplitted_recent_all_masses[i])
                recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront.append(rank_splitted_recent_condensed_gas_all_times_Myr[i])
                recent_condensed_gas_radii_just_from_near_satellite_infront.append(ranksplitted_recent_all_radii[i])
            
            
            
            
            
            
        recent_condensed_gas_all_times_Myr_just_from_near_bubbles = np.asarray(recent_condensed_gas_all_times_Myr_just_from_near_bubbles)
        recent_condensed_gas_masses_just_from_near_bubbles        = np.asarray(recent_condensed_gas_masses_just_from_near_bubbles       )
        recent_condensed_gas_radii_just_from_near_bubbles         = np.asarray(recent_condensed_gas_radii_just_from_near_bubbles)

        recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake = np.asarray(recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake)
        recent_condensed_gas_masses_just_from_near_satellite_wake        = np.asarray(recent_condensed_gas_masses_just_from_near_satellite_wake       )
        recent_condensed_gas_radii_just_from_near_satellite_wake         = np.asarray(recent_condensed_gas_radii_just_from_near_satellite_wake)
        
        recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront = np.asarray(recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront)
        recent_condensed_gas_masses_just_from_near_satellite_infront        = np.asarray(recent_condensed_gas_masses_just_from_near_satellite_infront       )
        recent_condensed_gas_radii_just_from_near_satellite_infront         = np.asarray(recent_condensed_gas_radii_just_from_near_satellite_infront)


        COMM.Barrier()
        if SIZE > 1:
            os.makedirs(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/TEMP_DATA_DIR/", exist_ok=True)
            with h5py.File(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/TEMP_DATA_DIR/{RANK}_analysis_temporary_file_for_triggering_of_condensation.h5",'w') as f:
                f.create_dataset(f"recent_condensed_gas_all_times_Myr_just_from_near_bubbles", data=recent_condensed_gas_all_times_Myr_just_from_near_bubbles)
                f.create_dataset(f"recent_condensed_gas_masses_just_from_near_bubbles", data=recent_condensed_gas_masses_just_from_near_bubbles)
                f.create_dataset(f"recent_condensed_gas_radii_just_from_near_bubbles", data=recent_condensed_gas_radii_just_from_near_bubbles)
                f.create_dataset(f"recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake", data=recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake)
                f.create_dataset(f"recent_condensed_gas_masses_just_from_near_satellite_wake", data=recent_condensed_gas_masses_just_from_near_satellite_wake)
                f.create_dataset(f"recent_condensed_gas_radii_just_from_near_satellite_wake", data=recent_condensed_gas_radii_just_from_near_satellite_wake)
                f.create_dataset(f"recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront", data=recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront)
                f.create_dataset(f"recent_condensed_gas_masses_just_from_near_satellite_infront", data=recent_condensed_gas_masses_just_from_near_satellite_infront)
                f.create_dataset(f"recent_condensed_gas_radii_just_from_near_satellite_infront", data=recent_condensed_gas_radii_just_from_near_satellite_infront)
            
            COMM.Barrier()
            
            if RANK == 0:
                recent_condensed_gas_all_times_Myr_just_from_near_bubbles = np.asarray([])
                recent_condensed_gas_masses_just_from_near_bubbles        = np.asarray([])
                recent_condensed_gas_radii_just_from_near_bubbles         = np.asarray([])
                recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake = np.asarray([])
                recent_condensed_gas_masses_just_from_near_satellite_wake        = np.asarray([])
                recent_condensed_gas_radii_just_from_near_satellite_wake         = np.asarray([])
                recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront = np.asarray([])
                recent_condensed_gas_masses_just_from_near_satellite_infront        = np.asarray([])
                recent_condensed_gas_radii_just_from_near_satellite_infront         = np.asarray([])
                for rank_i in range(SIZE):
                    file = f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/TEMP_DATA_DIR/{rank_i}_analysis_temporary_file_for_triggering_of_condensation.h5"
                    with h5py.File(file,'r') as f:
                        recent_condensed_gas_all_times_Myr_just_from_near_bubbles = np.concatenate((recent_condensed_gas_all_times_Myr_just_from_near_bubbles, np.asarray(f["recent_condensed_gas_all_times_Myr_just_from_near_bubbles"])))
                        recent_condensed_gas_masses_just_from_near_bubbles        = np.concatenate((recent_condensed_gas_masses_just_from_near_bubbles       , np.asarray(f["recent_condensed_gas_masses_just_from_near_bubbles"])))
                        recent_condensed_gas_radii_just_from_near_bubbles         = np.concatenate((recent_condensed_gas_radii_just_from_near_bubbles        , np.asarray(f["recent_condensed_gas_radii_just_from_near_bubbles"])))
                        recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake = np.concatenate((recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake, np.asarray(f["recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake"])))
                        recent_condensed_gas_masses_just_from_near_satellite_wake        = np.concatenate((recent_condensed_gas_masses_just_from_near_satellite_wake       , np.asarray(f["recent_condensed_gas_masses_just_from_near_satellite_wake"])))
                        recent_condensed_gas_radii_just_from_near_satellite_wake         = np.concatenate((recent_condensed_gas_radii_just_from_near_satellite_wake        , np.asarray(f["recent_condensed_gas_radii_just_from_near_satellite_wake"])))
                        recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront = np.concatenate((recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront, np.asarray(f["recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront"])))
                        recent_condensed_gas_masses_just_from_near_satellite_infront        = np.concatenate((recent_condensed_gas_masses_just_from_near_satellite_infront       , np.asarray(f["recent_condensed_gas_masses_just_from_near_satellite_infront"])))
                        recent_condensed_gas_radii_just_from_near_satellite_infront         = np.concatenate((recent_condensed_gas_radii_just_from_near_satellite_infront        , np.asarray(f["recent_condensed_gas_radii_just_from_near_satellite_infront"])))
                    
                    
                    os.remove(file)
            else:
                recent_condensed_gas_all_times_Myr_just_from_near_bubbles = None
                recent_condensed_gas_masses_just_from_near_bubbles        = None
                recent_condensed_gas_radii_just_from_near_bubbles         = None
                recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake = None
                recent_condensed_gas_masses_just_from_near_satellite_wake        = None
                recent_condensed_gas_radii_just_from_near_satellite_wake         = None
                recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront = None
                recent_condensed_gas_masses_just_from_near_satellite_infront        = None
                recent_condensed_gas_radii_just_from_near_satellite_infront         = None
        
            
    if RANK == 0:
        
        # print("recent_condensed_gas_all_times_Myr_just_from_near_bubbles", np.min(recent_condensed_gas_all_times_Myr_just_from_near_bubbles), np.max(recent_condensed_gas_all_times_Myr_just_from_near_bubbles))
        # print("recent_condensed_gas_all_times_Myr", np.min(recent_condensed_gas_all_times_Myr), np.max(recent_condensed_gas_all_times_Myr))
        # # exit()
        
        # print("recent_condensed_gas_all_times_Myr_just_from_near_bubbles", recent_condensed_gas_all_times_Myr_just_from_near_bubbles)
        # print("recent_condensed_gas_masses_just_from_near_bubbles", recent_condensed_gas_masses_just_from_near_bubbles)
        # # exit()
        if DO_INSTANTANEOUSLY_CONDENSED_PLOTS:
            plt.figure(figsize=(6,4))
            plt.plot(single_times_arr, np.log10(total_condensed_mass_in_snap))
            plt.xlabel("dt (Myr)")
            plt.ylabel("log M_condense")
            plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_time_vs_total_mass_condensing.png")
            plt.close()
            
            
            
            

        
        
        
        if PLOT_BHWIND_ENERGIES:

            
            
            
            KE_x_bins = np.linspace(min(recent_condensed_gas_all_times_Myr)-(DT/2),max(recent_condensed_gas_all_times_Myr)+(DT/2), int(len(set(recent_condensed_gas_all_times_Myr))//KE_x_bins_coarsening) + 1)
            binned_KE, KE_bin_edges, _ = scipy.stats.binned_statistic(x=jet_times, values = jet_KE, statistic="sum", bins=KE_x_bins,)
            KE_x_bin_midpoints = KE_bin_edges[:-1] + (np.diff(KE_bin_edges)[0]/2)
            
            
            
            
            for bin_i in [0.1,1,5,10,20]:
                fft_bin_size_Myr = bin_i
                KE_x_bins_fft = np.arange(min(recent_condensed_gas_all_times_Myr)-(DT/2),max(recent_condensed_gas_all_times_Myr)+(DT/2) + fft_bin_size_Myr, fft_bin_size_Myr)
                binned_KE_fft, KE_bin_edges_fft, _ = scipy.stats.binned_statistic(x=jet_times, values = jet_KE - np.median(jet_KE), statistic="sum", bins=KE_x_bins_fft,)
                KE_x_bin_midpoints_fft = KE_bin_edges_fft[:-1] + (np.diff(KE_bin_edges_fft)[0]/2)
                
                
                
                
                ''' fftfreq gives
                    f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
                    f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
                    in cycles/original-unit e.g. cycles per second
                    in rads multiply by 2pi to get 2npi/L   where L is the total width of the original data
                '''
                
                
                KE_fft_freqspace    = scipy.fft.fft(binned_KE_fft, n = len(KE_x_bins_fft))
                KE_fft_freqs        = scipy.fft.fftfreq(len(KE_x_bins_fft), d = fft_bin_size_Myr)
                
                
                plt.figure(figsize=(6,4))
                
                # plt.plot(2* np.pi * KE_fft_freqs, KE_fft_freqspace)
                plt.plot(2* np.pi * KE_fft_freqs[0:len(KE_x_bins_fft)//2], KE_fft_freqspace[0:len(KE_x_bins_fft)//2], color = "blue")
                plt.xlabel("2 pi/L  (rads/Myr)")
                plt.ylabel("KE Fourier Frequency")
                plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_Fourier_Frequency_JetKE_{bin_i}.png")
                plt.close()
            # exit()
            
            
            
            
            
            
            
            
            
            
            
            
            # print(jet_times)
            # print(jet_KE)
            # exit()
        
        '''' Entropy plots'''
        # my_cmap = plt.cm.jet
        # my_cmap.set_under('black',1)
        filtered_gas_entropy_stacked =  np.column_stack(filtered_gas_entropy_list_to_stack)
        # print("all_gas_entropy_stacked", np.log10(all_gas_entropy_stacked[:10,:]).round(decimals=1))
        
        # for i in range(3):
        #     # print(filtered_gas_entropy_stacked[i])


        # print("entropy_radial_profiles_xvals", entropy_radial_profiles_xvals)
        # for col_i in range(np.shape(filtered_gas_entropy_stacked)[1]):
        #     entropy_fig = plt.figure(figsize=(6,4))
        #     plt.plot( entropy_radial_profiles_xvals, np.log10(filtered_gas_entropy_stacked[:,col_i]), color = "green")
        #     plt.xlabel("Log R/R500")
        #     plt.ylabel("Log Sx")
        #     plt.xlim(-1.8, 0.5)
        #     plt.ylim(0.2, 6)
        #     plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/ENTROPY_PROFILES/zm{zoomnum}_Level{zoomlevel}_{str(snapnums_read[col_i]).zfill(3)}_xray_entropy_from_stacked_filtered_gas.png")
        #     plt.close(entropy_fig)
        





        
        # for ent_plot_i in range(2):
        #     if ent_plot_i == 1 and not PLOT_BHWIND_ENERGIES:
        #         break
        #     ent_plot_tag = ""
        #     fig, ax = plt.subplots()
        #     fig.set_size_inches((10, 5)) # 10 * nradiusbins/ntimebins))
            
        #     if ent_plot_i == 1:
        #         ent_plot_tag += f"_plus_jet_KE_coarsening_factor={KE_x_bins_coarsening}"
        #         fig = plt.figure()
        #         fig.set_size_inches((10, 8)) # 10 * nradiusbins/ntimebins))
        #         gs = fig.add_gridspec(2, hspace=0)
        #         ax, ax2 = gs.subplots(sharex=True, sharey=False,) #Note that if you do not use constrained layout, the pad command makes the parent Axes shrink:
        #         ax2.plot(KE_x_bin_midpoints, binned_KE)
        #         ax2.set_xlabel("dt (Myr)")
        #         ax2.set_ylabel("Jet KE (ergs)")
        #         ax.set_ylabel("log R/R500")
        #         plt.sca(ax)
                
        #     ent_cax = ax.inset_axes([1.01, 0, 0.05, 1])
        #     ax.set_facecolor('black')
                    
        #     entropy_plot_extent = [0, max(single_times_arr)-start_time, np.nanmin(np.log10(entropy_radial_profiles_xvals)), np.nanmax(np.log10(entropy_radial_profiles_xvals)) ]
        #     im = ax.imshow(np.log10(filtered_gas_entropy_stacked),  aspect='auto', vmin = 1, vmax = 2.2, extent = entropy_plot_extent,)
        #     cbar = fig.colorbar(im, label = "Log(Sx/keV cm**2)", cax = ent_cax)
            
        #     if ent_plot_i == 0: ax.set_xlabel("dt (Myr)")
        #     ax.set_ylabel("R/ R500")
        #     ax.set_ylim(entropy_plot_extent[2], -0.25)#entropy_plot_extent[3])
        #     plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_entropy_medians{ent_plot_tag}.png")
        #     plt.close()

        # print(entropy_radial_profiles_xvals)
        
        # print(np.count_nonzero(np.isnan(all_gas_entropy_stacked)))
        # # print(np.count_nonzero(np.isnan(all_gas_entropy_stacked)))
        # print(np.count_nonzero(np.isnan(np.log10(all_gas_entropy_stacked))))
        # exit()
        # os.makedirs(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/",exist_ok=True)
        
        for plot_j in range(2):
            for plot_i in range(16):
                histogram_plot_tag = ""
                
                if plot_j == 0:
                    fig, ax = plt.subplots()
                    fig.set_size_inches((10, 5)) # 10 * nradiusbins/ntimebins))
                    plt.xlabel("dt (Myr)")
                    plt.ylabel("log R/R500")
                    
                    
                
                if plot_j == 1:
                    if not PLOT_BHWIND_ENERGIES:
                        break
                    histogram_plot_tag += f"_plus_jet_KE_coarsening_factor={KE_x_bins_coarsening}"
                    fig = plt.figure()
                    fig.set_size_inches((10, 8)) # 10 * nradiusbins/ntimebins))
                    gs = fig.add_gridspec(2, hspace=0)
                    ax, ax2 = gs.subplots(sharex=True, sharey=False,) #Note that if you do not use constrained layout, the pad command makes the parent Axes shrink:


                    ax2.plot(KE_x_bin_midpoints, binned_KE)
                    ax2.set_xlabel("dt (Myr)")
                    ax2.set_ylabel("Jet KE (ergs)")
                    ax.set_ylabel("log R/R500")
                    plt.sca(ax)
                    
                # divider = make_axes_locatable(ax)
                # hist_cax = divider.append_axes('right', size='5%', pad=0.05)
                hist_cax = ax.inset_axes([1.01, 0, 0.05, 1])
                ax.set_facecolor('black')
                
                if plot_i == 0 and DO_INSTANTANEOUSLY_CONDENSED_PLOTS:
                    ntimebins = -1 + ( np.max(instantaneously_condensed_gas_all_times_Myr) - np.min(instantaneously_condensed_gas_all_times_Myr ))//DT
                    nradius_bins=30
                    h = ax.hist2d(x = instantaneously_condensed_gas_all_times_Myr, y = np.log10(instantaneous_all_radii/halo_r500), bins=(int(ntimebins), int(nradius_bins)), weights = instantaneous_all_masses, norm=matplotlib.colors.LogNorm())
                    fig.colorbar(h[3], label = "Condensing Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_instantaneous_{ntimebins}by{nradius_bins}_{histogram_plot_tag}.png")
                if plot_i == 1 and DO_INSTANTANEOUSLY_CONDENSED_PLOTS:
                    x_bins = np.linspace(min(instantaneously_condensed_gas_all_times_Myr)-(DT/2),max(instantaneously_condensed_gas_all_times_Myr)+(DT/2),len(set(instantaneously_condensed_gas_all_times_Myr)) + 1)
                    y_bins = np.linspace(-2.5,0.5,200)
                    h = ax.hist2d(x = instantaneously_condensed_gas_all_times_Myr, y = np.log10(instantaneous_all_radii/halo_r500), bins=(x_bins, y_bins), weights = instantaneous_all_masses, norm=matplotlib.colors.LogNorm())
                    fig.colorbar(h[3], label = "Condensing Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_instantaneous_custom_bins_{histogram_plot_tag}.png")
                    
                if plot_i == 15 and DO_INSTANTANEOUSLY_CONDENSED_PLOTS:
                    x_bins = np.linspace(min(instantaneously_condensed_gas_all_times_Myr)-(DT/2),max(instantaneously_condensed_gas_all_times_Myr)+(DT/2),len(set(instantaneously_condensed_gas_all_times_Myr)) + 1)
                    y_bins = np.linspace(-2.5,0.5,200)
                    h = ax.hist2d(x = instantaneously_condensed_gas_all_times_Myr, y = np.log10(instantaneous_all_radii/halo_r500), bins=(x_bins, y_bins), weights = instantaneous_all_masses, norm=matplotlib.colors.LogNorm(vmin = 1e6, vmax =1e8))
                    fig.colorbar(h[3], label = "Condensing Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_instantaneous_custom_bins_{histogram_plot_tag}_custom_vlims.png")
                    
                if plot_i == 2 and DO_INSTANTANEOUSLY_CONDENSED_PLOTS:
                    x_bins = np.linspace(min(instantaneously_condensed_gas_all_times_Myr)-(DT/2),max(instantaneously_condensed_gas_all_times_Myr)+(DT/2),int(0.5*len(set(instantaneously_condensed_gas_all_times_Myr))) + 1)
                    y_bins = np.linspace(-1.7,max(np.log10(recent_all_radii/halo_r500))+0.1,100)
                    h = ax.hist2d(x = instantaneously_condensed_gas_all_times_Myr, y = np.log10(instantaneous_all_radii/halo_r500), bins=(x_bins, y_bins), weights = instantaneous_all_masses, norm=matplotlib.colors.LogNorm())
                    fig.colorbar(h[3], label = "Condensing Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_instantaneous_custom_bins_beyond_0.01R500_{histogram_plot_tag}.png")    
                # if plot_i == 3:       
                #     h = ax.hist2d(x = recent_condensed_gas_all_times_Myr, y = np.log10(recent_all_radii/halo_r500), bins=(int(ntimebins), int(nradius_bins)), weights = recent_all_masses, norm=matplotlib.colors.LogNorm())
                #     fig.colorbar(h[3], ax=ax, label = "Condensing Mass")            
                #     plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_{ntimebins}by{nradius_bins}.png")
                if plot_i == 4:
                    x_bins = np.linspace(min(recent_condensed_gas_all_times_Myr)-(DT/2),max(recent_condensed_gas_all_times_Myr)+(DT/2),len(set(recent_condensed_gas_all_times_Myr)) + 1)
                    y_bins = np.linspace(-2.5,0.5,200)
                    h = ax.hist2d(x = recent_condensed_gas_all_times_Myr, y = np.log10(recent_all_radii/halo_r500), bins=(x_bins, y_bins), weights = recent_all_masses, norm=matplotlib.colors.LogNorm())
                    fig.colorbar(h[3], label = "Condensing Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_custom_bins_{histogram_plot_tag}.png")  
                if plot_i == 5:
                    x_bins = np.linspace(min(recent_condensed_gas_all_times_Myr)-(DT/2),max(recent_condensed_gas_all_times_Myr)+(DT/2),len(set(recent_condensed_gas_all_times_Myr)) + 1)
                    y_bins = np.linspace(-1.7,max(np.log10(recent_all_radii/halo_r500))+0.1,600)
                    h = ax.hist2d(x = recent_condensed_gas_all_times_Myr, y = np.log10(recent_all_radii/halo_r500), bins=(x_bins, y_bins), weights = recent_all_masses, norm=matplotlib.colors.LogNorm())
                    fig.colorbar(h[3], label = "Condensing Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_custom_bins_beyond_0.01R500_{histogram_plot_tag}.png")   
                if plot_i == 6:
                    x_bins = np.linspace(min(recent_condensed_gas_all_times_Myr)-(DT/2),max(recent_condensed_gas_all_times_Myr)+(DT/2),len(set(recent_condensed_gas_all_times_Myr)) + 1)
                    y_bins = np.linspace(-2.4,max(np.log10(recent_all_radii/halo_r500))+0.1,600)
                    h = ax.hist2d(x = recent_condensed_gas_all_times_Myr, y = np.log10(recent_all_radii/halo_r500), bins=(x_bins, y_bins), weights = recent_all_masses, norm=matplotlib.colors.LogNorm())
                    fig.colorbar(h[3],  label = "Condensing Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_custom_bins_beyond_-2.4_R500_{histogram_plot_tag}.png")
                # if plot_i == 7 and DO_CLOSE_TO_WAKE_AND_BUBBLES_PLOTS:       
                #     h = ax.hist2d(x = recent_condensed_gas_all_times_Myr_just_from_near_bubbles, y = np.log10(recent_condensed_gas_radii_just_from_near_bubbles/halo_r500), bins=(int(ntimebins), int(nradius_bins)), weights = recent_condensed_gas_masses_just_from_near_bubbles, norm=matplotlib.colors.LogNorm())
                #     fig.colorbar(h[3], ax=ax, label = "Condensing Mass")       
                #     plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_just_closetobubbles_condensed_gas_{ntimebins}by{nradius_bins}_{histogram_plot_tag}.png")
                if plot_i == 8 and DO_CLOSE_TO_WAKE_AND_BUBBLES_PLOTS:
                    x_bins = np.linspace(min(recent_condensed_gas_all_times_Myr_just_from_near_bubbles)-(DT/2),max(recent_condensed_gas_all_times_Myr_just_from_near_bubbles)+(DT/2),len(set(recent_condensed_gas_all_times_Myr_just_from_near_bubbles)) + 1)
                    y_bins = np.linspace(-2.5,0.5,200)
                    h = ax.hist2d(x = recent_condensed_gas_all_times_Myr_just_from_near_bubbles, y = np.log10(recent_condensed_gas_radii_just_from_near_bubbles/halo_r500), bins=(x_bins, y_bins), weights = recent_condensed_gas_masses_just_from_near_bubbles, norm=matplotlib.colors.LogNorm())
                    fig.colorbar(h[3], label = "Condensing Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_just_closetobubbles_condensed_gas_custom_bins_{histogram_plot_tag}.png")  
                if plot_i == 9 and DO_CLOSE_TO_WAKE_AND_BUBBLES_PLOTS:
                    x_bins = np.linspace(min(recent_condensed_gas_all_times_Myr_just_from_near_bubbles)-(DT/2),max(recent_condensed_gas_all_times_Myr_just_from_near_bubbles)+(DT/2),len(set(recent_condensed_gas_all_times_Myr_just_from_near_bubbles)) + 1)
                    y_bins = np.linspace(-1.7,-0.35,600)
                    h = ax.hist2d(x = recent_condensed_gas_all_times_Myr_just_from_near_bubbles, y = np.log10(recent_condensed_gas_radii_just_from_near_bubbles/halo_r500), bins=(x_bins, y_bins), weights = recent_condensed_gas_masses_just_from_near_bubbles, norm=matplotlib.colors.LogNorm())
                    fig.colorbar(h[3], label = "Condensing Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_just_closetobubbles_condensed_gas_custom_bins_beyond_0.01R500_{histogram_plot_tag}.png")   
                if plot_i == 10 and DO_CLOSE_TO_WAKE_AND_BUBBLES_PLOTS:
                    x_bins = np.linspace(min(recent_condensed_gas_all_times_Myr_just_from_near_bubbles)-(DT/2),max(recent_condensed_gas_all_times_Myr_just_from_near_bubbles)+(DT/2),len(set(recent_condensed_gas_all_times_Myr_just_from_near_bubbles)) + 1)
                    y_bins = np.linspace(-2.4,-0.35,600)
                    h = ax.hist2d(x = recent_condensed_gas_all_times_Myr_just_from_near_bubbles, y = np.log10(recent_condensed_gas_radii_just_from_near_bubbles/halo_r500), bins=(x_bins, y_bins), weights = recent_condensed_gas_masses_just_from_near_bubbles, norm=matplotlib.colors.LogNorm())
                    fig.colorbar(h[3], label = "Condensing Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_just_closetobubbles_condensed_gas_custom_bins_beyond_-2.4_R500_{histogram_plot_tag}.png")
                if plot_i == 11 and DO_CLOSE_TO_WAKE_AND_BUBBLES_PLOTS:
                    x_bins = np.linspace(min(recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake)-(DT/2),max(recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake)+(DT/2),len(set(recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake)) + 1)
                    y_bins = np.linspace(-2.4,-0.35,600)
                    h = ax.hist2d(x = recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake, y = np.log10(recent_condensed_gas_radii_just_from_near_satellite_wake/halo_r500), bins=(x_bins, y_bins), weights = recent_condensed_gas_masses_just_from_near_satellite_wake, norm=matplotlib.colors.LogNorm())
                    fig.colorbar(h[3], label = "Condensing Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_just_closetowake_condensed_gas_custom_bins_beyond_-2.4_R500_{histogram_plot_tag}.png")
                if plot_i == 12:
                    x_bins = np.linspace(min(recent_still_hot_bubble_particle_times_arr)-(DT/2),max(recent_still_hot_bubble_particle_times_arr)+(DT/2),len(set(recent_still_hot_bubble_particle_times_arr)) + 1)
                    y_bins = np.linspace(-2.4,0.5,600)
                    h = ax.hist2d(x = recent_still_hot_bubble_particle_times_arr, y = np.log10(recent_still_hot_bubble_particle_radii_arr/halo_r500), bins=(x_bins, y_bins), weights = recent_still_hot_bubble_particle_masses_arr, norm=matplotlib.colors.LogNorm(), cmap = "cmr.rainforest")
                    fig.colorbar(h[3], label = "Hot Bubble Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_mass_in_hot_bubbles_custom_bins_beyond_-2.4_R500_{histogram_plot_tag}.png")  
                if plot_i == 13:
                    x_bins = np.linspace(min(recent_cooled_bubble_particle_times_arr)-(DT/2),max(recent_cooled_bubble_particle_times_arr)+(DT/2),len(set(recent_cooled_bubble_particle_times_arr)) + 1)
                    y_bins = np.linspace(-2.4,0.5,600)
                    h = ax.hist2d(x = recent_cooled_bubble_particle_times_arr, y = np.log10(recent_cooled_bubble_particle_radii_arr/halo_r500), bins=(x_bins, y_bins), weights = recent_cooled_bubble_particle_masses_arr, norm=matplotlib.colors.LogNorm(), cmap = "Wistia_r")
                    fig.colorbar(h[3], label = "Cooled Bubble Mass", cax = hist_cax)            
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_mass_in_cooled_bubbles_custom_bins_beyond_-2.4_R500_{histogram_plot_tag}.png")
                    
                    
                if plot_i == 14 and DO_CLOSE_TO_WAKE_AND_BUBBLES_PLOTS:
                    x_bins = np.linspace(min(recent_condensed_gas_all_times_Myr)-(DT/2),max(recent_condensed_gas_all_times_Myr)+(DT/2),len(set(recent_condensed_gas_all_times_Myr)) + 1)
                    y_bins = np.linspace(-2.4,0.5,600)
                    all_condensation_hist ,_ , _= np.histogram2d(x = recent_condensed_gas_all_times_Myr, y = np.log10(recent_all_radii/halo_r500), bins=(x_bins, y_bins), weights = recent_all_masses,)
                    satwake_hist,_ , _ = np.histogram2d(x = recent_condensed_gas_all_times_Myr_just_from_near_satellite_wake, y = np.log10(recent_condensed_gas_radii_just_from_near_satellite_wake/halo_r500), bins=(x_bins, y_bins), weights = recent_condensed_gas_masses_just_from_near_satellite_wake,)
                    satinfront_hist,_ , _ = np.histogram2d(x = recent_condensed_gas_all_times_Myr_just_from_near_satellite_infront, y = np.log10(recent_condensed_gas_radii_just_from_near_satellite_infront/halo_r500), bins=(x_bins, y_bins), weights = recent_condensed_gas_masses_just_from_near_satellite_infront,)
                    nearbubbles_hist ,xedges, yedges= np.histogram2d(x = recent_condensed_gas_all_times_Myr_just_from_near_bubbles, y = np.log10(recent_condensed_gas_radii_just_from_near_bubbles/halo_r500), bins=(x_bins, y_bins), weights = recent_condensed_gas_masses_just_from_near_bubbles,)
                    unprovoked_condensation_hist = all_condensation_hist - (satwake_hist+satinfront_hist+nearbubbles_hist)
                    ax.imshow(unprovoked_condensation_hist, interpolation='nearest', origin='lower',
                            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])  
                    plt.savefig(f"../DATA_AND_PLOTS/{DATE_STR}/CONDENSED_GAS_and_BUBBLES_FINAL_DATA_PLOTS/zm{zoomnum}_Level{zoomlevel}_2dhist_recent_unprovoked_condensed_mass_in_cooled_bubbles_custom_bins_beyond_-2.4_R500_{histogram_plot_tag}.png")
                    
                    
                    # cooled_bubble_hist= np.histogram2d(x = recent_cooled_bubble_particle_times_arr, y = np.log10(recent_cooled_bubble_particle_radii_arr/halo_r500), bins=(x_bins, y_bins), weights = recent_cooled_bubble_particle_masses_arr, )

                    
                    
                    
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
            plt.close()
            
        # 

exit()
