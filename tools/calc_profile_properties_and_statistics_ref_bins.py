## Import libraries
import numpy as np
import pynbody as pnb
import os
import copy
import argparse
import dill
import pprint

from distutils.util import strtobool

import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

from gen_sim_data import weighted_mean_of_profiles, weighted_quantiles_of_profiles, calc_bin_overlap

from astro_constants import NA_no_units, NA, kB_no_units, kB, mH, mp_no_units, mp, mu_e, mu, G_no_units, G
from astro_constants import ABOUT_ZERO, SMALL_OVERLAP
from solar_abundances import info_Asplund09, info_Lodders09, Z_tot_Asplund09, Z_tot_Lodders09



## Command line keyword arguments
parser = argparse.ArgumentParser(description="Analyze the X-ray profiles of the IGrM.")
parser.add_argument('--snap', action='store', type=str, required=True,
                   help='Name of snapshot, eg. snapshot_151')
parser.add_argument('--dir', action='store', type=str, required=True, 
                   help='Path to profiles file to analyze')
parser.add_argument('--name', action='store', type=str, required=True, 
                   help='Name of profiles file')
parser.add_argument('--suffix', action='store', type=str, required=False, 
                   help='keywords to add to end of saved file names to make them recognizable')

parser.add_argument('--code', action='store', type=str, required=True, 
                   help='Simulation type (currently either just Simba or Simba-C')
parser.add_argument('--size', action='store', type=str, required=True, 
                   help='Simulation size, eg. N512L50 or N1024L100')

parser.add_argument('--ndim', action='store', type=int, required=True, 
                   help='Number of dimensions for profiles (2 or 3)')
parser.add_argument('--change_xscale', type=lambda x:bool(strtobool(x)), required=False, default=False, #nargs='?',const=True, 
                    help='Change scaling of radial axis (default=False). If changing xscale, \
                   note that radial bins will likely become different for each halo, so averaging will require \
                   you to set use_ref_bins=True.')
parser.add_argument('--new_xscale', action='store', type=str, required=False, default=None,
                   help='New units of radial axis (R500, R200, Physical)')

parser.add_argument('--avg_type', action='store', type=str, required=False, default='mean',
                   help='If using reference radial bins, must decide whether to employ weighted mean and standard deviation \
                   or weighted median and 16th/84th quantiles (default=mean, median)')
parser.add_argument('--use_errors', type=lambda x:bool(strtobool(x)), required=False, default=False, #nargs='?',const=True, 
                    help='If using reference radial bins, must decide whether to weight average by errors (default=False)')

parser.add_argument('--res_file', action='store', type=str, required=True, 
                   help='Name of profiles file to use for resolution cuts (if profiles you are analyzing are 3D then this should\
                   be the same as "name", if profiles are 2D, then this should be the corresponding 3D profile (must have same\
                   snapshot number, directory and xscale as "name", and be the same code and same size as "name"))')
parser.add_argument('--profiles_are_2d_but_res_file_is_3d', type=lambda x:bool(strtobool(x)), required=True,
                    help='Set True if "name" contains 2D profiles, but "res_file" contains 3D profiles')
parser.add_argument('--use_3d_for_2d_res', type=lambda x:bool(strtobool(x)), required=True,
                    help='Set True if resolution cuts of 2D profiles should use 3D profiles')

parser.add_argument('--three_d_file', action='store', type=str, required=False, 
                   help='Name of 3D profiles file to use for CC/NCC classifications (must have same\
                   snapshot number, directory and xscale as "name", and be the same code and same size as "name")')
parser.add_argument('--two_d_file', action='store', type=str, required=False, 
                   help='Name of 2D profiles file to use for CC/NCC classifications (must have same\
                   snapshot number, directory and xscale as "name", and be the same code and same size as "name")')

parser.add_argument('--include_resolution_cuts', type=lambda x:bool(strtobool(x)), required=False, default=False,
                    help='Include cuts that remove unresolved bins of individual halos from averaging (default=False)')
parser.add_argument('--min_resolved_counts_per_halo', action='store', type=int, required=False, default=64,
                   help='Minimum number of particles in a radial bin to be considered resolved \
                   (sets bins of halos not meeting this criteria to NaN; default=64)')
parser.add_argument('--min_resolved_counts_per_halo_core', action='store', type=int, required=False, default=16,
                   help='Minimum number of particles in core radial bin to be considered resolved \
                   (sets bins of halos not meeting this criteria to NaN; default=16)')
parser.add_argument('--min_resolved_counts_per_bin', action='store', type=int, required=False, default=64,
                   help='Minimum number of particles in a radial bin to be considered resolved, \
                   specifically for excluding bins with too few resolved halos (default=64)')
parser.add_argument('--min_total_halos', action='store', type=int, required=False, default=5,
                   help='Minimum number of total halos in a radial bin to be considered sufficient for \
                   averaging over a population (default=5)')
parser.add_argument('--min_resolved_halos', action='store', type=int, required=False, default=5,
                   help='Minimum number of resolved halos in a radial bin to be considered sufficient for \
                   averaging over a population (default=5)')
parser.add_argument('--min_frac_resolved_halos', action='store', type=float, required=False, default=0.16,
                   help='Minimum fraction of resolved halos in a radial bin (num_resolved_halos/num_total_halos) \
                   to be considered sufficient for averaging over a population (default=0.16)')

## Choose which extra properties to compute
parser.add_argument('--calc_thermo_props', type=lambda x:bool(strtobool(x)), required=False, default=False,
                    help='Calculate additional thermal profiles for each halo (default=False)')
parser.add_argument('--calc_metal_props', type=lambda x:bool(strtobool(x)), required=False, default=False,
                    help='Calculate additional abundance/metallicity profiles for each halo (default=False)')
parser.add_argument('--classify_halos', type=lambda x:bool(strtobool(x)), required=False, default=False,
                    help='Classify each halo, eg. CC vs NCC (default=False)')
parser.add_argument('--calc_log_props', type=lambda x:bool(strtobool(x)), required=False, default=False,
                    help='Calculate logarithmic profiles for each halo (default=False)')

args = parser.parse_args()


## Save input parameters to output with results
parameters = {}
for arg, arg_val in vars(args).items():
    parameters[arg] = arg_val


## Input parameters as global variables
SNAP = args.snap
PROFILES_DIR = args.dir
PROFILES_SAVE_NAME = args.name
SUFFIX = args.suffix

CODE = args.code
SIZE = args.size

NDIM = args.ndim
CHANGE_XSCALE = args.change_xscale
NEW_XSCALE = args.new_xscale

AVG_TYPE = args.avg_type
USE_ERRORS = args.use_errors

RES_FILE = args.res_file
PROFILES_ARE_2D_BUT_RES_FILE_IS_3D = args.profiles_are_2d_but_res_file_is_3d
USE_3D_FOR_2D_RES = args.use_3d_for_2d_res

THREE_D_FILE = args.three_d_file
TWO_D_FILE = args.two_d_file

INCLUDE_RESOLUTION_CUTS = args.include_resolution_cuts
MIN_RESOLVED_COUNTS_PER_HALO = args.min_resolved_counts_per_halo
MIN_RESOLVED_COUNTS_PER_HALO_CORE = args.min_resolved_counts_per_halo_core
MIN_RESOLVED_COUNTS_PER_BIN = args.min_resolved_counts_per_bin
MIN_TOTAL_HALOS = args.min_total_halos
MIN_RESOLVED_HALOS = args.min_resolved_halos
MIN_FRAC_RESOLVED_HALOS = args.min_frac_resolved_halos

CALC_THERMO_PROPS = args.calc_thermo_props
CALC_METAL_PROPS = args.calc_metal_props
CLASSIFY_HALOS = args.classify_halos
CALC_LOG_PROPS = args.calc_log_props

# if CHANGE_XSCALE:
#     SUFFIX = 'new_xscale=%s-%s' % (NEW_XSCALE, SUFFIX)
    
# SUFFIX = '%s-use_errors=%s-min_resolved_counts_per_halo=%s-min_resolved_counts_per_bin=%s-min_total_halos=%s-min_resolved_halos=%s-min_frac_resolved_halos=%s-%s' % (AVG_TYPE, USE_ERRORS, SUFFIX, MIN_RESOLVED_COUNTS_PER_HALO, MIN_RESOLVED_COUNTS_PER_BIN, MIN_TOTAL_HALOS, MIN_RESOLVED_HALOS, MIN_FRAC_RESOLVED_HALOS)



## Axis of all arrays of all properties along which the radii increase
RADIAL_AXIS = 0



## Set base units and units of x-axis depending on scaling
R_UNITS = 'kpc'
M_UNITS = 'Msol'

XSCALE_UNITS_DICT = {
    'R500':'',
    'R200':'',
    'Physical':R_UNITS,
}


## Simulation cosmological properties
COSMOLOGY_SIM = FlatLambdaCDM(name=r'Flat $\Lambda$CDM', H0=68 * u.km / u.s / u.Mpc, Om0=0.3)
h_sim = COSMOLOGY_SIM.h
# h70_sim = h_sim * 100 / 70

if '100' in SIZE:
    GRAV_SOFTENING_LENGTH = 0.5
elif '50' in SIZE:
    GRAV_SOFTENING_LENGTH = 0.25
elif '25' in SIZE:
    GRAV_SOFTENING_LENGTH = 0.125
elif '12.5' in SIZE:
    GRAV_SOFTENING_LENGTH = 0.0625

GRAV_SOFTENING_LENGTH *= 2.8/h_sim  ## 2.8 for assumed cubic spline kernel conversion from 
                                    ## Plummer equivalent to actual smoothing kernel radius
GRAV_SOFTENING_LENGTH = pnb.array.SimArray(GRAV_SOFTENING_LENGTH, units='kpc')



## Indexes of different metals in Simba and Simba-C
metals_idx_simbac = {'H': 0, 'He': 1, 'Li': 2, 'Be': 3, 'B': 4, 'C': 5, 'N': 6, 'O': 7, 'F': 8, 'Ne': 9, 'Na': 10, 'Mg': 11, 'Al':12,
              'Si': 13, 'P': 14, 'S': 15, 'Cl': 16, 'Ar': 17, 'K': 18,'Ca': 19, 'Sc': 20, 'Ti': 21, 'V': 22, 'Cr': 23,
              'Mn': 24, 'Fe': 25, 'Co': 26, 'Ni': 27, 'Cu': 28, 'Zn': 29}
metals_idx_simba = {'H':0, 'He': 1, 'C': 2, 'N': 3, 'O': 4, 'Ne': 5, 'Mg': 6, 'Si': 7, 'S': 8, 'Ca': 9, 'Fe': 10}
metals_idx = {'Simba':metals_idx_simba, 'Simba-C':metals_idx_simbac}



## For getting correct units of scaling properties in 2d vs 3d
geometric_units = {
    ## Num dimensions
    2:{
        'volume':'cm**2',
        '1/volume':'cm**-2',
    },
    3:{
        'volume':'cm**3',
        '1/volume':'cm**-3',
    },
}


    
## Load profiles file
def save_object_with_dill(obj, filename, protocol=dill.HIGHEST_PROTOCOL):
    with open(filename, 'wb') as f:  # Overwrites any existing file.
        dill.dump(obj, f, protocol)

profiles_file = os.path.join(PROFILES_DIR, SNAP + '_' + PROFILES_SAVE_NAME + '.pkl')
with open(profiles_file, 'rb') as f:
    profiles_ = dill.load(f)
    
res_file = os.path.join(PROFILES_DIR, SNAP + '_' + RES_FILE + '.pkl')
with open(res_file, 'rb') as f:
    res_profiles_ = dill.load(f) ## don't need to calculate anything extra, since these are just for resolution cuts
    
three_d_file = os.path.join(PROFILES_DIR, SNAP + '_' + THREE_D_FILE + '.pkl')
with open(three_d_file, 'rb') as f:
    three_d_profiles_ = dill.load(f)
    
two_d_file = os.path.join(PROFILES_DIR, SNAP + '_' + TWO_D_FILE + '.pkl')
with open(two_d_file, 'rb') as f:
    two_d_profiles_ = dill.load(f)
    
    
## Save original x-axis if changing its scaling
profiles_['xaxis_og'] = copy.deepcopy(profiles_['xaxis'])
res_profiles_['xaxis_og'] = copy.deepcopy(res_profiles_['xaxis'])
three_d_profiles_['xaxis_og'] = copy.deepcopy(three_d_profiles_['xaxis'])
two_d_profiles_['xaxis_og'] = copy.deepcopy(two_d_profiles_['xaxis'])
print()
print('Original:')
print(profiles_['xaxis']['all'])
print()

## Original scaling of xaxis
XSCALE = profiles_['xaxis']['all']['scale']
if CHANGE_XSCALE:
    XSCALE = NEW_XSCALE

## Change scaling of x-axis to new one
if XSCALE.lower() == 'r500':
    # Intermediate Bins for R500
    _rbin1 = np.append([0], np.logspace(-2, 0, 10))
    _rbin2 = np.logspace(0, np.log10(3), 5)
    _rbin3 = np.logspace(np.log10(3), 1, 10)
    _rbins = np.unique(np.concatenate((_rbin1, _rbin2, _rbin3)))
elif XSCALE.lower() == 'physical':
    # Intermediate Bins for physical radius
    _rbin1 = np.append([0], np.logspace(0.5, 3, 20))
    _rbin2 = np.logspace(3, np.log10(3000), 5)
    _rbin3 = np.logspace(np.log10(3000), 4, 10)
    _rbins = np.unique(np.concatenate((_rbin1, _rbin2, _rbin3)))
else:
    print('ERROR: XSCALE IS NOT ONE OF THE VALID OPTIONS')
    raise SystemExit

REF_BINS = pnb.array.SimArray(_rbins, units=XSCALE_UNITS_DICT[XSCALE])
log_bins = np.log10(REF_BINS)
bin_centres = (REF_BINS[:-1] + REF_BINS[1:])/2.
log_bin_centres = (log_bins[:-1] + log_bins[1:])/2.

profiles_['xaxis']['all']['scale'] = XSCALE
profiles_['xaxis']['all']['bins'] = REF_BINS
profiles_['xaxis']['all']['log_bins'] = log_bins
profiles_['xaxis']['all']['bin_centres'] = bin_centres
profiles_['xaxis']['all']['log_bin_centres'] = log_bin_centres

res_profiles_['xaxis']['all']['scale'] = XSCALE
res_profiles_['xaxis']['all']['bins'] = REF_BINS
res_profiles_['xaxis']['all']['log_bins'] = log_bins
res_profiles_['xaxis']['all']['bin_centres'] = bin_centres
res_profiles_['xaxis']['all']['log_bin_centres'] = log_bin_centres

three_d_profiles_['xaxis']['all']['scale'] = XSCALE
three_d_profiles_['xaxis']['all']['bins'] = REF_BINS
three_d_profiles_['xaxis']['all']['log_bins'] = log_bins
three_d_profiles_['xaxis']['all']['bin_centres'] = bin_centres
three_d_profiles_['xaxis']['all']['log_bin_centres'] = log_bin_centres

two_d_profiles_['xaxis']['all']['scale'] = XSCALE
two_d_profiles_['xaxis']['all']['bins'] = REF_BINS
two_d_profiles_['xaxis']['all']['log_bins'] = log_bins
two_d_profiles_['xaxis']['all']['bin_centres'] = bin_centres
two_d_profiles_['xaxis']['all']['log_bin_centres'] = log_bin_centres

## Since the core profiles are just being used to classify the groups as CC or NCC,
## it doesn't matter if they get averaged correctly or not

print()
print('New:')
print(profiles_['xaxis']['all'])
print()




## deepcopy of profiles_ to make sure profiles_ remains intact for saving later
profiles = copy.deepcopy(profiles_)
res_profiles = copy.deepcopy(res_profiles_)
three_d_profiles = copy.deepcopy(three_d_profiles_)
two_d_profiles = copy.deepcopy(two_d_profiles_)

## All x and y properties
print()
print()
print(profiles['all_profiles'][profiles_['all_id'][0]]['x'].keys())
print()
print(profiles['all_profiles'][profiles_['all_id'][0]]['y'].keys())
print()
print()



## Xray bands to use
# xray_bands = ['0.5-1.8keV_total', '0.5-2.0keV_total', '0.1-2.4keV_total', '0.5-7.0keV_total', '0.5-10.0keV_total', 
#               '0.5-1.8keV_cont', '0.5-2.0keV_cont', '0.1-2.4keV_cont', '0.5-7.0keV_cont', '0.5-10.0keV_cont']
# xray_bands = ['0.5-10.0keV_total', '0.5-10.0keV_cont']
xray_bands = ['0.5-10.0keV_total']


    
## Classical Poisson error to use when averaging
def calc_poisson_err(x, log=False):
    if log:
        return (10**x)**0.5
    else:
        return x**0.5


def set_nonfinite_to_nan(values):
    values[~np.isfinite(values)] = np.nan
    
    
    



## Properties to remove to save memory and storage space

x_props_to_remove = []


# weighted_prop_weights_og = ['median', 'mass', 'volume', 'ne', 'T'] + ['Lx_' + xband for xband in xray_bands]
weighted_prop_weights_og = ['median', 'mass', 'volume'] + ['Lx_' + xband for xband in xray_bands]

if CALC_THERMO_PROPS:
    y_props_to_remove = []  ## good (thermo, all)
else:
#     y_props_to_remove = ['volume', 'Ne', 'Nh', 'mass', 'T_spec'] + ['Lx_' + xband for xband in xray_bands] ## good (metallicity)
    y_props_to_remove = ['volume', 'Ne', 'Nh']

if CALC_THERMO_PROPS and CALC_METAL_PROPS:
    weighted_props_to_remove = []  ## good (all)
elif CALC_THERMO_PROPS:
    weighted_props_to_remove = ['Z', 'Ztot']  ## good (thermo)
elif CALC_METAL_PROPS:
    weighted_props_to_remove = ['T', 'ne', 'K', 'P'] + ['tcool_' + xband for xband in xray_bands] ## good (metallicity)
#     weighted_props_to_remove = []
else:
    weighted_props_to_remove = ['Z', 'Ztot', 'T', 'ne', 'K', 'P'] + ['tcool_' + xband for xband in xray_bands] ## good (none)


for weighted_prop in weighted_props_to_remove:
    for weight in weighted_prop_weights_og:
        if weighted_prop == 'Ztot' and weight == 'median':
            continue
        if weighted_prop != weight:
            y_props_to_remove.append(weighted_prop + '_' + weight)

# print()
# print()
# print(y_props_to_remove)
# print()
# print()




## Properties to weight properly, sum properly, and metallicity types to compute

# weighted_props = ['Z', 'nX/nH', 'T', 'ne', 'K', 'P'] + ['tcool_' + xband for xband in xray_bands]
weighted_prop_weights = ['median', 'mass', 'volume'] + ['Lx_' + xband for xband in xray_bands]
# weighted_prop_weights = ['mass', 'volume'] + ['Lx_' + xband for xband in xray_bands]

sum_props = ['n', 'Ne', 'Nh', 'mass'] + ['Lx_' + xband for xband in xray_bands]
mean_props = ['mass'] + ['Lx_' + xband for xband in xray_bands]

# metallicity_types = ['Z', 'nX/nH']
metallicity_types = ['Z']


## dict for scaling thermal properties by characteristic values, eg. M500
thermo_scaling_dict = {
    'mass':{
        'scale_prop':'M',
        'scale_types':['500', 'igrm500'],# 'bar500', 'star500', 'self_star', 'total_star'],
    },
    'mass_enc':{
        'scale_prop':'M',
        'scale_types':['500', 'igrm500'],# 'bar500', 'star500', 'self_star', 'total_star'],
    },
}

for T_weight in ['median', 'spec'] + weighted_prop_weights:
    if T_weight != 'T':
        thermo_scaling_dict['T_' + T_weight] = {
            'scale_prop':'T',
            'scale_types':['500ideal']
        }
        
for kT_weight in ['median', 'spec'] + weighted_prop_weights:
    if kT_weight != 'T':
        thermo_scaling_dict['kT_' + kT_weight] = {
            'scale_prop':'T',
            'scale_types':['spec_corr', 'spec500']
        }
        
# Can't use this for 2d profiles since ne500 will always be in units of cm**-3
# Anything that has ne in it
for ne_weight in ['bin', 'median'] + weighted_prop_weights:
    if ne_weight not in ['rho', 'ne', 'nh']:
        thermo_scaling_dict['ne_' + ne_weight] = {
            'scale_prop':'ne',
            'scale_types':['500ideal']
        }
        
# for K_weight in ['median'] + weighted_prop_weights + [weight_+'_bin_T' for weight_ in ['median', 'spec'] + weighted_prop_weights] + [weight_+'_bin_T_ne' for weight_ in ['median', 'spec'] + weighted_prop_weights]:
for K_weight in ['median'] + weighted_prop_weights:
    if K_weight not in ['T_bin_T', 'T_bin_T_ne', 'spec_bin_T_ne', 'ne_bin_T_ne']:
        thermo_scaling_dict['K_' + K_weight] = {
            'scale_prop':'S',
            'scale_types':['500ideal', '500']
        }
        
# for P_weight in ['median'] + weighted_prop_weights + [weight_+'_bin_T' for weight_ in ['median', 'spec'] + weighted_prop_weights] + [weight_+'_bin_T_ne' for weight_ in ['median', 'spec'] + weighted_prop_weights]:
for P_weight in ['median'] + weighted_prop_weights:
    if P_weight not in ['T_bin_T', 'T_bin_T_ne', 'spec_bin_T_ne', 'ne_bin_T_ne']:
        thermo_scaling_dict['P_' + P_weight] = {
            'scale_prop':'P',
            'scale_types':['500ideal']
        }

        
        

# dict for scaling metallicities by solar abundances
metallicity_scaling_dict = {
    'asplund09':{
        'info':info_Asplund09,
        'Z_name':'ZX',
        'nX/nH_name':'NX_over_NH',
        'total':Z_tot_Asplund09,
    },
#     'lodders09':{
#         'info':info_Lodders09,
#         'Z_name':'ZX',
#         'nX/nH_name':'NX_over_NH',
#         'total':Z_tot_Lodders09,
#     },
}





## Compute all extra properties for each halo
print()
print('Computing extra profile properties')
print()

## need to make sure 3d and 2d profile values calculated before profiles for CC/NCC classification
profiles_to_compute_names = ['three_d_profiles', 'two_d_profiles', 'profiles', 'res_profiles']
profiles_to_compute = [three_d_profiles, two_d_profiles, profiles, res_profiles]
profiles_ndims = [3, 2, NDIM, 3 if PROFILES_ARE_2D_BUT_RES_FILE_IS_3D else 2] ## last one won't always work
compute_properties = [True, True, True, False]
compute_metallicities = [False, False, True, False]
# profiles_to_compute = ['all_profiles', 'all_core_profiles']
profile_types_to_compute = ['all_core_profiles', 'all_profiles']
profile_ids_to_compute = ['all_core_id', 'all_id']

for profile_name, profile_to_compute, profile_ndim, compute_props, compute_metallicity in zip(profiles_to_compute_names, profiles_to_compute, profiles_ndims, compute_properties, compute_metallicities):
    print()
    print()
    print(profile_name)
    print()
#     for halo in profiles['all_id']:
    for profile_type_to_compute, profile_id_to_compute in zip(profile_types_to_compute, profile_ids_to_compute):
        print()
        print()
        print(profile_type_to_compute)
        print()
#     for halo in profile_to_compute['all_id']:
        for halo in profile_to_compute[profile_id_to_compute]:
            print('halo', halo)
    #         halo_idx = list(profiles['temp_halo_info']['halo_id']).index(halo)
            halo_idx = list(profile_to_compute['temp_halo_info']['halo_id']).index(halo)
            print('halo_idx:', halo_idx)
#             print()


#         for profile_type_to_compute in profile_types_to_compute:
            try:
#                 x_profile = profiles[profile_to_compute][halo]['x']
#                 x_res_profile = res_profiles[profile_to_compute][halo]['x']
#                 profile = profiles[profile_to_compute][halo]['y']
    
                x_profile = profile_to_compute[profile_type_to_compute][halo]['x']
                profile = profile_to_compute[profile_type_to_compute][halo]['y']
            except:
                print(profile_type_to_compute, 'does not exist for halo', halo)
                continue


            ## Do x quantity calculations
            xscale_values = {
#                 'R500':pnb.array.SimArray(profiles['temp_halo_info']['R']['500'].in_units(R_UNITS)[halo_idx], units=R_UNITS),
#                 'R200':pnb.array.SimArray(profiles['temp_halo_info']['R']['200'].in_units(R_UNITS)[halo_idx], units=R_UNITS),
                'R500':pnb.array.SimArray(profile_to_compute['temp_halo_info']['R']['500'].in_units(R_UNITS)[halo_idx], units=R_UNITS),
                'R200':pnb.array.SimArray(profile_to_compute['temp_halo_info']['R']['200'].in_units(R_UNITS)[halo_idx], units=R_UNITS),
                'Physical':pnb.array.SimArray(1, units=''),
            }

            ## I made it so that the radial bins that come with the individual halos are always physical (but this could be different!)
            x_profile['physical_rbins'] = copy.deepcopy(x_profile['rbins'])# * old_xscale_value/xscale_values['Physical']
#             x_res_profile['physical_rbins'] = copy.deepcopy(x_res_profile['rbins'])

            xscale_value = xscale_values[XSCALE]
            print('XSCALE:', XSCALE, ':', xscale_value, xscale_value.units)
            x_profile['rbins'] *= xscale_values['Physical']/xscale_value
#             x_res_profile['rbins'] *= xscale_values['Physical']/xscale_value

            x_profile['rbins_r500'] = x_profile['physical_rbins'].in_units(R_UNITS)/xscale_values['R500']


            ## don't think I need the below info ##
            x_profile['log_rbins'] = np.log10(x_profile['rbins'])
            x_profile['bin_centres'] = (x_profile['rbins'][:-1] + x_profile['rbins'][1:])/2.
            x_profile['log_bin_centres'] = (x_profile['log_rbins'][:-1] + x_profile['log_rbins'][1:])/2.



            ## Do y quantity calculations

            # Give number of particles a unit
            profile['n'] = pnb.array.SimArray(profile['n'].astype(float), units='')
            
            if compute_props:

                ## Always calculate some basic quantities
                profile['volume_enc'] = np.nancumsum(profile['volume'])

#                 print('volume units:', profile['volume'].units)

                profile['ne_bin'] = (profile['Ne']/profile['volume']).in_units(geometric_units[profile_ndim]['1/volume'])
                profile['nh_bin'] = (profile['Nh']/profile['volume']).in_units(geometric_units[profile_ndim]['1/volume'])
                profile['fe'] = profile['Ne']/profile['Nh']

                for sum_prop in sum_props:
                    profile[sum_prop + '_enc'] = np.nancumsum(profile[sum_prop], axis=RADIAL_AXIS)
                    profile[sum_prop + '_density'] = profile[sum_prop] / profile['volume']
                    profile[sum_prop + '_density_enc'] = profile[sum_prop + '_enc'] / profile['volume_enc']

                    ## Add extra axis so can be plotted with other metals as if it is an element, by giving its element as H (ie. index 0)
                    profile[sum_prop + '_for_metals'] = profile[sum_prop][:,np.newaxis]
                    profile[sum_prop + '_enc_for_metals'] = profile[sum_prop + '_enc'][:,np.newaxis]
                    profile[sum_prop + '_density_for_metals'] = profile[sum_prop + '_density'][:,np.newaxis]
                    profile[sum_prop + '_density_enc_for_metals'] = profile[sum_prop + '_density_enc'][:,np.newaxis]


    #             profile['mass_density_fb_rhocrit'] = profile['mass_density'] / profiles['cosmo_props']['rho_gas_scaling']
    #             profile['mass_density_enc_fb_rhocrit'] = profile['mass_density_enc'] / profiles['cosmo_props']['rho_gas_scaling']
    #             profile['mass_density_rhocrit'] = profile['mass_density'] / profiles['cosmo_props']['rho_crit']
    #             profile['mass_density_enc_rhocrit'] = profile['mass_density_enc'] / profiles['cosmo_props']['rho_crit']

                profile['mass_density_fb_rhocrit'] = profile['mass_density'] / profile_to_compute['cosmo_props']['rho_gas_scaling']
                profile['mass_density_enc_fb_rhocrit'] = profile['mass_density_enc'] / profile_to_compute['cosmo_props']['rho_gas_scaling']
                profile['mass_density_rhocrit'] = profile['mass_density'] / profile_to_compute['cosmo_props']['rho_crit']
                profile['mass_density_enc_rhocrit'] = profile['mass_density_enc'] / profile_to_compute['cosmo_props']['rho_crit']

                ## Add extra axis so can be plotted with other metals as if it is an element, by giving its element as H (ie. index 0)
                profile['mass_density_fb_rhocrit_for_metals'] = profile['mass_density_fb_rhocrit'][:,np.newaxis]
                profile['mass_density_enc_fb_rhocrit_for_metals'] = profile['mass_density_enc_fb_rhocrit'][:,np.newaxis]
                profile['mass_density_rhocrit_for_metals'] = profile['mass_density_rhocrit'][:,np.newaxis]
                profile['mass_density_enc_rhocrit_for_metals'] = profile['mass_density_enc_rhocrit'][:,np.newaxis]
                
                
                
                for mean_prop in mean_props:
                    profile[mean_prop + '_mean'] = profile[mean_prop]/profile['n']


                # Calculate normalized density, temperature, pressure and entropy profiles
                if CALC_THERMO_PROPS or CLASSIFY_HALOS:

                    ## Entropy and pressure profiles calculated directly from the temperature and ne profiles
                    for weight in weighted_prop_weights + ['spec']:
                        if weight != 'T':
                            try:
                                profile['kT_' + weight] = kB * profile['T_' + weight]

                                profile['K_' + weight + '_bin_T'] = profile['kT_' + weight] * profile['ne_bin']**(-2,3)
                                profile['P_' + weight + '_bin_T'] = profile['kT_' + weight] * profile['ne_bin']

                                if weight not in ['ne', 'spec']:
                                    profile['K_' + weight + '_bin_T_ne'] = profile['kT_' + weight] * profile['ne_' + weight]**(-2,3) 
                                    profile['P_' + weight + '_bin_T_ne'] = profile['kT_' + weight] * profile['ne_' + weight]
                                    
                            except:
                                continue



                    ## Cooling time profiles calculated directly from the ne, T, Lx, volume, fe profiles
                    for weight in weighted_prop_weights + ['spec']:
                        if weight != 'T':
                            try:
                                for xband in xray_bands:
                                    profile['tcool_' + xband + '_' + weight + '_bin'] = 6 * profile['ne_bin'] * profile['kT_' + weight] / (2 * profile['fe'] * profile['Lx_'+ xband] / profile['volume'])
                                    
                            except:
                                continue


                    
                    if CALC_THERMO_PROPS:
                        ## Scaling properties by their relevant scaling factors
                        for prop, scaling_factors in thermo_scaling_dict.items():
                            scale_prop = scaling_factors['scale_prop']
                            scale_types = scaling_factors['scale_types']

                            for scale_type in scale_types:
                                scaled_prop_name = prop + '_' + scale_prop + scale_type

                                ## CHANGE THIS SO CAN SCALE BY VALUES THAT HAVE dIFFERENT UNITS THAN PROP
                                try:
            #                         profile[scaled_prop_name] = profile[prop] / profiles['temp_halo_info'][scale_prop][scale_type].in_units(profile[prop].units)[halo_idx]
                                    profile[scaled_prop_name] = profile[prop] / profile_to_compute['temp_halo_info'][scale_prop][scale_type].in_units(profile[prop].units)[halo_idx]

                                    profile[scaled_prop_name].units = ''

                                    if scale_prop in ['T', 'K', 'P']:
                                        profile[scaled_prop_name+'_v2'] = 2*profile[scaled_prop_name]

                                except:
                                    continue






                ## Calculate extra properties related to metallicity
                if CALC_METAL_PROPS and compute_metallicity:

                    for metallicity_type in metallicity_types:
                        for weight in weighted_prop_weights:
                            
                            try:
#                                 print()
#                                 print()
#                                 print(weight.upper() + ' start')
#                                 print()
#                                 print()

                                ## Correct metallicity units from NoUnit to 1
                                profile[metallicity_type + '_' + weight].units = ''


                                ## Calculate metal ratios
                                metal_profile = profile[metallicity_type + '_' + weight]
                                Z_shape = np.shape(metal_profile)
                                Z_num = np.repeat(np.expand_dims(metal_profile, axis=2), repeats=Z_shape[1], axis=2)
                                Z_denom = np.swapaxes(Z_num, axis1=1, axis2=2)
                                Z_ratio = Z_num/Z_denom
                                Z_ratio = pnb.array.SimArray(Z_ratio, units='')
                                profile[metallicity_type + '_' + weight + '_ratio'] = Z_ratio

                                ## Check metal ratio calculations
            #                     Z_profile = profile[metallicity_type + '_' + weight]
            #                     Z_shape = np.shape(metal_profile)
            #                     n_rbins = Z_shape[0]
            #                     n_metals = Z_shape[1]
            #                     Z_ratio = pnb.array.SimArray(np.zeros((n_rbins, n_metals, n_metals)), units='')
            #                     for element_num in range(n_metals):
            #                         for element_denom in range(n_metals):
            #                             Z_ratio[:, element_num, element_denom] = Z_profile[:,element_num]/Z_profile[:,element_denom]
            #                     profile[metallicity_type + '_' + weight + '_ratio_v2'] = Z_ratio


                                ## Metallicity and metallicity ratios relative to solar metallicity
                                for normalization, norm_info in metallicity_scaling_dict.items():
                                    profile[metallicity_type+'_'+weight+'_solar_'+normalization] = copy.deepcopy(profile[metallicity_type+'_' +weight])
                                    profile[metallicity_type+'_'+weight+'_ratio_solar_'+normalization] = copy.deepcopy(profile[metallicity_type+'_' +weight+'_ratio'])
            #                         profile[metallicity_type+'_'+weight+'_ratio_v2_solar_'+normalization] = copy.deepcopy(profile[metallicity_type+'_' +weight+'_ratio_v2'])

                                    if metallicity_type.lower() == 'z' and weight != 'median':
                                        profile['Ztot_' + weight + '_solar_' + normalization] = profile['Ztot_' + weight] / norm_info['total']

                                    for element, element_idx in metals_idx[CODE].items():
                                        profile[metallicity_type+'_'+weight+'_solar_'+normalization][:,element_idx] /= norm_info['info'][element][norm_info[metallicity_type+'_name']]

                                        for element_denom, element_denom_idx in metals_idx[CODE].items():
                                            profile[metallicity_type+'_'+weight+'_ratio_solar_'+normalization][:,element_idx,element_denom_idx] /= norm_info['info'][element][norm_info[metallicity_type+'_name']]/norm_info['info'][element_denom][norm_info[metallicity_type+'_name']]
            #                                 profile[metallicity_type+'_'+weight+'_ratio_v2_solar_'+normalization][:,element_idx,element_denom_idx] /= norm_info['info'][element][norm_info[metallicity_type+'_name']]/norm_info['info'][element_denom][norm_info[metallicity_type+'_name']]
            
#                                 print()
#                                 print()
#                                 print(weight.upper() + ' finish')
#                                 print()
#                                 print()
            
            
                            except:
                                continue




                    ## Calculate metal mass and metal mass volume density profiles

        #             igrm_mass_profile = profile['mass']
        #             metal_mass_profile = copy.deepcopy(profile['Z_mass'])
        # #             metal_mass_enc_profile = copy.deepcopy(profile['Z_mass'])
        # #             metal_mass_density_profile = copy.deepcopy(profile['Z_mass'])
        # #             metal_mass_density_enc_profile = copy.deepcopy(profile['Z_mass'])
        # #             metal_mass_density_rhoIGrM_profile = copy.deepcopy(profile['Z_mass'])
        # #             metal_mass_density_enc_rhoIGrMenc_profile = copy.deepcopy(profile['Z_mass'])
        #             n_rbins, n_elements = np.shape(metal_mass_profile)
        #             for rbin in range(n_rbins):
        #                 metal_mass_profile[rbin] *= igrm_mass_profile[rbin]
        # #                 metal_mass_density_profile[rbin] = metal_mass_profile[rbin]/profile['volume']
        # #                 metal_mass_density_rhoIGrM_profile[rbin] = metal_mass_density_profile[rbin]/profile['mass_density'][rbin]
        #             profile['metal_mass_v2'] = pnb.array.SimArray(metal_mass_profile, units=igrm_mass_profile.units)
        # #             profile['metal_mass_density'] = pnb.array.SimArray(metal_mass_density_profile, 
        # #                                                                units=profile['metal_mass'].units/profile['volume'].units)
        # #             profile['metal_mass_enc'] = np.nancumsum(profile['metal_mass'], axis=0)
        # #             for rbin in range(n_rbins):
        # #                 metal_mass_density_enc_profile[rbin] = profile['metal_mass_enc'][rbin]/profile['volume_enc'][rbin]


                    expanded_mass = np.expand_dims(profile['mass'], 
                                                   tuple(i for i in range(profile['Z_mass'].ndim) if i != RADIAL_AXIS))
                    expanded_mass_enc = np.expand_dims(profile['mass_enc'], 
                                                       tuple(i for i in range(profile['Z_mass'].ndim) if i != RADIAL_AXIS))
                    expanded_mass_density = np.expand_dims(profile['mass_density'], 
                                                           tuple(i for i in range(profile['Z_mass'].ndim) if i != RADIAL_AXIS))
                    expanded_mass_density_enc = np.expand_dims(profile['mass_density_enc'], 
                                                               tuple(i for i in range(profile['Z_mass'].ndim) if i != RADIAL_AXIS))
                    expanded_volume = np.expand_dims(profile['volume'], 
                                                     tuple(i for i in range(profile['Z_mass'].ndim) if i != RADIAL_AXIS))
                    expanded_volume_enc = np.expand_dims(profile['volume_enc'], 
                                                         tuple(i for i in range(profile['Z_mass'].ndim) if i != RADIAL_AXIS))

                    profile['metal_mass'] = profile['Z_mass'] * expanded_mass
                    profile['metal_mass_enc'] = np.nancumsum(profile['metal_mass'], axis=RADIAL_AXIS)
                    profile['metal_mass_density'] = profile['metal_mass']/expanded_volume
                    profile['metal_mass_density_enc'] = profile['metal_mass_enc']/expanded_volume_enc

                    profile['metal_mass_Migrm'] = profile['metal_mass']/expanded_mass ## should be same as Z_mass
                    profile['metal_mass_enc_Migrmenc'] = profile['metal_mass_enc']/expanded_mass_enc

                    profile['metal_mass_density_rhoigrm'] = profile['metal_mass_density']/expanded_mass_density ## also should be same as Z_mass
                    profile['metal_mass_density_enc_rhoigrmenc'] = profile['metal_mass_density_enc']/expanded_mass_density_enc

                    profile['metal_mass_density_fb_rhocrit'] = profile['metal_mass_density']/(profile_to_compute['cosmo_props']['rho_gas_scaling'])
                    profile['metal_mass_density_enc_fb_rhocrit'] = profile['metal_mass_density_enc']/(profile_to_compute['cosmo_props']['rho_gas_scaling'])

                    profile['metal_mass_density_Zmass_fb_rhocrit'] = profile['metal_mass_density']/(profile['Z_mass'] * profile_to_compute['cosmo_props']['rho_gas_scaling'])
                    profile['metal_mass_density_enc_Zmassenc_fb_rhocrit'] = profile['metal_mass_density_enc']/(profile['metal_mass_enc_Migrmenc'] * profile_to_compute['cosmo_props']['rho_gas_scaling'])

                    profile['mass_density_X_Y_fb_rhocrit'] = profile['mass_density'] / (profile['Z_mass'][:,metals_idx[CODE]['H']] * profile['Z_mass'][:,metals_idx[CODE]['He']] * profile_to_compute['cosmo_props']['rho_gas_scaling'])
                    profile['mass_density_enc_Xenc_Yenc_fb_rhocrit'] = profile['mass_density_enc'] / (profile['metal_mass_enc_Migrmenc'][:,metals_idx[CODE]['H']] * profile['metal_mass_enc_Migrmenc'][:,metals_idx[CODE]['He']] * profile_to_compute['cosmo_props']['rho_gas_scaling'])


                    profile['total_metal_mass'] = np.nansum(profile['metal_mass'][:,2:], axis=1) ## exclude H and He
                    profile['total_metal_mass_enc'] = np.nancumsum(profile['total_metal_mass'], axis=0)
                    profile['total_metal_mass_density'] = profile['total_metal_mass']/expanded_volume
                    profile['total_metal_mass_density_enc'] = profile['total_metal_mass_enc']/expanded_volume_enc

                    profile['total_metal_mass_density_rhoigrm'] = profile['total_metal_mass_density']/expanded_mass_density ## also should be same as Z_mass_tot
                    profile['total_metal_mass_density_enc_rhoigrmenc'] = profile['total_metal_mass_density_enc']/expanded_mass_density_enc

                    profile['total_metal_mass_density_fb_rhocrit'] = profile['total_metal_mass_density']/(profile_to_compute['cosmo_props']['rho_gas_scaling'])
                    profile['total_metal_mass_density_enc_fb_rhocrit'] = profile['total_metal_mass_density_enc']/(profile_to_compute['cosmo_props']['rho_gas_scaling'])

                    ## Add extra axis so can be plotted with other metals as if it is an element, by giving its element as H (ie. index 0)
                    profile['total_metal_mass_for_metals'] = profile['total_metal_mass'][:,np.newaxis]
                    profile['total_metal_mass_enc_for_metals'] = profile['total_metal_mass_enc'][:,np.newaxis]
                    profile['total_metal_mass_density_for_metals'] = profile['total_metal_mass_density'][:,np.newaxis]
                    profile['total_metal_mass_density_enc_for_metals'] = profile['total_metal_mass_density_enc'][:,np.newaxis]
                    profile['total_metal_mass_density_rhoigrm_for_metals'] = profile['total_metal_mass_density_rhoigrm'][:,np.newaxis]
                    profile['total_metal_mass_density_enc_rhoigrmenc_for_metals'] = profile['total_metal_mass_density_enc_rhoigrmenc'][:,np.newaxis]
                    profile['total_metal_mass_density_fb_rhocrit_for_metals'] = profile['total_metal_mass_density_fb_rhocrit'][:,np.newaxis]
                    profile['total_metal_mass_density_enc_fb_rhocrit_for_metals'] = profile['total_metal_mass_density_enc_fb_rhocrit'][:,np.newaxis]





                # Add in log versions of all properties
                if CALC_LOG_PROPS:
                    log_prop_names = []
                    log_prop_vals = []
                    for prop in profile.keys():
                #             if prop not in x_axis_props:
                        log_prop_names.append('log_'+prop)
                        log_prop_vals.append(np.log10(profile[prop]))

                    for log_prop_name, log_prop_val in zip(log_prop_names, log_prop_vals):
                        profile[log_prop_name] = log_prop_val


    #             print()




    #             ## Any classification, including cool core (CC) vs non-cool core (NCC)
    #             if CLASSIFY_HALOS and profile_type_to_compute.lower() == 'all_profiles':

    #         #         profile = profiles['all_profiles'][halo]['y']

    #                 for classification_type, classification_info in classifications_dict.items():

    #                     profile_to_classify = classification_info['profile_to_use'][classification_info['profile_type']][halo]['y']

    #                     for ii in range(classification_info['num']):
    #                         kwargs = {kwarg:classification_info['kwargs'][kwarg][ii] for kwarg in classification_info['kwargs'].keys()}
    #         #                 print(kwargs)

    #                         try:
    #         #                     if 'tcool' in classification_type:
    #         #     #                     kwargs['profile'] = profiles['all_core_profiles'][halo]['y']
    #         #                         this_classification = classification_info['classification_function'](profiles['all_core_profiles'][halo]['y'], 
    #         #                                                                                              kwargs)
    #         #                     elif 'temp' in classification_type:
    #         #                         this_classification = classification_info['classification_function'](profile, 
    #         #                                                                                              profiles['all_core_profiles'][halo]['y'], 
    #         #                                                                                              kwargs)

    #                             this_classification = classification_info['classification_function'](profile_to_classify, kwargs)

    #                         except:
    #                             print('No core profile for halo', halo, 'for CC/NCC', classifiation_type, 'classification')
    #                             continue

    #                         profile_to_compute[this_classification + '_' + 
    #                                  classification_type + '_' + 
    #                                  kwargs['prop_type'] + '_profiles'][halo] = copy.deepcopy(profile_to_compute['all_profiles'][halo])
    #                         profile_to_compute[this_classification + '_' + 
    #                                  classification_type + '_' + 
    #                                  kwargs['prop_type'] + '_id'] = np.append(profile_to_compute[this_classification + '_' + 
    #                                                                                    classification_type + '_' + 
    #                                                                                    kwargs['prop_type'] + '_id'], halo)


#                 ## delete certain properties if want to save space

#                 for prop_to_remove in x_props_to_remove:
#                     try:
#                         del profile_to_compute[profile_type_to_compute][halo]['x'][prop_to_remove]
#                     except:
#                         print('Either', profile_type_to_compute, 'does not exist in', profile_to_compute, 'or', prop_to_remove, 'does not exist in', profile_type_to_compute)

#                 for prop_to_remove in y_props_to_remove:
#                     try:
#                         del profile_to_compute[profile_type_to_compute][halo]['y'][prop_to_remove]
#                     except:
#                         print('Either', profile_type_to_compute, 'does not exist in', profile_to_compute, 'or', prop_to_remove, 'does not exist in', profile_type_to_compute)


        #         for prop_to_remove in x_props_to_remove:
        #             if halo in profiles['all_core_id']:
        #                 try:
        #                     del profiles['all_core_profiles'][halo]['x'][prop_to_remove]
        #                 except:
        #                     print(prop_to_remove, 'does not exist in all core profiles')

        #             for classification in profiles['classifications']:
        #                 if halo in profiles[classification+'_id']:
        #                     try:
        #                         del profiles[classification+'_profiles'][halo]['x'][prop_to_remove]
        #                     except:
        #                         print(prop_to_remove, 'does not exist in', classification, 'profiles')


    #             for prop_to_remove in y_props_to_remove:
    #                 if halo in profiles['all_core_id']:
    #                     try:
    #                         del profiles['all_core_profiles'][halo]['y'][prop_to_remove]
    #                         if CALC_LOG_PROPS and 'all_core_profiles' in profile_types_to_compute:
    #                             del profiles['all_core_profiles'][halo]['y']['log_'+prop_to_remove]
    #                     except:
    #                         print(prop_to_remove, 'does not exist in all core profiles')

    #                 for classification in profiles['classifications']:
    #                     if halo in profiles[classification+'_id']:
    #                         try:
    #                             del profiles[classification+'_profiles'][halo]['y'][prop_to_remove]
    #                             if CALC_LOG_PROPS and 'all_profiles' in profile_types_to_compute:
    #                                 del profiles[classification+'_profiles'][halo]['y']['log_'+prop_to_remove]
    #                         except:
    #                             print(prop_to_remove, 'does not exist in', classification, 'profiles')


            print()
            
            
print()
print('done')
print()
print()
            
            
        
        
        
        
        
        
        
def determine_resolved(counts, resolved_counts=1, iterate_over_bins=False):
    if iterate_over_bins:
        resolved_filter = np.full(counts.shape, True)
        for ii in range(counts.shape[0]):
            for jj in range(ii, counts.shape[0]):
                curr_filter = np.expand_dims(counts[jj], axis=0) >= resolved_counts
                resolved_filter[ii] = resolved_filter[ii] & curr_filter
    else:
        resolved_filter = counts >= resolved_counts
    return resolved_filter


        
## For resolution cuts in individual radial bins of individual halos
if INCLUDE_RESOLUTION_CUTS:

#     def apply_func_to_dict_values(d, key_name, func, func_kwargs={}, extra_keys=[]):
#         for k, v in d.items():
#             if isinstance(v, dict):
#                 apply_func_to_dict_values(v, key_name, func, 
#                                           func_kwargs=func_kwargs, extra_keys=extra_keys)
#             else:
#                 if k == key_name:
#                     extra_kwargs = {extra_key:d[extra_key] for extra_key in extra_keys}
#                     func(v, **func_kwargs, **extra_kwargs)

#     def remove_values_from_averaging(dict_of_values, kwargs):
#         ## Remove negative values from averaging by setting them to NaN
#         ## Check - for each radial bin - that at least 16% of halos have at least N particles --> if not, set all values in that bin to NaN
#         ## Also make sure there at at least M halos in each radial bin
#         ## When a radial bin meets the above two requirements, do I include the unresolved halos in the averaging still?
#         ## --> yes, since they will be weighted less in the averaging due to have fewer particles (so increased poisson error)    

#         print('Removing values from averaging')
#         for function, info in kwargs['functions'].items():
#             print()
#             print(function)
#             print()
#             apply_func_to_dict_values(dict_of_values, kwargs['key_name'], function, 
#                                       func_kwargs=info['func_kwargs'], extra_keys=info['extra_keys'])


    def exclude_negatives(values, replace_with=np.nan):
        values[values<0] = replace_with

    def exclude_unresolved(values, n=None, resolved_counts=1, iterate_over_bins=False, replace_with=np.nan):
        resolved_filter = determine_resolved(n, resolved_counts=resolved_counts, iterate_over_bins=iterate_over_bins)
        values[~resolved_filter] = replace_with
    #     if iterate_over_bins:
    #         for ii in range(values.shape[0]):
    #             for jj in range(ii, values.shape[0]):
    # #                 if np.any(np.expand_dims(counts[jj], axis=0) < resolved_counts, axis=0):
    #                 unresolved_filter = np.expand_dims(counts[jj], axis=0) < resolved_counts
    #                 values[unresolved_filter] = replace_with
    #     else:
    #         values[counts<resolved_counts] = replace_with

    def exclude_less_than_length(values, physical_rbins=None, length=0, replace_with=np.nan):
        for xbin in range(len(values)):
            ## This way only excludes the innermost bin whose outer edge is less than the grav softening length
    #         if physical_xbins[xbin+1] <= length:
    #             values[xbin] = replace_with
            ## This way excludes all bins that are thinner than the grav softening length
            if physical_rbins[xbin+1] - physical_rbins[xbin] <= length:
                values[xbin] = replace_with


# if INCLUDE_RESOLUTION_CUTS:
#     ind_res_kwargs = {
# #         'key_name':'values',
#         'functions':{
#             exclude_negatives:{
#                 'func_kwargs':{'replace_with':np.nan},
#                 'extra_keys':{},
#             },
#             exclude_unresolved:{
#                 'func_kwargs':{
#                     'resolved_counts':MIN_RESOLVED_COUNTS_PER_HALO, 
#                     'iterate_over_bins':PROFILES_ARE_2D_BUT_RES_FILE_IS_3D, 
#                     'replace_with':np.nan,
#                 },
#                 'extra_keys':{
#                     'res_profile':'n',
#                 },
#             },
#             exclude_less_than_length:{
#                 'func_kwargs':{
#                     'length':GRAV_SOFTENING_LENGTH.in_units(XSCALE_UNITS_DICT['Physical']),
#                     'replace_with':np.nan
#                 },
#                 'extra_keys':{
#                     'x_profile':'physical_xbins',
#                 },
#             },
#         }
#     }


    
    ## Apply resolution cuts to selected profiles
    print()
    print('Applying resolution cuts to selected halos')
    print()

    profiles_to_compute_names = ['three_d_profiles', 'two_d_profiles', 'profiles']
    profiles_to_compute = [three_d_profiles, two_d_profiles, profiles]
    profiles_ndims = [3, 2, NDIM]

    profile_types_to_compute = ['all_core_profiles', 'all_profiles']
    profile_ids_to_compute = ['all_core_id', 'all_id']

    for profile_name, profile_to_compute, profile_ndim in zip(profiles_to_compute_names, profiles_to_compute, profiles_ndims):
        print()
        print()
        print(profile_name)
        print()

        for profile_type_to_compute, profile_id_to_compute in zip(profile_types_to_compute, profile_ids_to_compute):
            print()
            print()
            print(profile_type_to_compute)
            print()

            for halo in profile_to_compute[profile_id_to_compute]:
#                 print('halo', halo)
                halo_idx = list(profile_to_compute['temp_halo_info']['halo_id']).index(halo)
#                 print('halo_idx:', halo_idx)
                
                curr_profile = {}
                
                try:
                    x_profile = profile_to_compute[profile_type_to_compute][halo]['x']
                    y_profile = profile_to_compute[profile_type_to_compute][halo]['y']
                except:
                    print(profile_type_to_compute, 'does not exist for halo', halo)
                    continue
                    
                curr_profile['x_profile'] = x_profile
                curr_profile['y_profile'] = y_profile
                
#                 has_res_profile = True
#                 try:
#                     curr_profile['res_profile'] = res_profiles[profile_type_to_compute][halo]['y']#['n']
#                 except:
#                     has_res_profile = False
# #                 n_particles = y_profile['n']

                has_three_d_profile = True
                try:
                    curr_profile['three_d_profile'] = three_d_profiles[profile_type_to_compute][halo]['y']#['n']
                except:
                    has_three_d_profile = False


                ind_res_kwargs = {
                    'functions':{
                        exclude_negatives:{
                            'func_kwargs':{'replace_with':np.nan},
                            'extra_keys':{},
                        },
                        exclude_unresolved:{
                            'func_kwargs':{
                                'resolved_counts':MIN_RESOLVED_COUNTS_PER_HALO if profile_type_to_compute!='all_core_profiles' else MIN_RESOLVED_COUNTS_PER_HALO_CORE, 
                                'iterate_over_bins': False if profile_ndim==3 else USE_3D_FOR_2D_RES, 
                                'replace_with':np.nan,
                            },
                            'extra_keys':{
                                'y_profile' if (profile_ndim==3 or not USE_3D_FOR_2D_RES) else 'three_d_profile':'n',
#                                 'res_profile':'n',
                            },
                        },
#                         exclude_less_than_length:{
#                             'func_kwargs':{
#                                 'length':GRAV_SOFTENING_LENGTH.in_units(XSCALE_UNITS_DICT['Physical']),
#                                 'replace_with':np.nan
#                             },
#                             'extra_keys':{
#                                 'x_profile':'physical_rbins',
#                             },
#                         },
                    }
                }
                
                if profile_type_to_compute!='all_core_profiles':
                    ind_res_kwargs[exclude_less_than_length] = {
                        'func_kwargs':{
                            'length':GRAV_SOFTENING_LENGTH.in_units(XSCALE_UNITS_DICT['Physical']),
                            'replace_with':np.nan
                        },
                        'extra_keys':{
                            'x_profile':'physical_rbins',
                        },
                    }
    
#                 print()
#                 print(ind_res_kwargs)
#                 print()


                for prop, vals in curr_profile['y_profile'].items():
#                     print(prop)
                    set_nonfinite_to_nan(vals)
                    if prop == 'n': continue
                    for func, args in ind_res_kwargs['functions'].items():
                        if 'three_d_profile' in args['extra_keys'].keys() and not has_three_d_profile:
                            continue
                            
                        extra_kwargs = {extra_key:curr_profile[x_or_y_profile][extra_key] for x_or_y_profile, extra_key in args['extra_keys'].items()}
#                         print()
#                         print(profile_to_compute[profile_type_to_compute][halo]['y'][prop])
                        func(vals, **args['func_kwargs'], **extra_kwargs)
#                         print(profile_to_compute[profile_type_to_compute][halo]['y'][prop])
#                         print()
                        
#                     print()
                
                
                
        
        
        
            
        
        
        
        
        
        
## different ways of classifying groups (currently just CC/NCC)

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## do I need to apply resolution cuts to individual radial bins in individual groups before performing CC/NCC classification? -> yes?
## This includes core profiles --> if too few points, group does not get classified by that classification method
## Also need to deal with case where there are 0 groups in a classification (produces error farther down) --> dealt with
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def CC_classification_temp_frac(profile, kwargs):
    ## My adaptation of O'Sullivan+2017 method for groups
    
    temp_type = kwargs['prop_type']
    CC_core_temp_frac_of_peak_temp = kwargs['core_temp_frac_of_peak_temp']
    
    peak_temp = np.nanmax(profile['y'][temp_type])
    core_temp = profile['y'][temp_type][np.isfinite(profile['y'][temp_type])][0]
    
    frac_temp = core_temp/peak_temp
    
#     return 'CC' if frac_temp <= CC_core_temp_frac_of_peak_temp else 'NCC'
    if frac_temp <= CC_core_temp_frac_of_peak_temp:
        return 'CC'
    elif frac_temp > CC_core_temp_frac_of_peak_temp:
        return 'NCC'
    else:
        return None

# def CC_classification_temp_frac_v2(profile_peak, profile_core, kwargs):
#     ## My adaptation of O'Sullivan+2017 method for groups
    
#     temp_type = kwargs['prop_type']
#     CC_core_temp_frac_of_peak_temp = kwargs['core_temp_frac_of_peak_temp']
    
#     peak_temp = np.nanmax(profile_peak[temp_type])
#     core_temp = profile_core[temp_type][np.isfinite(profile_core[temp_type])][0]
    
#     frac_temp = core_temp/peak_temp
    
#     return 'CC' if frac_temp <= CC_core_temp_frac_of_peak_temp else 'NCC'


def CC_classification_temp_sigma(profile, kwargs):
    ## O'Sullivan+2017 method for groups
    
    temp_type = kwargs['prop_type']
    CC_sigma = kwargs['sigma']
    
    peak_temp = np.nanmax(profile['y'][temp_type])
    core_temp = profile['y'][temp_type][np.isfinite(profile['y'][temp_type])][0]
    
#     poisson_err_peak_temp = calc_poisson_err(peak_temp, profile['n'][np.nanargmax(profile[temp_type])], log=False)
#     poisson_err_core_temp = calc_poisson_err(core_temp, profile['n'][np.isfinite(profile[temp_type])][0], log=False)
    
    poisson_err_peak_temp = peak_temp/calc_poisson_err(profile['y']['n'][np.nanargmax(profile['y'][temp_type])], log=False)
    poisson_err_core_temp = core_temp/calc_poisson_err(profile['y']['n'][np.isfinite(profile['y'][temp_type])][0], log=False)
    
    sigma_temp = np.abs(peak_temp - core_temp)/np.sqrt(poisson_err_peak_temp**2 + poisson_err_core_temp**2)
    
#     return 'CC' if sigma_temp >= CC_sigma else 'NCC'
    if sigma_temp >= CC_sigma:
        return 'CC'
    elif sigma_temp < CC_sigma:
        return 'NCC'
    else:
        return None

# def CC_classification_temp_sigma_v2(profile_peak, profile_core, kwargs):
#     ## O'Sullivan+2017 method for groups
    
#     temp_type = kwargs['prop_type']
#     CC_sigma = kwargs['sigma']
    
#     peak_temp = np.nanmax(profile_peak[temp_type])
#     core_temp = profile_core[temp_type][np.isfinite(profile_core[temp_type])][0]
    
#     poisson_err_peak_temp = calc_poisson_err(peak_temp, profile_peak['n'][np.nanargmax(profile_peak[temp_type])], log=False)
#     poisson_err_core_temp = calc_poisson_err(core_temp, profile_core['n'][np.isfinite(profile_core[temp_type])][0], log=False)
    
#     sigma_temp = np.abs(peak_temp - core_temp)/np.sqrt(poisson_err_peak_temp**2 + poisson_err_core_temp**2)
    
#     return 'CC' if sigma_temp >= CC_sigma else 'NCC'


def CC_classification_tcool(profile, kwargs):
    ## Braspenning+2024 method for clusters
    ## 1 Gyr, 5 Gyr from Braspenning+2024
    ## 1 Gyr, 7.7 Gyr from Lehle+2024
    
    tcool_type = kwargs['prop_type']
    radial_bin_for_tcool = kwargs['radial_bin']
    CC_tcool = kwargs['CC_tcool']
    NCC_tcool = kwargs['NCC_tcool']
    
    tcool = profile['y'][tcool_type].in_units(CC_tcool.units)[radial_bin_for_tcool]
    
    if tcool <= CC_tcool:
        return 'CC'
    elif tcool >= NCC_tcool:
        return 'NCC'
    elif tcool > CC_tcool and tcool < NCC_tcool:
        return 'Weak_CC'
    else:
        return None
    

def CC_classification_entropy(profile, kwargs):
    ## Lehle+2024
    
    entropy_type = kwargs['prop_type']
    radial_bin_for_entropy = kwargs['radial_bin']
    CC_entropy = kwargs['CC_entropy']
    NCC_entropy = kwargs['NCC_entropy']
    
    entropy = profile['y'][entropy_type].in_units(CC_entropy.units)[radial_bin_for_entropy]
    
    if entropy <= CC_entropy:
        return 'CC'
    elif entropy > NCC_entropy:
        return 'NCC'
    elif entropy > CC_entropy and entropy <= NCC_entropy:
        return 'Weak_CC'
    else:
        return None
    
    
def CC_classification_ne(profile, kwargs):
    ## Lehle+2024
    
    ne_type = kwargs['prop_type']
    radial_bin_for_ne = kwargs['radial_bin']
    CC_ne = kwargs['CC_ne']
    NCC_ne = kwargs['NCC_ne']
    
    ne = profile['y'][ne_type].in_units(CC_ne.units)[radial_bin_for_ne]
    
    if ne > CC_ne:
        return 'CC'
    elif ne <= NCC_ne:
        return 'NCC'
    elif ne <= CC_ne and ne > NCC_ne:
        return 'Weak_CC'
    else:
        return None
    
    
def CC_classification_Cphys(profile, kwargs):
    ## Lehle+2024
    
    Lx_type = kwargs['prop_type']
#     inner_radial_bin = kwargs['inner_radial_bin']
#     outer_radial_bin = kwargs['outer_radial_bin']
    inner_radius = kwargs['inner_radius']
    outer_radius = kwargs['outer_radius']
    CC_Cphys = kwargs['CC_Cphys']
    NCC_Cphys = kwargs['NCC_Cphys']
    
    rbins = profile['x']['physical_rbins'].in_units(inner_radius.units)
    
    ## Get radial bins -> always round up to next outer bin (better than using fewer bins)
    ## unless the bin radius is larger than the max radius of the group
    if inner_radius >= np.nanmax(rbins):
        inner_radial_bin = len(rbins)
    else:
        closest_inner_radius = rbins[(rbins - inner_radius) >= 0][0]
        inner_radial_bin = list(rbins).index(closest_inner_radius)
    
    if outer_radius >= np.nanmax(rbins):
        outer_radial_bin = len(rbins)
    else:
        closest_outer_radius = rbins[(rbins - outer_radius) >= 0][0]
        outer_radial_bin = list(rbins).index(closest_outer_radius)
    
#     Lx_inner = np.nancumsum(profile[Lx_type][:inner_radial_bin])
#     Lx_outer = np.nancumsum(profile[Lx_type][:outer_radial_bin])
    
    Lx_inner = profile['y'][Lx_type][inner_radial_bin]
    Lx_outer = profile['y'][Lx_type][outer_radial_bin]
    
    Cphys = Lx_inner/Lx_outer
    
    if Cphys > CC_Cphys:
        return 'CC'
    elif Cphys <= NCC_Cphys:
        return 'NCC'
    elif Cphys <= CC_Cphys and Cphys > NCC_Cphys:
        return 'Weak_CC'
    else:
        return None
    
    
def CC_classification_Cscaled(profile, kwargs):
    ## Lehle+2024
    
    Lx_type = kwargs['prop_type']
#     inner_radial_bin = kwargs['inner_radial_bin']
#     outer_radial_bin = kwargs['outer_radial_bin']
    inner_radius = kwargs['inner_radius']
    outer_radius = kwargs['outer_radius']
    CC_Cscaled = kwargs['CC_Cscaled']
    NCC_Cscaled = kwargs['NCC_Cscaled']
    
    rbins = profile['x']['rbins_r500']
    
    ## Get radial bins -> always round up to next outer bin (better than using fewer bins)
    ## unless the bin radius is larger than the max radius of the group
    if inner_radius >= np.nanmax(rbins):
        inner_radial_bin = len(rbins)
    else:
        closest_inner_radius = rbins[(rbins - inner_radius) >= 0][0]
        inner_radial_bin = list(rbins).index(closest_inner_radius)
    
    if outer_radius >= np.nanmax(rbins):
        outer_radial_bin = len(rbins)
    else:
        closest_outer_radius = rbins[(rbins - outer_radius) >= 0][0]
        outer_radial_bin = list(rbins).index(closest_outer_radius)
    
#     Lx_inner = np.nancumsum(profile[Lx_type][:inner_radial_bin])
#     Lx_outer = np.nancumsum(profile[Lx_type][:outer_radial_bin])
    
    Lx_inner = profile['y'][Lx_type][inner_radial_bin]
    Lx_outer = profile['y'][Lx_type][outer_radial_bin]
    
    Cscaled = Lx_inner/Lx_outer
    
    if Cscaled > CC_Cscaled:
        return 'CC'
    elif Cscaled <= NCC_Cscaled:
        return 'NCC'
    elif Cscaled <= CC_Cscaled and Cscaled > NCC_Cscaled:
        return 'Weak_CC'
    else:
        return None


classifications_dict = {
    ## Temperature-based CC/NCC classifications as in CLoGS O'Sullivan+2017
    'temp_frac':{
        'classifications':['CC', 'NCC'],
        'classification_function':CC_classification_temp_frac,
        'profile_to_use':three_d_profiles,
        'profile_type':'all_profiles',
        'num':3,
        'kwargs':{
            'profile':['all_profiles', 'all_profiles', 'all_profiles',],
            'prop_type':['T_spec', 'T_Lx_0.5-10.0keV_total', 'T_mass'],
            'core_temp_frac_of_peak_temp':[0.75, 0.75, 0.75],
        },
    },
    'temp_sigma':{
        'classifications':['CC', 'NCC'],
        'classification_function':CC_classification_temp_sigma,
        'profile_to_use':three_d_profiles,
        'profile_type':'all_profiles',
        'num':3,
        'kwargs':{
            'profile':['all_profiles', 'all_profiles', 'all_profiles',],
            'prop_type':['T_spec', 'T_Lx_0.5-10.0keV_total', 'T_mass'],
            'sigma':[3, 3, 3],
        },
    },
    ## Cooling time-based CC/NCC classification as in Braspenning+2024
    'tcool':{
        'classifications':['CC', 'Weak_CC', 'NCC'],
        'classification_function':CC_classification_tcool,
        'profile_to_use':three_d_profiles,
        'profile_type':'all_core_profiles',    ## <0.048R500
        'num':2,   ## Can't have multiple of same classification methov with same prop type currently!
        'kwargs':{
            'prop_type':['tcool_0.5-10.0keV_total_Lx_0.5-10.0keV_total', 'tcool_0.5-10.0keV_total_Lx_0.5-10.0keV_total_bin',
                         'tcool_0.5-10.0keV_total_Lx_0.5-10.0keV_total', 'tcool_0.5-10.0keV_total_Lx_0.5-10.0keV_total_bin'],
            'radial_bin':[0, 0, 0, 0],
            'CC_tcool':[pnb.array.SimArray(1, units='Gyr'), pnb.array.SimArray(1, units='Gyr'), 
                        pnb.array.SimArray(1, units='Gyr'), pnb.array.SimArray(1, units='Gyr')],
            'NCC_tcool':[pnb.array.SimArray(7.7, units='Gyr'), pnb.array.SimArray(7.7, units='Gyr'),  ## Lehle+2024
                         pnb.array.SimArray(5, units='Gyr'), pnb.array.SimArray(5, units='Gyr')],  ## Braspenning+2024
        },
    },
    ## Core entropy-based CC/NCC classification as in Lehle+2024
    'entropy':{
        'classifications':['CC', 'Weak_CC', 'NCC'],
        'classification_function':CC_classification_entropy,
        'profile_to_use':three_d_profiles,
        'profile_type':'all_core_profiles',
        'num':4,
        'kwargs':{
            'prop_type':['K_mass', 'K_Lx_0.5-10.0keV_total', 'K_mass_bin_T_ne', 'K_Lx_0.5-10.0keV_total_bin_T_ne'],
            'radial_bin':[0, 0, 0, 0],  ## Lehle+2024 have all their cores as <0.012R500c
            'CC_entropy':[pnb.array.SimArray(22, units='keV cm**2'), pnb.array.SimArray(22, units='keV cm**2'), 
                          pnb.array.SimArray(22, units='keV cm**2'), pnb.array.SimArray(22, units='keV cm**2')],
            'NCC_entropy':[pnb.array.SimArray(150, units='keV cm**2'), pnb.array.SimArray(150, units='keV cm**2'), 
                           pnb.array.SimArray(150, units='keV cm**2'), pnb.array.SimArray(150, units='keV cm**2')],  ## Lehle+2024
            ## Vogelsberger+2018 use core entropy < 30 keV cm^2 for CC
            ## Barnes+2018 use 30 & 60 keV cm^2, but define K0 excess by fitting formula to entropy profile
            ## Hahn+2017 use entropy at r=10 kpc < 40 keV cm^2 for CC
        },
    },
    ## Core pseudo entropy-based CC/NCC classification as in Truong+2019, Pearce+2021 (have not implemented)
    ## --> matches CHEERS selection criteria!
#     'pseudo_entropy':{
#         'classifications':['CC', 'Weak_CC', 'NCC'],
#         'classification_function':CC_classification_entropy,
#         'profile_to_use':three_d_profiles,
#         'profile_type':'all_core_profiles',
#         'num':4,
#         'kwargs':{
#             'prop_type':['K_mass', 'K_Lx_0.5-10.0keV_total', 'K_mass_bin_T_ne', 'K_Lx_0.5-10.0keV_total_bin_T_ne'],
#             'radial_bin':[0, 0, 0, 0],
#             'CC_entropy':[pnb.array.SimArray(22, units='keV cm**2'), pnb.array.SimArray(22, units='keV cm**2'), 
#                           pnb.array.SimArray(22, units='keV cm**2'), pnb.array.SimArray(22, units='keV cm**2')],
#             'NCC_entropy':[pnb.array.SimArray(150, units='keV cm**2'), pnb.array.SimArray(150, units='keV cm**2'), 
#                            pnb.array.SimArray(150, units='keV cm**2'), pnb.array.SimArray(150, units='keV cm**2')],  ## Lehle+2024
#         },
#     },
    ## Core entropy excess-based CC/NCC classification as in Barnes+2018 (have not implemented)
#     'entropy_excess':{
#         'classifications':['CC', 'Weak_CC', 'NCC'],
#         'classification_function':CC_classification_entropy,
#         'profile_to_use':three_d_profiles,
#         'profile_type':'all_core_profiles',
#         'num':4,
#         'kwargs':{
#             'prop_type':['K_mass', 'K_Lx_0.5-10.0keV_total', 'K_mass_bin_T_ne', 'K_Lx_0.5-10.0keV_total_bin_T_ne'],
#             'radial_bin':[0, 0, 0, 0],
#             'CC_entropy':[pnb.array.SimArray(22, units='keV cm**2'), pnb.array.SimArray(22, units='keV cm**2'), 
#                           pnb.array.SimArray(22, units='keV cm**2'), pnb.array.SimArray(22, units='keV cm**2')],
#             'NCC_entropy':[pnb.array.SimArray(150, units='keV cm**2'), pnb.array.SimArray(150, units='keV cm**2'), 
#                            pnb.array.SimArray(150, units='keV cm**2'), pnb.array.SimArray(150, units='keV cm**2')],  ## Lehle+2024
#         },
#     },
    ## Core ne-based CC/NCC classification as in Lehle+2024
    'ne':{
        'classifications':['CC', 'Weak_CC', 'NCC'],
        'classification_function':CC_classification_ne,
        'profile_to_use':three_d_profiles,
        'profile_type':'all_core_profiles',
        'num':1,
        'kwargs':{
            'prop_type':['ne_bin',],
            'radial_bin':[0,],
            'CC_ne':[pnb.array.SimArray(1.5e-2, units='cm**-3')],
            'NCC_ne':[pnb.array.SimArray(0.5e-2, units='cm**-3')],  ## Lehle+2024
        },
    },
    ## Physical concentration parameter-based CC/NCC classification as in Lehle+2024
    'Cphys':{
        'classifications':['CC', 'Weak_CC', 'NCC'],
        'classification_function':CC_classification_Cphys,
        'profile_to_use':two_d_profiles,
        'profile_type':'all_profiles',
        'num':1,
        'kwargs':{
            'prop_type':['Lx_0.5-10.0keV_total_enc',],   ## Is supposed to use soft X-ray band (0.5-5 keV)
            'inner_radius':[pnb.array.SimArray(40, units='kpc'),],
            'outer_radius':[pnb.array.SimArray(400, units='kpc'),],
            'CC_Cphys':[0.155],
            'NCC_Cphys':[0.075],  ## Lehle+2024
        },
    },
    ## Scaled concentration parameter-based CC/NCC classification as in Lehle+2024
    'Cscaled':{
        'classifications':['CC', 'Weak_CC', 'NCC'],
        'classification_function':CC_classification_Cscaled,
        'profile_to_use':two_d_profiles,
        'profile_type':'all_profiles',
        'num':1,
        'kwargs':{
            'prop_type':['Lx_0.5-10.0keV_total_enc',],   ## Is supposed to use soft X-ray band (0.5-5 keV)
            'inner_radius':[0.15,],   # Units of R/R500
            'outer_radius':[1,],
            'CC_Cscaled':[0.5],
            'NCC_Cscaled':[0.2],  ## Lehle+2024
        },
    },
}

classifications = []
## Add in classification types
if CLASSIFY_HALOS:
    for classification_type, classification_info in classifications_dict.items():
        for classification_label in classification_info['classifications']:
            for prop_type in classification_info['kwargs']['prop_type']:
                classifications.append(classification_label + '_' + classification_type + '_' + prop_type)


profiles['classifications'] = copy.deepcopy(classifications)
for classification in classifications:
    profiles[classification + '_id'] = np.array([], int)
    profiles[classification + '_profiles'] = {}
    
# classifications.append('all')
profiles['classifications'].append('all')
# print(profiles['classifications'])

print()
print('Classifications:')
print(profiles['classifications'])
print()
        
            


## Classify halos (e.g. CC/NCC)    
print()
print('Classifying halos as CC/NCC')
print()
if CLASSIFY_HALOS:
    for halo in profiles['all_id']:
        print('halo', halo)
        profile = profiles['all_profiles'][halo]['y']

        for classification_type, classification_info in classifications_dict.items():
            print(classification_type)

            try:
                profile_to_classify = classification_info['profile_to_use'][classification_info['profile_type']][halo]#['y']
            except:
                print('Profile', classification_info['profile_type'], 'not available for halo', halo)
                continue

            for ii in range(classification_info['num']):
                kwargs = {kwarg:classification_info['kwargs'][kwarg][ii] for kwarg in classification_info['kwargs'].keys()}
                
                this_classification = classification_info['classification_function'](profile_to_classify, kwargs)

#                 try:
#                     this_classification = classification_info['classification_function'](profile_to_classify, kwargs)

#                 except:
#                     print('No core profile for halo', halo, 'for CC/NCC', classifiation_type, 'classification')
#                     continue

                if this_classification is not None:
                    profiles[this_classification + '_' + 
                             classification_type + '_' + 
                             kwargs['prop_type'] + '_profiles'][halo] = copy.deepcopy(profiles['all_profiles'][halo])
                    profiles[this_classification + '_' + 
                             classification_type + '_' + 
                             kwargs['prop_type'] + '_id'] = np.append(profiles[this_classification + '_' + 
                                                                               classification_type + '_' + 
                                                                               kwargs['prop_type'] + '_id'], halo)
                
                
            
#             ## delete certain properties if want to save space

#             for prop_to_remove in x_props_to_remove:
#                 try:
#                     del profile_to_compute[profile_type_to_compute][halo]['x'][prop_to_remove]
#                     del profiles[
#                 except:
#                     print('Either', profile_type_to_compute, 'does not exist in', profile_to_compute, 'or', prop_to_remove, 'does not exist in', profile_type_to_compute)

#             for prop_to_remove in y_props_to_remove:
#                 try:
#                     del profile_to_compute[profile_type_to_compute][halo]['y'][prop_to_remove]
#                 except:
#                     print('Either', profile_type_to_compute, 'does not exist in', profile_to_compute, 'or', prop_to_remove, 'does not exist in', profile_type_to_compute)
            
                
        
        print()
        
        
        
        
print()
print('done')
print()
        
        
## delete certain properties if want to save space
print()
print('Removing properties to save space')
for classification in profiles['classifications']:
    for halo in profiles[classification+'_id']:
        for prop_to_remove in x_props_to_remove:
            try:
                del profiles[classification+'_profiles'][halo]['x'][prop_to_remove]
            except:
                print('Either', classification, 'or', prop_to_remove, 'does not exist')

        for prop_to_remove in y_props_to_remove:
            try:
                del profiles[classification+'_profiles'][halo]['y'][prop_to_remove]
                try:
                    del profiles[classification+'_profiles'][halo]['y']['log_'+prop_to_remove]
                except:
                    print('log(',+prop_to_remove+') does not exist')
            except:
                print('Either', classification, 'or', prop_to_remove, 'does not exist')
        
        
print('done')
print()



            

print()
print('Number of halos in each classification')
print()


for classification in profiles['classifications']:
    print(classification, ':', len(profiles[classification +'_id']))
    print()

print()





    
# Get only halo info that corresponds to profile ids
print('Get only halo info that corresponds to profile ids')
print()

profiles['halo_info'] = copy.deepcopy(profiles['temp_halo_info'])
# del profiles['temp_halo_info']
print('halo_info (same as temp_halo_info) keys:')
print(profiles['halo_info'].keys())
print()

bad_ids = np.setdiff1d(profiles['temp_halo_info']['halo_id'], profiles['all_id'])
bad_idxs = np.isin(profiles['temp_halo_info']['halo_id'], bad_ids).nonzero()[0]


for quan in profiles['halo_info'].keys():
    print(quan)
    if quan in ['group_ids', 'group_idx', 'group_idxs', 'rho_gas_scaling']:
        continue
    if quan in ['group_list', 'halo_id', 'N_lgal', 'CC_classification']:
        profiles['halo_info'][quan] = np.delete(profiles['halo_info'][quan], bad_idxs)
    else:
        detailed_quantity_dict = profiles['halo_info'][quan].keys()
        
        for d_quan in detailed_quantity_dict:
            if np.ndim(profiles['halo_info'][quan][d_quan]) > 0:
                profiles['halo_info'][quan][d_quan] = np.delete(profiles['halo_info'][quan][d_quan], bad_idxs)
    print('ok')

print()
print()
print()




## Bin profiles by halo properties

print('Bin profiles')
print()

bin_dict = {
#     'Unbinned':{
#         'prop':['M', '500'],
#         'prop_printable_name':r'$M_{500}\,/\,M_{\odot}$',
#         'bin_edges':pnb.array.SimArray(np.array([0, 1e16]), units='Msol'),
#     },
    'M500_bins':{
        'prop':['M', '500'],
        'prop_printable_name':r'$M_{500}\,/\,M_{\odot}$',
        'bin_edges':pnb.array.SimArray(np.array([1e11, 1e12, 10**12.5, 1e13, 10**13.5, 1e14, 1e15, 1e16]), units='Msol'),
    },
#     'thin_M500_bins':{
#         'prop':['M', '500'],
#         'prop_printable_name':r'$M_{500}\,/\,M_{\odot}$',
#         'bin_edges':pnb.array.SimArray(np.array([1e13, 10**13.1, 10**13.2, 10**13.3, 10**13.4, 10**13.5, 10**13.6, 10**13.7, 10**13.8, 
#                                                  10**13.9, 1e14, 10**14.2, 10**14.4, 10**14.6, 10**14.8, 1e15, 10**15.5, 1e16]), 
#                                        units='Msol'),
#     },
#     'Tspec_corr_bins':{
#         'prop':['T', 'spec_corr'],
#         'prop_printable_name':r'$k_B\,T_{spec,corr},[keV]$',
#         'bin_edges':pnb.array.SimArray(np.array([1e-1, 10**(-0.5), 10**(-0.2), 10**0, 10**1]), units='keV'),
#     },
#     'Tspec_corr_split':{
#         'prop':['T', 'spec_corr'],
#         'prop_printable_name':r'$k_B\,T_{spec,corr},[keV]$',
#         'bin_edges':pnb.array.SimArray(np.array([0, 1.1, 20]), units='keV')
#     },
}


profiles_prop_binned = {}
count = 0
for bin_prop_name_, bin_prop_info_ in bin_dict.items():
    print(bin_prop_name_)
    
    bin_prop_ = bin_prop_info_['prop']
    bin_prop_printable_name_ = bin_prop_info_['prop_printable_name']
    bin_edges_ = bin_prop_info_['bin_edges']
    bins_ = range(len(bin_edges_)-1)


    profiles_prop_binned[bin_prop_name_] = {i:copy.deepcopy(profiles) for i in bins_}
    profiles_prop_binned[bin_prop_name_]['bin_name'] = bin_prop_
    profiles_prop_binned[bin_prop_name_]['bin_printable_name'] = bin_prop_printable_name_
    profiles_prop_binned[bin_prop_name_]['bins'] = bins_
    profiles_prop_binned[bin_prop_name_]['bin_edges'] = bin_edges_
    profiles_prop_binned[bin_prop_name_]['metals_idx'] = metals_idx[CODE]
    profiles_prop_binned[bin_prop_name_]['sim'] = CODE + '(' + SIZE + ')'
    profiles_prop_binned[bin_prop_name_]['halo_info'] = copy.deepcopy(profiles['halo_info'])
    profiles_prop_binned[bin_prop_name_]['classifications'] = copy.deepcopy(profiles['classifications'])
    profiles_prop_binned[bin_prop_name_]['xaxis'] = copy.deepcopy(profiles['xaxis'])
    profiles_prop_binned[bin_prop_name_]['cosmo_props'] = copy.deepcopy(profiles['cosmo_props'])
    
    for i in bins_:
        for info_to_remove in ['xaxis', 'cosmo_props']:
            del profiles_prop_binned[bin_prop_name_][i][info_to_remove]

    # Find indices where elements should be inserted to maintain order
    # Need minus 1 to get correct bin
    bin_idxs_ = np.searchsorted(bin_edges_, profiles['halo_info'][bin_prop_[0]][bin_prop_[1]], side='left') - 1

    ## Check that this is correct!!!
    # go through each halo with the index of its bin for the current property
    for ii in range(len(bin_idxs_)):
        # Find the bins that this halo should be removed from
        if (bin_idxs_[ii]>=0 and bin_idxs_[ii]<len(bins_)):
            bins_to_remove_from_ = np.delete(bins_, bin_idxs_[ii])
        else:
            # If the halo's property falls outside any of the bins, remove it from all bins
            bins_to_remove_from_ = bins_
            
        # remove the halo/profile info from the bins that this halo does not fall in
        for jj in bins_to_remove_from_:
            halo_to_delete = profiles['all_id'][ii]
            idx_to_delete = np.where(profiles_prop_binned[bin_prop_name_][jj]['all_id']==halo_to_delete)
            
            profiles_prop_binned[bin_prop_name_][jj]['all_id'] = np.setdiff1d(profiles_prop_binned[bin_prop_name_][jj]['all_id'], halo_to_delete)
            profiles_prop_binned[bin_prop_name_][jj]['all_profiles'].pop(profiles['all_id'][ii])  # remove profile info
            
            for classification_ in profiles_prop_binned[bin_prop_name_]['classifications']:
#             for classification_ in profiles_prop_binned['classifications']:
                if halo_to_delete in profiles_prop_binned[bin_prop_name_][jj][classification_+'_profiles'].keys():
                    profiles_prop_binned[bin_prop_name_][jj][classification_+'_profiles'].pop(profiles['all_id'][ii])  # remove profile info
                profiles_prop_binned[bin_prop_name_][jj][classification_+'_id'] = np.setdiff1d(profiles_prop_binned[bin_prop_name_][jj][classification_+'_id'], halo_to_delete)

            # Remove halo info
            for quan in profiles_prop_binned[bin_prop_name_][jj]['halo_info'].keys():
                if quan in ['group_ids', 'group_idx', 'group_idxs', 'rho_gas_scaling', 'classifications']:
                    continue
                if quan in ['group_list', 'halo_id', 'N_lgal', 'CC_classification']:
                    profiles_prop_binned[bin_prop_name_][jj]['halo_info'][quan] = np.delete(profiles_prop_binned[bin_prop_name_][jj]['halo_info'][quan], idx_to_delete)
                else:
                    detailed_quantity_dict = profiles_prop_binned[bin_prop_name_][jj]['halo_info'][quan].keys()
                    for d_quan in detailed_quantity_dict:
                        if np.ndim(profiles_prop_binned[bin_prop_name_][jj]['halo_info'][quan][d_quan]) > 0:
                            profiles_prop_binned[bin_prop_name_][jj]['halo_info'][quan][d_quan] = np.delete(profiles_prop_binned[bin_prop_name_][jj]['halo_info'][quan][d_quan], idx_to_delete)


    ## Remove bins that have no halos in them
    for mm in range(len(bins_)):
        bin_ = bins_[mm]
        print(bin_, ':', len(profiles_prop_binned[bin_prop_name_][bin_]['all_id']))
        if len(profiles_prop_binned[bin_prop_name_][bin_]['all_id'])==0:
            del(profiles_prop_binned[bin_prop_name_][bin_])
            profiles_prop_binned[bin_prop_name_]['bins'] = np.delete(profiles_prop_binned[bin_prop_name_]['bins'], list(profiles_prop_binned[bin_prop_name_]['bins']).index(bin_))
    print()
    
    count += 1
    
print()



# Make sure halo info has correct units
print('Make sure halo info has correct units:')
print()

for bin_type, info_ in profiles_prop_binned.items():
    print(bin_type)
    for bin_ in info_['bins']:
        for quan in info_[bin_]['halo_info'].keys():
            if quan in ['rho_gas_scaling']:
                info_[bin_]['halo_info'][quan].units = info_[bin_]['temp_halo_info'][quan].units
            elif quan not in ['group_list', 'group_ids', 'group_idx', 'group_idxs', 'halo_id', 'N_lgal', 'CC_classification']:
                detailed_quantity_dict = info_[bin_]['halo_info'][quan].keys()
                for d_quan in detailed_quantity_dict:
                    info_[bin_]['halo_info'][quan][d_quan].units = info_[bin_]['temp_halo_info'][quan][d_quan].units
    print()
    
    
    
    
    
    
    
####### Calculating average profiles ###################################################
    
## Average together all radial profiles
print()
print()
print('Average together all radial profiles')
print()
print()


# def apply_func_to_dict_values(d, key_name, func, func_kwargs={}, extra_keys=[]):
#     for k, v in d.items():
#         if isinstance(v, dict):
#             apply_func_to_dict_values(v, key_name, func, 
#                                       func_kwargs=func_kwargs, extra_keys=extra_keys)
#         else:
#             if k == key_name:
#                 extra_kwargs = {extra_key:d[extra_key] for extra_key in extra_keys}
#                 func(v, **func_kwargs, **extra_kwargs)
                
# def remove_values_from_averaging(dict_of_values, kwargs):
#     ## Remove negative values from averaging by setting them to NaN
#     ## Check - for each radial bin - that at least 16% of halos have at least N particles --> if not, set all values in that bin to NaN
#     ## Also make sure there at at least M halos in each radial bin
#     ## When a radial bin meets the above two requirements, do I include the unresolved halos in the averaging still?
#     ## --> yes, since they will be weighted less in the averaging due to have fewer particles (so increased poisson error)    
    
#     print('Removing values from averaging')
#     for function, info in kwargs['functions'].items():
#         print()
#         print(function)
#         print()
#         apply_func_to_dict_values(dict_of_values, kwargs['key_name'], function, 
#                                   func_kwargs=info['func_kwargs'], extra_keys=info['extra_keys'])
                

# def exclude_negatives(values, replace_with=np.nan):
#     values[values<0] = replace_with
    
# def determine_resolved(counts, resolved_counts=1, iterate_over_bins=False):
#     if iterate_over_bins:
#         resolved_filter = np.full(counts.shape, True)
#         for ii in range(counts.shape[0]):
#             for jj in range(ii, counts.shape[0]):
#                 curr_filter = np.expand_dims(counts[jj], axis=0) >= resolved_counts
#                 resolved_filter[ii] = resolved_filter[ii] & curr_filter
#     else:
#         resolved_filter = counts >= resolved_counts
#     return resolved_filter

# def exclude_unresolved(values, counts=None, resolved_counts=1, iterate_over_bins=False, replace_with=np.nan):
#     resolved_filter = determine_resolved(counts, resolved_counts=resolved_counts, iterate_over_bins=iterate_over_bins)
#     values[~resolved_filter] = replace_with
# #     if iterate_over_bins:
# #         for ii in range(values.shape[0]):
# #             for jj in range(ii, values.shape[0]):
# # #                 if np.any(np.expand_dims(counts[jj], axis=0) < resolved_counts, axis=0):
# #                 unresolved_filter = np.expand_dims(counts[jj], axis=0) < resolved_counts
# #                 values[unresolved_filter] = replace_with
# #     else:
# #         values[counts<resolved_counts] = replace_with
    
# def exclude_less_than_length(values, physical_xbins=None, length=0, replace_with=np.nan):
#     for xbin in range(len(values)):
#         ## This way only excludes the innermost bin whose outer edge is less than the grav softening length
# #         if physical_xbins[xbin+1] <= length:
# #             values[xbin] = replace_with
#         ## This way excludes all bins that are thinner than the grav softening length
#         if physical_xbins[xbin+1] - physical_xbins[xbin] <= length:
#             values[xbin] = replace_with

            
# if INCLUDE_RESOLUTION_CUTS:
#     avg_kwargs['functions'] = {
#         exclude_negatives:{
#             'func_kwargs':{'replace_with':np.nan},
#             'extra_keys':[],
#         },
#         exclude_unresolved:{
#             'func_kwargs':{
#                 'resolved_counts':MIN_RESOLVED_COUNTS_PER_HALO, 
#                 'iterate_over_bins':PROFILES_ARE_2D_BUT_RES_FILE_IS_3D, 
#                 'replace_with':np.nan,
#             },
#             'extra_keys':['counts'],
#         },
#         exclude_less_than_length:{
#             'func_kwargs':{
#                 'length':GRAV_SOFTENING_LENGTH.in_units(XSCALE_UNITS_DICT['Physical']),
#                 'replace_with':np.nan
#             },
#             'extra_keys':['physical_xbins'],
#         },
#     }
    
    
            
avg_kwargs = {
    'median':{
        'quantiles':[0.25,0.5,0.75],#[0.16,0.5,0.84],
        'interpolate':True,
    },
    'key_name':'values',
    'functions':{},
}
        

resolution_kwargs = {
    'ref_bins':REF_BINS,
    'min_resolved_counts':MIN_RESOLVED_COUNTS_PER_BIN,
    'min_total_halos':MIN_TOTAL_HALOS,
    'min_resolved_halos':MIN_RESOLVED_HALOS,
    'min_frac_resolved_halos':MIN_FRAC_RESOLVED_HALOS,
}


# def exclude_bins_with_too_few_resolved_halos(profiles, ref_bins=None, physical_ref_bins=None, min_resolved_counts=1,
#                                              min_total_halos=1, min_resolved_halos=0, min_frac_resolved_halos=0, 
#                                              grav_softening_length=0):
def exclude_bins_with_too_few_resolved_halos(profiles, ref_bins=None, min_resolved_counts=1,
                                             min_total_halos=1, min_resolved_halos=0, min_frac_resolved_halos=0):
    profile_length = len(ref_bins)-1
    
    resolution_filter = np.zeros(profile_length, bool)
    resolution_info = {}

    for ref_bin in range(profile_length):
        ref_bin_lo = ref_bins[ref_bin]
        ref_bin_hi = ref_bins[ref_bin+1]
#         ref_bin_middle = (ref_bin_lo + ref_bin_hi)/2.

        resolution_filter[ref_bin] = False
        resolution_info[ref_bin] = {
            'total_halos':0,
            'resolved_halos':0,
        }

        for halo, info in profiles.items():
            rbins = info['xbins'].copy()
            counts = info['counts'].copy()
            n_rbins = len(rbins)-1
            
            resolved_filter = determine_resolved(counts, resolved_counts=min_resolved_counts, 
                                                 iterate_over_bins=PROFILES_ARE_2D_BUT_RES_FILE_IS_3D)
            

            for r_idx in range(n_rbins):
                
                ## Weight factor for geometric overlap of radial bin of observed halo with reference radial bin
                overlap = calc_bin_overlap(rbins[r_idx], rbins[r_idx+1], ref_bin_lo, ref_bin_hi)
                
                if overlap < SMALL_OVERLAP:
                    overlap = 0
                    
                resolution_info[ref_bin]['total_halos'] += overlap
                
#                 if counts[r_idx] >= min_resolved_counts:
                if resolved_filter[r_idx]:
                    resolution_info[ref_bin]['resolved_halos'] += overlap

                    
        try:
            resolution_info[ref_bin]['frac_resolved_halos'] = resolution_info[ref_bin]['resolved_halos']/resolution_info[ref_bin]['total_halos']
        except:
            resolution_info[ref_bin]['frac_resolved_halos'] = 0

        if (resolution_info[ref_bin]['total_halos'] >= min_total_halos and 
            resolution_info[ref_bin]['resolved_halos'] >= min_resolved_halos and 
            resolution_info[ref_bin]['frac_resolved_halos'] >= min_frac_resolved_halos):
            resolution_filter[ref_bin] = True
            
    return resolution_filter, resolution_info
    

                                             

profiles_binned = copy.deepcopy(profiles_prop_binned)
profiles_stats = {}
for bin_type, info in profiles_binned.items():
    print(bin_type)
    print()
    print()
    
    
    ## Storing all profiles in each binning property, bin and classification
    
    profiles_stats[bin_type] = {}
    
    for bin_ in info['bins']:
        print('Bin', bin_)
        print()
        
        profiles_stats[bin_type][bin_] = {}
        
        for classification in info[bin_]['classifications']:
            print(classification)
            
#             if len(dict_of_profiles.keys()) < 0.5:
            if len(info[bin_][classification+'_id']) < 0.5:
                print()
                print('No halos in classification')
                print()
                print()
                continue
            
            profiles_stats[bin_type][bin_][classification+'_profiles'] = {
                'avg':{},  ## avg
                'limits':{},  ## limits
            }
            
            dict_of_profiles = {}
    
            for halo in info[bin_][classification+'_id']:
                halo_idx = list(info[bin_][classification+'_id']).index(halo)
                print('Saving halo', halo, 'profile, idx', halo_idx)
                

                # Add to saved profiles
                for prop in info[bin_][classification+'_profiles'][halo]['y'].keys():
                    
                    if prop not in dict_of_profiles.keys():
                        dict_of_profiles[prop] = {}

                    dict_of_profiles[prop][halo] = {
                        'xbins':info[bin_][classification+'_profiles'][halo]['x']['rbins'],
                        'physical_xbins':info[bin_][classification+'_profiles'][halo]['x']['physical_rbins'],
                        'values':info[bin_][classification+'_profiles'][halo]['y'][prop],
#                         'counts':info[bin_][classification+'_profiles'][halo]['y']['n'],
                        'counts':res_profiles['all_profiles'][halo]['y']['n'],
                    }

                    if USE_ERRORS:
                        errors = 1./calc_poisson_err(info[bin_][classification+'_profiles'][halo]['y']['n'], log=False)
                    else:
                        errors = np.ones(len(info[bin_][classification+'_profiles'][halo]['y'][prop]))
                    dict_of_profiles[prop][halo]['errs'] = errors
                    
            print()
            print()
    
            
            ######## Calculate relevant statistics ###########################################################
            
            ## Get filter to apply to results to remove bins that do not meet resolution criteria
#             print()
#             print(dict_of_profiles.keys())
#             print(dict_of_profiles[list(dict_of_profiles.keys())[0]])
#             print()
            resolution_filter, resolution_info = exclude_bins_with_too_few_resolved_halos(dict_of_profiles['n'], **resolution_kwargs)
            print()
            print('resolution_filter:')
            print(resolution_filter)
            print()
            print('resolution_info:')
            pprint.pp(resolution_info)
            print()
#             print()
#             print(dict_of_profiles.keys())
#             print(dict_of_profiles[list(dict_of_profiles.keys())[0]])
#             print()
            
    
#             prop0 = list(dict_of_profiles.keys())[0]
#             halo0 = list(dict_of_profiles[prop0].keys())[0]
#             print('prop0:', prop0)
#             print('halo0:', halo0)
# #             halo0 = prop0_info[list(prop0_info.keys())[0]]
#             print("dict_of_profiles[prop0][halo0]['values'].shape:", dict_of_profiles[prop0][halo0]['values'].shape)
#             print("dict_of_profiles[prop0][halo0]['values']:")
#             print(dict_of_profiles[prop0][halo0]['values'])
#             print()
    
#             ## Remove specified values from averaging
#             remove_values_from_averaging(dict_of_profiles, avg_kwargs) ## WHY NO WORK??? -> works :)
            
# #             print('dict_of_profiles:')
# #             print(dict_of_profiles)
# #             print()
            
#             prop0 = list(dict_of_profiles.keys())[0]
#             halo0 = list(dict_of_profiles[prop0].keys())[0]
#             print('prop0:', prop0)
#             print('halo0:', halo0)
#             print("dict_of_profiles[prop0][halo0]['values'].shape:", dict_of_profiles[prop0][halo0]['values'].shape)
#             print("dict_of_profiles[prop0][halo0]['values']:")
#             print(dict_of_profiles[prop0][halo0]['values'])
#             print()
            

            for prop, prop_values in dict_of_profiles.items():
                print('bin', bin_,',', prop)
                
#                 if 'Z' not in prop:
#                 if prop != 'Z_median_ratio':
#                     continue
                
#                 print(prop_values)
#                 print()

#                 print('prop_values[values].shape:', prop_values['values'].shape)
#                 print('prop_values[values]:')
#                 print(prop_values['values'])
#                 print()

#                 ## Remove specified values from averaging
#                 prop_values = remove_values_from_averaging(prop_values, avg_kwargs) ## WHY NO WORK???
                
                prop_units = profiles['all_profiles'][profiles['all_id'][0]]['y'][prop].units
#                 print('prop_values[values].shape:', prop_values['values'].shape)
#                 print('prop_values[values]:')
#                 print(prop_values['values'])
#                 print()
#                 print()
#                 print(prop_values)
                
                if AVG_TYPE.lower() == 'mean':
                    ## Average profiles according to Mernier+2017 equation 3
                    result = weighted_mean_of_profiles(info['xaxis']['all']['bins'], prop_values)
#                     print('result.shape:', result.shape)
#                     print('result:')
#                     print(result)
#                     print()
                    averaged_profile = result[:,0]
                    std_profile = result[:,1]
#                     print('averaged_profile.shape:', averaged_profile.shape)
#                     print('averaged_profile:')
#                     print(averaged_profile)
#                     print()
#                     print('std_profile.shape:', std_profile.shape)
#                     print('std_profile:')
#                     print(std_profile)
#                     print()
                    if prop.startswith('log'):
                        averaged_limits = [np.log10(10**averaged_profile - 10**std_profile), 
                                           np.log10(10**averaged_profile + 10**std_profile)]
                    else:
                        averaged_limits = [averaged_profile - std_profile, 
                                           averaged_profile + std_profile]
                elif AVG_TYPE.lower() == 'median':
                    ## Same thing as above, but with weighted quantiles
                    result = weighted_quantiles_of_profiles(info['xaxis']['all']['bins'], prop_values, 
                                                            calc_kwargs=avg_kwargs['median'])
#                     print('final result.shape:', result.shape)
#                     print('final result:')
#                     print(result)
#                     print()
                    closest_idx_to_middle = np.argmin(np.abs(np.array(avg_kwargs['median']['quantiles']) - 0.5))
#                     print('closest_idx_to_middle:', closest_idx_to_middle)
                    averaged_profile = result[:,closest_idx_to_middle]
                    averaged_limits = [result[:,0], result[:,-1]]
                else:
                    raise Exception('avg_type must be mean or median')
            
            
                ## Remove bins that do not meet resolution criteria
#                 averaged_profile = averaged_profile[resolution_filter]
#                 averaged_limits = [limit[resolution_filter] for limit in averaged_limits]
                averaged_profile[~resolution_filter] = np.nan
                for limit in averaged_limits:
                    limit[~resolution_filter] = np.nan
            
                
                averaged_profile = pnb.array.SimArray(averaged_profile, units=prop_units)
                averaged_limits = {
                    'lo':pnb.array.SimArray(averaged_limits[0], units=prop_units), ## lo
                    'hi':pnb.array.SimArray(averaged_limits[1], units=prop_units), ## hi
                }
            
                profiles_stats[bin_type][bin_][classification+'_profiles']['avg'][prop] = averaged_profile
                profiles_stats[bin_type][bin_][classification+'_profiles']['limits'][prop] = averaged_limits
            

                print('bin', bin_,',', prop, 'done')
                print()

        
            print()
            print()
        print()
        print()
    print()
    print()
    
    
profiles_binned['all_profiles'] = copy.deepcopy(profiles)



## Save extended profiles and profiles stats to file
print('Saving extended profiles and profiles stats to file')

# defaults to allow_pickle=True
save_file = os.path.join(PROFILES_DIR, SNAP + '_' + PROFILES_SAVE_NAME + '-' + SUFFIX + '_binned_extended' + '.pkl')
print(save_file)
save_object_with_dill(profiles_binned, save_file, protocol=3)

save_file = os.path.join(PROFILES_DIR, SNAP + '_' + PROFILES_SAVE_NAME + '-' + SUFFIX + '_binned_stats' + '.pkl')
print(save_file)
save_object_with_dill(profiles_stats, save_file, protocol=3)

save_file = os.path.join(PROFILES_DIR, SNAP + '_' + PROFILES_SAVE_NAME + '-' + SUFFIX + '_parameters' + '.pkl')
print(save_file)
save_object_with_dill(parameters, save_file, protocol=3)

print()
print('DONE')