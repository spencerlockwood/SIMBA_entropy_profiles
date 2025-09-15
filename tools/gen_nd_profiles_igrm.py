import numpy as np
import pynbody as pnb
import os
import h5py
import XIGrM.gas_properties as g_p
from XIGrM.gas_properties import m_p, k_B, default_elements
import XIGrM.X_properties as x_p
import XIGrM.halo_analysis as h_a
import XIGrM.prepare_pyatomdb as ppat
from pynbody.analysis.profile import Profile, VerticalProfile, InclinedProfile, QuantileProfile
import argparse
import dill

from astro_constants import NA_no_units, NA, kB_no_units, kB, mH, mp_no_units, mp, mu_e, mu, G_no_units, G


parser = argparse.ArgumentParser(description="Create radial profiles of general X-ray properties of the IGrM.")
parser.add_argument('--code', action='store', type=str, required=True, 
                   help='Simulation type (currently either just Simba or Simba-C)')
parser.add_argument('--redshift', action='store', type=float, required=True, 
                   help='Simulation redshift')
parser.add_argument('--snap', action='store', type=str, required=True, 
                   help='Path to snapshot to analyze')
parser.add_argument('--ahf_basename', action='store', type=str, required=True, 
                   help='Path+basename to AHF files for snapshot')
parser.add_argument('--nthreads', action='store', type=int, required=True, 
                   help='Number of threads to use to calculate ahf halos')
parser.add_argument('--prop', action='store', type=str, required=True, 
                   help='Path to snapshot host prop file')
parser.add_argument('--save', action='store', type=str, required=True, 
                   help='Path to file to save profiles in')
parser.add_argument('--filter', action='store', type=str, required=True, 
                   help="Filter to place on each halo: 'Sphere' (for 3d or 2d profiles), 'Cylinder' (for 2d profiles)")
parser.add_argument('--profile_type', action='store', type=str, required=True, 
                   help="Type of profile to create ('Profile', 'VerticalProfile', 'InclinedProfile', 'QuantileProfile')")
parser.add_argument('--ndim', action='store', type=int, required=True, 
                   help='Number of dimensions for profile, ie. 2 is projected, 3 is spherically averaged (can only be 2 or 3)')
parser.add_argument('--weight_by', action='store', type=str, required=True, 
                   help='Property to weight averages in each radial bin by, defaults to mass')
parser.add_argument('--xscale', action='store', type=str, required=True, 
                   help="Units of radial axis ('R500', 'R200', or 'Physical')")
parser.add_argument('--temp_cut', action='store', type=str, required=True, 
                   help='Temperature in K above which to keep gas particles for the hot diffuse gas')
parser.add_argument('--nh_cut', action='store', type=str, required=True, 
                   help='Hydrogen number density in cm**-3 below which to keep gas particles for the hot diffuse gas')
parser.add_argument('--N_lgal', action='store', type=int, required=True, 
                   help='Minimum number of luminous (resolved) galaxies in a halo to be considered a group')
parser.add_argument('--log_M500_cut', action='store', type=float, required=True, 
                   help='Only keep groups above this value (in Msun)')
args = parser.parse_args()



profile_types = {
    'Profile':Profile,
    'VerticalProfile':VerticalProfile,
    'InclinedProfile':InclinedProfile,
    'QuantileProfile':QuantileProfile,
}



code = args.code
sim_z = args.redshift
data_file = args.snap
ahf_basename = args.ahf_basename
nthreads = args.nthreads
host_prop_file = args.prop
save_file = args.save + '_temp_cut='+args.temp_cut+'_nh_cut='+args.nh_cut+'.pkl'
halo_filt = args.filter
profile_type = profile_types[args.profile_type]
profile_ndims = args.ndim
weight_by = args.weight_by
xscale = args.xscale
temp_cut = args.temp_cut
nh_cut = args.nh_cut
N_lgal = args.N_lgal
log_M500_cut_val = args.log_M500_cut
    
    
    
    

def save_object_with_dill(obj, filename):
    with open(filename, 'wb') as f:  # Overwrites any existing file.
        dill.dump(obj, f, dill.HIGHEST_PROTOCOL)
    
    

## Function to turn mass fraction Z to abundance nX/nH for any element
def mass_fraction_to_relative_abundance(mass_fraction, elements=default_elements):
    ## Number density of each element relative to H number density
    
    atomicNumbers = ppat.elsymbs_to_z0s(elements)
    aMasses = ppat.get_atomic_masses(atomicNumbers)
    
    relative_massfraction = mass_fraction / mass_fraction[:, metals_idx_xigrm['H']].reshape(-1,1) # mass fraction relative to hydrogen
    relative_abundance = relative_massfraction * aMasses[0] / aMasses
    relative_abundance = pnb.array.SimArray(relative_abundance)
    relative_abundance.units = ''
    return relative_abundance
    



###### Profile Properties #####################################################################

@Profile.profile_property
def prop_weighted_avg(self, prop, weight):
    with self.sim.immediate_mode:
        pprop_ = self.sim[prop]
        pprop = pprop_.view(np.ndarray)
        pweight_ = self.sim[weight]
        pweight = pweight_.view(np.ndarray)

    result = pnb.array.SimArray(np.zeros((self.nbins,) + pprop.shape[1:]), units=pprop_.units)
    for i in range(self.nbins):
        if len(self.binind[i]) == 0:
            result[i] = np.nan
        else:
            result[i] = np.average(pprop[self.binind[i]], weights=pweight[self.binind[i]], axis=0)
    return result

@Profile.profile_property
def prop_median(self, prop):
    with self.sim.immediate_mode:
        pprop_ = self.sim[prop]
        pprop = pprop_.view(np.ndarray)
    
    median_prop = pnb.array.SimArray(np.zeros((self.nbins,) + pprop.shape[1:]), units=pprop_.units)
    for i in range(self.nbins):
        if len(self.binind[i]) == 0:
            median_prop[i] = np.nan
        else:
            median_prop[i] = np.nanmedian(pprop[self.binind[i]], axis=0)
    return median_prop



## Geometric Properties

@Profile.profile_property
def volume(self):
    return self._binsize



## Total metallicity

@Profile.profile_property
def Ztot_weight(self, weight):
    with self.sim.immediate_mode:
        pmf_ = self.sim['metals'][:,metals_idx['Z']]
        pmf = pmf_.view(np.ndarray)
        pweight_ = self.sim[weight]
        pweight = pweight_.view(np.ndarray)

    metallicity = pnb.array.SimArray(np.zeros(self.nbins), units='')
    for i in range(self.nbins):
        if len(self.binind[i]) == 0:
            metallicity[i] = np.nan
        else:
            metallicity[i] = np.average(pmf[self.binind[i]], weights=pweight[self.binind[i]], axis=0)
    return metallicity



## Number density/abundance profiles

@Profile.profile_property
def Ne(self):
    with self.sim.immediate_mode:
        pvol_ = self.sim['volume'].in_units('cm**3')
        pvol = pvol_.view(np.ndarray)
        pne_ = self.sim['ne'].in_units('cm**-3')
        pne = pne_.view(np.ndarray)
    
    Ne_ = pnb.array.SimArray(np.zeros(self.nbins), units=pvol_.units * pne_.units)
    for i in range(self.nbins):
        if len(self.binind[i]) == 0:
            Ne_[i] = 0
        else:
            Ne_[i] = np.nansum(pvol[self.binind[i]] * pne[self.binind[i]], axis=0)#.sum()
    return Ne_

@Profile.profile_property
def Nh(self):
    with self.sim.immediate_mode:
        pvol_ = self.sim['volume'].in_units('cm**3')
        pvol = pvol_.view(np.ndarray)
        pnh_ = self.sim['nh'].in_units('cm**-3')
        pnh = pnh_.view(np.ndarray)
    
    Nh_ = pnb.array.SimArray(np.zeros(self.nbins), units=pvol_.units * pnh_.units)
    for i in range(self.nbins):
        if len(self.binind[i]) == 0:
            Nh_[i] = 0
        else:
            Nh_[i] = np.nansum(pvol[self.binind[i]] * pnh[self.binind[i]], axis=0)#.sum()
    return Nh_

@Profile.profile_property
def NX(self):
    with self.sim.immediate_mode:
        pnXnH_ = self.sim['nX/nH'].in_units('')
        pnXnH = pnXnH_.view(np.ndarray)
        pnH_ = self.sim['nh'].in_units('cm**-3')
        pnH = pnH_.view(np.ndarray).reshape(-1,1)    
        pvol_ = self.sim['volume'].in_units('cm**3')
        pvol = pvol_.view(np.ndarray).reshape(-1,1)
    
    n_part, n_elements = pnXnH.shape
    NX_ = pnb.array.SimArray(np.zeros((self.nbins, n_elements)), units=pnXnH_ * pnH_.units * pvol_.units)
    for i in range(self.nbins):
        if len(self.binind[i]) == 0:
            NX_[i] = 0
        else:
            NX_[i] = np.nansum(pnXnH[self.binind[i]] * pnH[self.binind[i]] * pvol[self.binind[i]], axis=0)
    return NX_



## Thermal Properties

@Profile.profile_property
def Lx_band(self, band):
    '''
    Luminosity within each bin, calculated via dirrectly summing.
    '''
    with self.sim.immediate_mode:
        pLx_ = self.sim['Lx_' + band].in_units('erg s**-1')
        pLx = pLx_.view(np.ndarray)
        
    Lx = pnb.array.SimArray(np.zeros(self.nbins), units=pLx_.units)
    for i in range(self.nbins):
        if len(self.binind[i]) == 0:
            Lx[i] = 0.0 #np.nan
        else:
            Lx[i] = np.nansum(pLx[self.binind[i]])#.sum()
    return Lx

@Profile.profile_property
def T_spec(self):
    kT = pnb.array.SimArray(np.zeros(self.nbins), units='keV')
    subsnap = self.sim
    for i in range(self.nbins):
        if len(self.binind[i]) == 0:
            kT[i] = np.nan
        else:
            kT[i] = x_p.cal_tspec(hdgas=subsnap[self.binind[i]], cal_f=tspec_cal_file, datatype='gizmo')
    
    T = (kT/kB).in_units('K')
    return T



def gas_profile(snapshot, center, rbins, xray_bands, 
                profile_type=Profile, profile_ndims=3, weight_by='mass'):
    tx = pnb.transformation.inverse_translate(snapshot, center)
    boxsize = snapshot.properties['boxsize']
    
#     log_rbins = np.log10(rbins)
#     bin_centres = (rbins[:-1] + rbins[1:])/2.
#     log_bin_centres = (log_rbins[:-1] + log_rbins[1:])/2.
    
    
    result = {
        'x':{},
        'y':{},
    }
    with tx:
        original_pos = snapshot['pos'].copy()
        
        # Correct the position of particles crossing the box periodical boundary.
#         h_a.correct_pos(snapshot['pos'], boxsize)
        for i in range(3): # Correct the position of particles crossing the box periodical boundary.
            index1, = np.where(snapshot['pos'][:, i] < -boxsize/2)
            snapshot['pos'][index1, i] += boxsize
            index2, = np.where(snapshot['pos'][:, i] > boxsize/2)
            snapshot['pos'][index2, i] -= boxsize


        ## By default, for ndim=2, particles are projected onto x-y axis (ie. with z-axis),
        ## so disc filter is fine
        pg = profile_type(snapshot, ndim=profile_ndims, bins=rbins, weight_by=weight_by)
        
        
        
        ## x-axis values
        
        result['x']['rbins'] = rbins
#         result['x']['log_rbins'] = log_rbins
#         result['x']['bin_centres'] = bin_centres
#         result['x']['log_bin_centres'] = log_bin_centres
        
        
        
        ## y-axis values
        
        result['y']['n'] = pnb.array.SimArray(pg['n'], units='') # KEEP, this is the number of particles, we need this
        # MAKE A FUNCTION FO
        result['y']['volume'] = pg['volume']
        
        result['y']['Ne'] = pg['Ne']
        result['y']['Nh'] = pg['Nh']
#         result['y']['NX'] = pg['NX']
        
        result['y']['mass'] = pg['mass']
        for xband in xray_bands:
            result['y']['Lx_' + xband] = Lx_band(pg, xband)
            
        median_props = ['mass'] + ['Lx_' + xband for xband in xray_bands]
        for prop in median_props:
            result['y'][prop + '_median'] = prop_median(pg, prop)
#             result['y'][prop + '_mean'] = prop_weighted_avg(pg, prop, prop)
        
#         props = ['Z', 'nX/nH', 'T', 'ne', 'K', 'P'] + ['tcool_' + xband for xband in xray_bands]
        props = ['Z', 'T', 'ne', 'K', 'P'] + ['tcool_' + xband for xband in xray_bands]
#         weights = ['mass', 'volume', 'ne', 'T'] + ['Lx_' + xband for xband in xray_bands]
        weights = ['mass', 'volume'] + ['Lx_' + xband for xband in xray_bands]
        for prop in props:
            result['y'][prop + '_median'] = prop_median(pg, prop)
            for weight in weights:
                if prop != weight:
                    result['y'][prop + '_' + weight] = prop_weighted_avg(pg, prop, weight)
                    
        for weight in weights:
            result['y']['Ztot_' + weight] = Ztot_weight(pg, weight)
        
        result['y']['T_spec'] = pg['T_spec']
        
#         result['y']['tdyn'] = pg['dyntime']
        
        snapshot['pos'] = original_pos
    
    return result



################################################# Start of Analysis #########################################################

## Call simulation
print()
print('Reading in simulation snapshot')
s = pnb.load(data_file)
s.physical_units()
print('done')
print()


## Cosmological Parameters (figure out how to change for different redshifts)
age = pnb.array.SimArray(pnb.analysis.cosmology.age(s, z=None, unit='Gyr'), units='Gyr')
z = pnb.array.SimArray(s.properties['Redshift'], units='1')
rho_crit = pnb.array.SimArray(pnb.analysis.cosmology.rho_crit(s, z=None, unit='Msol kpc^-3'), units='Msol kpc**-3')
h0 = pnb.array.SimArray(s.properties['h'], units='1')
omL0 = pnb.array.SimArray(s.properties['omegaL0'], units='1')
omB0 = pnb.array.SimArray(s.properties['omegaB0'], units='1')
omM0 = pnb.array.SimArray(s.properties['omegaM0'], units='1')

## Scaling density for gas (universal for all halos)
f_b = pnb.array.SimArray(omB0/omM0, units='1')  # cosmic baryon fraction
rho_gas_scaling = rho_crit * f_b


## X-ray emissivity files with correct bands

# Aviv
# _05_18_band_file = '/project/b/babul/aspadawe/data/atomdb/Data/0.5-1.8keV_emissivity_all_elements_1000bins.hdf5'
# _05_2_band_file = '/project/b/babul/aspadawe/data/atomdb/Data/0.5-2.0keV_emissivity_all_elements_1000bins.hdf5'
# _01_24_band_file = '/project/b/babul/aspadawe/data/atomdb/Data/0.1-2.4keV_emissivity_all_elements_1000bins.hdf5'
# _05_7_band_file = '/project/b/babul/aspadawe/data/atomdb/Data/0.5-7.0keV_emissivity_all_elements_1000bins.hdf5'
_05_10_band_file = '/project/b/babul/aspadawe/data/atomdb/Data/0.5-10.0keV_emissivity_all_elements_1000bins.hdf5'

# xray_bands = ['0.5-1.8keV_total', '0.5-2.0keV_total', '0.1-2.4keV_total', '0.5-7.0keV_total', '0.5-10.0keV_total', 
#               '0.5-1.8keV_cont', '0.5-2.0keV_cont', '0.1-2.4keV_cont', '0.5-7.0keV_cont', '0.5-10.0keV_cont']

# xray_bands = ['0.5-10.0keV_total', '0.5-10.0keV_cont']
xray_bands = ['0.5-10.0keV_total']

# xray_nbins = [1000, 1000, 1000, 1000, 1000, 
#               1000, 1000, 1000, 1000, 1000]
# xray_nbins = [1000, 1000]
xray_nbins = [1000]
              
              
# xray_band_limits = [[0.5,1.8], [0.5,2.0], [0.1,2.4], [0.5,7.0], [0.5,10.0], 
#                     [0.5,1.8], [0.5,2.0], [0.1,2.4], [0.5,7.0], [0.5,10.0]]
# xray_band_limits = [[0.5,10.0], [0.5,10.0]]
xray_band_limits = [[0.5,10.0]]

# xray_emissivity_files = [_05_18_band_file, _05_2_band_file, _01_24_band_file, _05_7_band_file, _05_10_band_file, 
#                          _05_18_band_file, _05_2_band_file, _01_24_band_file, _05_7_band_file, _05_10_band_file]
# xray_emissivity_files = [_05_10_band_file, _05_10_band_file]
xray_emissivity_files = [_05_10_band_file]

# xray_modes = ['total', 'total', 'total', 'total', 'total', 
#              'cont', 'cont', 'cont', 'cont', 'cont']
# xray_modes = ['total', 'cont']
xray_modes = ['total']



## Find proper pytspec calibration file. See pytspec documentation for details.
# cal_dat_dir = '/home/b/babul/aspadawe/project/data/tspec_calibration/'
cal_dat_dir = '/project/b/babul/aspadawe/data/tspec_calibration/'
cal_dat_files = os.listdir(cal_dat_dir)
cal_redshift = []
for calfile in cal_dat_files:
    cal_redshift += [eval(calfile[10:17])]
calfile_idx = np.abs(np.array(cal_redshift) - s.properties['z']).argmin()
tspec_cal_file = (cal_dat_dir + cal_dat_files[calfile_idx]) ## good



## Metal indices
if code.lower() == 'simba':
    # Simba
    metals_idx = {'Z':0, 'He': 1, 'C': 2, 'N': 3, 'O': 4, 'Ne': 5, 'Mg': 6, 'Si': 7, 'S': 8, 'Ca': 9, 'Fe': 10}
    metals_idx_xigrm = {'H':0, 'He': 1, 'C': 2, 'N': 3, 'O': 4, 'Ne': 5, 'Mg': 6, 'Si': 7, 'S': 8, 'Ca': 9, 'Fe': 10}
    elems = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe']
elif code.lower() == 'simba-c':
    # Simba-C
    metals_idx = {'Z':0, 'H':1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al':13,
                  'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24,
                  'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30}
    metals_idx_xigrm = {'H':0, 'He': 1, 'Li': 2, 'Be': 3, 'B': 4, 'C': 5, 'N': 6, 'O': 7, 'F': 8, 'Ne': 9, 'Na': 10, 'Mg': 11, 
                        'Al':12, 'Si': 13, 'P': 14, 'S': 15, 'Cl': 16, 'Ar': 17, 'K': 18,'Ca': 19, 'Sc': 20, 'Ti': 21, 'V': 22, 'Cr': 23,
                      'Mn': 24, 'Fe': 25, 'Co': 26, 'Ni': 27, 'Cu': 28, 'Zn': 29}
    elems=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
           'P', 'S', 'Cl', 'Ar', 'K','Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']#, 'Ga', 'Ge']\n",
else:
    raise SystemExit('Code is not Simba or Simba-C')



## Calculate gas particle properties required for isolating hot diffuse gas (IGrM)
print()
print('Calculating gas properties reguired for isolating IGrM')

s.gas['u'] /= pnb.array.SimArray(s.properties['a'], units='1')

# Hydrogen mass fraction
if code.lower() == 'simba':
    s.gas['X_H'] = 1-s.gas['metals'][:,metals_idx['Z']]-s.gas['metals'][:,metals_idx['He']]
elif code.lower() == 'simba-c':
    s.gas['X_H'] = s.gas['metals'][:,metals_idx['H']]

s.gas['nh'] = g_p.nh(s) # Hydrogen number density
s.gas['ElectronAbundance'].units = '1'
s.gas['T'] = g_p.temp(s) # From internal energy to temperature
s.gas['temp'] = s.gas['T']

print('DONE')
print()


## Get hot diffuse gas using provided cuts
print()
print('Isolating IGrM (temperate, hydrogen number density, and decoupled wind particle cuts)')
## default: temp_cut = '5e5 K'
## default: nh_cut = '0.13 cm**-3'
temp_filt = pnb.filt.HighPass('T', temp_cut)
density_filt = pnb.filt.LowPass('nh', nh_cut)
wind_filt = ~pnb.filt.HighPass('DelayTime', 0)  # Cut out hot decoupled wind particles
total_filt = temp_filt & density_filt & wind_filt
hotdiffusegas = s.gas[total_filt]
print('DONE')
print()



print()
print('Calculating IGrM properties')

hotdiffusegas['ne'] = hotdiffusegas['ElectronAbundance'] * hotdiffusegas['nh'].in_units('cm**-3') # eletron number density
hotdiffusegas['volume'] = hotdiffusegas['mass']/hotdiffusegas['rho'] # effective volume of gas particles
hotdiffusegas['K'] = kB * hotdiffusegas['T'] * hotdiffusegas['ne']**(-2,3) # entropy
hotdiffusegas['P'] = kB * hotdiffusegas['T'] * hotdiffusegas['ne'] # Pressure


# Calculate mass fraction of the included elements. The first column must be Hydrogen mass fraction
print()
print('Calculating mass fractions')
if code.lower() == 'simba':
    # Simba
    hotdiffusegas['mass_fraction'] = np.zeros(hotdiffusegas['metals'].shape)
    hotdiffusegas['mass_fraction'][:, 0] = hotdiffusegas['X_H']
    hotdiffusegas['mass_fraction'][:, 1:] = hotdiffusegas['metals'][:, 1:]
elif code.lower() == 'simba-c':
    # Simba-C
    hotdiffusegas['mass_fraction'] = hotdiffusegas['metals'][:,1:-3]

hotdiffusegas['Z'] = hotdiffusegas['mass_fraction']

print('DONE')
print()

## Gas property initialization
print()
print('Initializing Lx')
Emission_type = []
for xband_ in xray_bands:
    Emission_type.append('Lx_' + xband_)
for i in Emission_type:
    hotdiffusegas[i] = 0
    hotdiffusegas[i].units = 'erg s**-1'

print('DONE')
print()


## Abundances relative to solar abundance, don't have to be in gizmo format. 

## New method to calculate abundances and luminosities that uses less memory
print()
# print('Calculating abundances, cooling rates, & X-ray luminosities')
print('Calculating abundances & X-ray luminosities')

nsteps = 500000
tmpgas = hotdiffusegas[0:10]

## Abundances in solar units of Anders & Grevesse (1989)
tmpgas['abundance'] = g_p.abundance_to_solar(tmpgas['mass_fraction'], elements=elems)

# Convert mass fractions of each element to their densities relative to H density (ie. abundance X/H)
# tmpgas['nX/nH'] = mass_fraction_to_relative_abundance(tmpgas['mass_fraction'], elements=elems)

## Is it important to calculate L_X with abundance in solar units which I will be using later, or does it not matter?
## Or should L_X be calculated without any solar normalization? --> pyatomdb needs AG89 solar normalized abundances
for xband_, xlims_, xfile_, xmode_, xnbins_ in zip(xray_bands, xray_band_limits, xray_emissivity_files, xray_modes, xray_nbins):
    print(xband_)
    ## Total X-ray luminosity
    tmpgas['Lx_'+ xband_] = g_p.calcu_luminosity(gas=tmpgas, filename=xfile_, mode=xmode_, band=xlims_, elements=elems, bins=xnbins_)
    
    ## Cooling time calculation from Braspenning+2024 eqn 18
    tmpgas['tcool_' + xband_] = 6 * tmpgas['ne'] * kB * tmpgas['T'] / (2 * tmpgas['ElectronAbundance'] * tmpgas['Lx_'+ xband_] / 
                                                                       tmpgas['volume'].in_units('cm**3'))

for i in range(nsteps, len(hotdiffusegas['abundance'][:,0]), nsteps):
    tmpgas = hotdiffusegas[i-nsteps:i]
    tmpgas['abundance'] = g_p.abundance_to_solar(tmpgas['mass_fraction'], elements=elems)
#     tmpgas['nX/nH'] = mass_fraction_to_relative_abundance(tmpgas['mass_fraction'], elements=elems)
    for xband_, xlims_, xfile_, xmode_, xnbins_ in zip(xray_bands, xray_band_limits, xray_emissivity_files, xray_modes, xray_nbins):
        tmpgas['Lx_'+ xband_] = g_p.calcu_luminosity(gas=tmpgas, filename=xfile_, mode=xmode_, band=xlims_, elements=elems, bins=xnbins_)
        tmpgas['tcool_' + xband_] = 6 * tmpgas['ne'] * kB * tmpgas['T'] / (2 * tmpgas['ElectronAbundance'] * tmpgas['Lx_'+ xband_] / 
                                                                           tmpgas['volume'].in_units('cm**3'))
        
    count = i

tmpgas = hotdiffusegas[count:]
tmpgas['abundance'] = g_p.abundance_to_solar(tmpgas['mass_fraction'],elements=elems)
# tmpgas['nX/nH'] = mass_fraction_to_relative_abundance(tmpgas['mass_fraction'], elements=elems)
for xband_, xlims_, xfile_, xmode_, xnbins_ in zip(xray_bands, xray_band_limits, xray_emissivity_files, xray_modes, xray_nbins):
    tmpgas['Lx_'+ xband_] = g_p.calcu_luminosity(gas=tmpgas, filename=xfile_, mode=xmode_, band=xlims_, elements=elems, bins=xnbins_)
    tmpgas['tcool_' + xband_] = 6 * tmpgas['ne'] * kB * tmpgas['T'] / (2 * tmpgas['ElectronAbundance'] * tmpgas['Lx_'+ xband_] / 
                                                                       tmpgas['volume'].in_units('cm**3'))

print('DONE')
print()



## Halo Analysis from saved halo prop files
## Get info from halo prop file
# temp_halo_info = {'M':{}, 'L':{}, 'T':{}, 'R':{}, 'S':{}, 'ne':{}, 'nh':{}, 'N_lgal':{}, 'halo_id':{}, 'metals':{}}
temp_halo_info = {}

print()
print('Reading in saved halo prop file')
with h5py.File(host_prop_file, "r") as f:
    quantity_dict = f.keys()
    for quan in quantity_dict:
        temp_halo_info[quan] = {}
        if quan in ['group_list', 'halo_id', 'N_lgal']:
            temp_halo_info[quan] = f[quan][:]
        else:
            detailed_quantity_dict = f[quan].keys()
            for d_quan in detailed_quantity_dict:
                temp_halo_info[quan][d_quan] = {}
                temp_halo_info[quan][d_quan] = \
                        pnb.array.SimArray(f[quan][d_quan][:], units=f[quan][d_quan].attrs['units'])

print('DONE')
print()

temp_halo_info_nhalos = len(temp_halo_info['halo_id'])
           
print('temp_halo_info keys:')
print(temp_halo_info.keys())
print()
print()
print('temp_halo_info halo_id:')
print(temp_halo_info['halo_id'][:50])
print(temp_halo_info['halo_id'][-51:])
print('end id:', temp_halo_info['halo_id'][-1])
print('max id:', max(temp_halo_info['halo_id']))
print()
print('Number of halos:', temp_halo_info_nhalos)
print('Min M500:', min(temp_halo_info['M']['500']))
print('Max M500:', max(temp_halo_info['M']['500']))
print('Min N_lgal:', min(temp_halo_info['N_lgal']))
print('Max N_lgal:', max(temp_halo_info['N_lgal']))
print()

print()
print('temp_halo_info:')
print(temp_halo_info)
print()
print()
                

print()
print('Applying N_lgal and M500 cuts')
Nlgal_cut = temp_halo_info['N_lgal'] >= N_lgal
Migrm_cut = temp_halo_info['M']['igrm500'] > 0
M500_cut = temp_halo_info['M']['500'] >= 10**(log_M500_cut_val) #10**(13.0)
group_idxs, = np.nonzero(Nlgal_cut & Migrm_cut & M500_cut)
group_ids = temp_halo_info['halo_id'][group_idxs]
print('Number of groups:', len(group_ids))
print()
print('group_ids:')
print(group_ids)
print()
print('group_idxs:')
print(group_idxs)
print()
temp_halo_info['group_idxs'] = group_idxs
temp_halo_info['group_ids'] = group_ids
print('DONE')
print()



## Halo analysis here so that can centre on halos
print()
print('Halo analysis to centre halos:')
print()
print('Getting ahf halos')
# h = s.halos(ahf_mpi=True)
h = s.halos(ahf_mpi=True, use_iord=True, nthreads=nthreads, ahf_basename=ahf_basename)
print('\nDONE')
print()
print()
print('Getting XIGrM halos')
halo = h_a.halo_props(h, datatype='gizmo_ahf')
print('\nDONE')
print()
print()
print('Calculating lower mass limit of galaxies')
galaxy_low_limit_ = 64 * s.star['mass'].mean()
print('\nDONE')
print()
print('Initializing galaxy_low_limit')
halo.init_relationship(galaxy_low_limit=galaxy_low_limit_)
print('\nDONE')
# halo.init_metallicity(elements=elems, radii=['500'], weight_types=['mass', 'Lx']) ## only if want to re-open prop file to make changes
print()
print('Calculating halo centres')
halo.get_center()
print('\nDONE')
print()
print('Halo Error List:')
print(halo.errorlist)
print()

print()
print('h:')
print(h)
print()
print()
print('halo.prop:')
print(halo.prop)
print()


h_nhalos = len(h)
halo_prop_nhalos = np.shape(halo.prop['M'])[0]
print()
print('end temp_halo_info id:', temp_halo_info['halo_id'][-1])
print('h number of halos:', h_nhalos)
print('temp_halo_info number of halos:', temp_halo_info_nhalos)
print('halo.prop number of halos:', halo_prop_nhalos)
print()
print('DONE')
print()



## Calculate extra global properties
print()
print('Calculating extra global properties')
print()

r_units = 'kpc'
M_units = 'Msol'

extra_props_ = ['T', 'ne', 'S', 'P']
for extra_prop_ in extra_props_:
    if extra_prop_ not in temp_halo_info.keys():
        temp_halo_info[extra_prop_] = {}
        

## temp_halo_info
## Test definitions of T500, P500, K500 from Braspenning+2024 (Flamingo paper)?
print('Adding them to temp_halo_info')
print()
M200_ = temp_halo_info['M']['200'].in_units(M_units)
M500_ = temp_halo_info['M']['500'].in_units(M_units)
R200_ = temp_halo_info['R']['200'].in_units(r_units)
R500_ = temp_halo_info['R']['500'].in_units(r_units)
T200_ = (G * M200_ * mu * mp / (kB * R200_)).in_units('K')
T500_ = (G * M500_ * mu * mp / (kB * R500_)).in_units('K')
ne200_ = (200. * rho_gas_scaling / (mu_e * mp)).in_units('cm**-3')
ne500_ = (500. * rho_gas_scaling / (mu_e * mp)).in_units('cm**-3')
ne200_ = pnb.array.SimArray(np.full(len(M200_), ne200_), units=ne200_.units)
ne500_ = pnb.array.SimArray(np.full(len(M500_), ne500_), units=ne500_.units)
K200_ = pnb.array.SimArray(kB * T200_ * ne200_**(-2./3.), units=temp_halo_info['S']['500'].units)
K500_ = pnb.array.SimArray(kB * T500_ * ne500_**(-2./3.), units=temp_halo_info['S']['500'].units)
P200_ = kB * T200_ * ne200_
P500_ = kB * T500_ * ne500_

temp_halo_info['T']['200ideal'] = T200_
temp_halo_info['T']['500ideal'] = T500_
temp_halo_info['ne']['200ideal'] = ne200_
temp_halo_info['ne']['500ideal'] = ne500_
temp_halo_info['S']['200ideal'] = K200_
temp_halo_info['S']['500ideal'] = K500_
temp_halo_info['P']['200ideal'] = P200_
temp_halo_info['P']['500ideal'] = P500_

print('DONE')
print()



## Calculate profiles

print()
print('Calculating profiles:')
print()


xaxis_units_dict = {
    'R500':'',
    'M500':r_units + ' ' + M_units,
    'R200':'',
    'Physical':r_units,
}

# xaxis_units_dict = {
#     'R500':r_units,
#     'log(M500)':''
#     'R500*log(M500)':r_units
#     'R200':r_units,
#     'log(M200)':'',
#     'R200*M200':r_units,
#     'Physical':'',
# }



if xscale.lower() == 'r500':
    # Small Bins
    # _rbin1 = np.append([0], np.logspace(-2, 0, 37))
    # _rbin2 = np.logspace(0, np.log10(2), 7)

    # Wide Bins
    # _rbin1 = np.append([0], np.logspace(-1, 0, 10))
    # _rbin2 = np.logspace(0, np.log10(3), 10)

    # Intermediate Bins
    _rbin1 = np.append([0], np.logspace(-2, 0, 10))
    _rbin2 = np.logspace(0, np.log10(3), 5)
    _rbin3 = np.logspace(np.log10(3), 1, 10)
    _rbins = np.unique(np.concatenate((_rbin1, _rbin2, _rbin3)))

    # Fine core Bins
    # _rbin1 = np.append([0], np.logspace(-3, 0, 30))
    # _rbin2 = np.logspace(0, np.log10(3), 5)

elif xscale.lower() == 'physical':
    # Intermediate Bins for physical radius
#     _rbin1 = np.append([0], np.logspace(0.5, 3, 10)) # 20?  ## Physical softening length is ~2 kpc
    _rbin0 = np.append([0], np.logspace(0.5, 2, 10))
    _rbin1 = np.logspace(2, 3, 10)
    _rbin2 = np.logspace(3, np.log10(3000), 5)
    _rbin3 = np.logspace(np.log10(3000), 4, 10)
#     _rbins = np.unique(np.concatenate((_rbin1, _rbin2, _rbin3)))
    _rbins = np.unique(np.concatenate((_rbin0, _rbin1, _rbin2, _rbin3)))


# _rbins = np.unique(np.append(_rbin1, _rbin2))
rbins = pnb.array.SimArray(_rbins, units=xaxis_units_dict[xscale])
log_rbins = np.log10(rbins)
bin_centres = (rbins[:-1] + rbins[1:])/2.
log_bin_centres = (log_rbins[:-1] + log_rbins[1:])/2.


xunit_converter = {
    'R500':'R500',
    'R200':'R200',
    'Physical':r_units,
}

print('rbins [%s]:' % xunit_converter[xscale])
print(rbins)
print()


## Radial bins to get the core for CC/NCC classification (set to be in units of R500)
_rbins_core = np.array([0, 0.048, 0.05])
rbins_core = pnb.array.SimArray(_rbins_core, units='')
log_rbins_core = np.log10(rbins_core)
bin_centres_core = (rbins_core[:-1] + rbins_core[1:])/2.
log_bin_centres_core = (log_rbins_core[:-1] + log_rbins_core[1:])/2.


## Larger radial bins to get the core for CC/NCC classification (set to be in units of R500)
_rbins_big_core = np.array([0, 0.1, 0.15])
rbins_big_core = pnb.array.SimArray(_rbins_core, units='')
log_rbins_big_core = np.log10(rbins_core)
bin_centres_big_core = (rbins_core[:-1] + rbins_core[1:])/2.
log_bin_centres_big_core = (log_rbins_core[:-1] + log_rbins_core[1:])/2.


## Add in option to create profiles with various different scalings of x axis (R500, R200, Physical, ...)

final_result = {}
final_result['xaxis'] = {
    'all':{
        'scale':xscale,
        'bins':rbins,
        'log_bins':log_rbins,
        'bin_centres':bin_centres,
        'log_bin_centres':log_bin_centres,
    },
    'core':{
        'scale':'R500',
        'bins':rbins_core,
        'log_bins':log_rbins_core,
        'bin_centres':bin_centres_core,
        'log_bin_centres':log_bin_centres_core,
    },
    'core_big':{
        'scale':'R500',
        'bins':rbins_big_core,
        'log_bins':log_rbins_big_core,
        'bin_centres':bin_centres_big_core,
        'log_bin_centres':log_bin_centres_big_core,
    },
}
final_result['all_profiles'] = {}
final_result['all_core_profiles'] = {}
final_result['all_big_core_profiles'] = {}
group_ids_to_remove = np.array([])
group_core_ids = []
group_big_core_ids = []

for i in range(len(group_ids)):
    halo_id = group_ids[i]
    halo_idx = group_idxs[i]
    
    print('Group', i+1, '/', len(group_ids))
    print('Halo', halo_id, 'with idx', halo_idx)
    
    xscale_values = {
        'R500':pnb.array.SimArray(temp_halo_info['R']['500'].in_units(r_units)[halo_idx], units=r_units),
        'R200':pnb.array.SimArray(temp_halo_info['R']['200'].in_units(r_units)[halo_idx], units=r_units),
        'Physical':pnb.array.SimArray(1, units=''),
    }
    
    xscale_value = xscale_values[xscale]
    print(xscale, ':', xscale_value, xscale_value.units)
    
    
    ## For binning
    temp_rbins = rbins * xscale_value
    temp_radius = pnb.array.SimArray(temp_rbins[-1], units=r_units)

    ## For determining CC/NCC classification from tcool of core
    temp_rbins_core = rbins_core * xscale_values['R500']
    temp_radius_core = pnb.array.SimArray(temp_rbins_core[-1], units=r_units)
    
    temp_rbins_big_core = rbins_big_core * xscale_values['R500']
    temp_radius_big_core = pnb.array.SimArray(temp_rbins_big_core[-1], units=r_units)
    
    
    center = halo.center[halo_id - 1]
    print('Center:', center)
    
    ## Extra factor of 1.1 just to make sure all particles are included
    if halo_filt.lower() == 'sphere':
        filter_all = pnb.filt.Sphere(radius=1.1*temp_radius, cen=center)
        filter_core = pnb.filt.Sphere(radius=1.1*temp_radius_core, cen=center)
        filter_big_core = pnb.filt.Sphere(radius=1.1*temp_radius_big_core, cen=center)
    elif halo_filt.lower() == 'cylinder':
        ## By default, oriented along the z-axis
        filter_all = pnb.filt.Disc(radius=1.1*temp_radius, height=3*xscale_values['R500'], cen=center)
        filter_core = pnb.filt.Disc(radius=1.1*temp_radius_core, height=3*xscale_values['R500'], cen=center)
        filter_big_core = pnb.filt.Disc(radius=1.1*temp_radius_big_core, height=3*xscale_values['R500'], cen=center)
    
   
    temp_gas = hotdiffusegas[filter_all]
    print('Length of temp_gas:', len(temp_gas))
    
    temp_gas_core = hotdiffusegas[filter_core]
    print('Length of temp_gas_core:', len(temp_gas_core))
    
    temp_gas_big_core = hotdiffusegas[filter_big_core]
    print('Length of temp_gas_big_core:', len(temp_gas_big_core))
    
    # Need halo to not be empty of gas
    if xscale_value>0 and len(temp_gas)>0:
        temp_result = gas_profile(temp_gas, center, temp_rbins, xray_bands, 
                                  profile_type=profile_type, profile_ndims=profile_ndims, weight_by=weight_by)
        final_result['all_profiles'][halo_id] = temp_result
        if len(temp_gas_core)>0:
            temp_result_core = gas_profile(temp_gas_core, center, temp_rbins_core, xray_bands, 
                                           profile_type=profile_type, profile_ndims=profile_ndims, weight_by=weight_by)
            final_result['all_core_profiles'][halo_id] = temp_result_core
            group_core_ids.append(halo_id)
        if len(temp_gas_big_core)>0:
            temp_result_big_core = gas_profile(temp_gas_big_core, center, temp_rbins_big_core, xray_bands, 
                                           profile_type=profile_type, profile_ndims=profile_ndims, weight_by=weight_by)
            final_result['all_big_core_profiles'][halo_id] = temp_result_big_core
            group_big_core_ids.append(halo_id)
    else:
        print("removing idx", halo_idx)
        group_ids_to_remove = np.append(group_ids_to_remove, halo_id)
        
    print()
print()


# Remove the indices of halos that were empty of gas
group_ids = np.array([group_id for group_id in group_ids if group_id not in group_ids_to_remove])
group_core_ids = np.array(group_core_ids)
group_big_core_ids = np.array(group_core_ids)

final_result['all_id'] = group_ids
final_result['all_core_id'] = group_core_ids
final_result['all_big_core_id'] = group_core_ids
final_result['temp_halo_info'] = temp_halo_info
final_result['cosmo_props'] = {
    'age':age,
    'z':z,
    'h0':h0,
    'OmegaL0':omL0,
    'OmegaM0':omM0,
    'OmegaB0':omB0,
    'rho_crit':rho_crit,
    'rho_gas_scaling':rho_gas_scaling,
}
print(final_result)
print('DONE')
print()


# Save profiles
print()
print('Saving profiles to file')
save_object_with_dill(final_result, save_file)
print(save_file)
print('\nDONE\n')