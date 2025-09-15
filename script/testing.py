import numpy as np
import pynbody as pnb
import os
import h5py
from pynbody.analysis.profile import Profile, VerticalProfile, InclinedProfile, QuantileProfile
import argparse
import dill
import caesar
import yt
import astropy.constants as astroc

m_p = astroc.m_p.cgs.value # mass of proton in cgs units
k_B = astroc.k_B.cgs.value # Boltzmann constant in cgs units
from astro_constants import NA_no_units, NA, kB_no_units, kB, mH, mp_no_units, mu_e, mu, G_no_units, G

profile_types = {
    'Profile':Profile,
    'VerticalProfile':VerticalProfile,
    'InclinedProfile':InclinedProfile,
    'QuantileProfile':QuantileProfile,
}

code = 'Simba-C'
name = 'entropy_prof_L1_zoom_252_644' # for output file name
#snap_file = args.snap
#output_file = args.output_filepath + '_temp_cut='+args.temp_cut+'_nh_cut='+args.nh_cut+'.pkl'
output_filename = name + '_temp_cut=5e5K_nh_cut=0.13_cm-3.pkl'
profile_type = profile_types['Profile']
weight_by = 'mass'

# define args
snap_file = '/scratch/b/babul/fjenning/xray_bubbles/PAPER2_Level1_zoom_252/snaps/FT2Myr_snap_halo_252_644.hdf5'
output_path = '/home/b/babul/stlock/entropy_profiles/outputs/'
temp_cut = '5e5 K' # default for now
nh_cut = '0.13 cm**-3' # default for now

print()
print('Reading in simulation snapshot')
s = pnb.load(snap_file)
s.physical_units()
print('done')
print()

#useful functions, XIGrM, dill, etc

def save_object_with_dill(obj, filename):
    with open(filename, 'wb') as f:  # Overwrites any existing file.
        dill.dump(obj, f, dill.HIGHEST_PROTOCOL)

# XIGrM Functions
def nh(sim):
    '''
    Calculating hydrogen number density from density.
    '''
    #X_H = 1-sim.gas['metals'][:,0]-sim.gas['metals'][:,1] #0 for metal and 1 for He
    result = sim.gas['rho'].in_units('g cm**-3').view(np.ndarray)/m_p*sim.gas['X_H']
    result = pnb.array.SimArray(result)
    result.units = 'cm**-3'
    return result


def temp(sim):
    '''
    Convert internal energy to temperature, following 
    the instructions from `GIZMO documentation <http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html#snaps>`_
    '''
    #XH = 1-sim.gas['metals'][:,0]-sim.gas['metals'][:,1]
    GAMMA = 5.0/3
    # If ne = 0 gas is unionized
    mu = 4.0 / (3.0 * sim.gas['X_H'] + 4.0 * sim.gas['X_H'] * sim.gas['ElectronAbundance'] + 1.0)
    # mu = 0.6
    result = pnb.array.SimArray(sim.gas['u'].in_units('cm**2 s**-2').view(np.ndarray) * (GAMMA - 1) *\
                                 mu * (m_p / k_B))
    result.units = 'K'
    return result

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
def volume(self):
    return self._binsize

@Profile.profile_property # call as something like result['y'][prop + '_' + weight] = prop_weighted_avg(pg, prop, weight)
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



massive = halo_data[0].virial_quantities['m500c']
target_halo = halo_data[0]
for halo in halo_data: # find most massive halo
    if halo.virial_quantities['m500c'] > massive:
        target_halo = halo

# get the center of that halo!
center = target_halo.minpotpos
r500 = target_halo.virial_quantities['r500c']

# convert to kpc
a = s.properties['a']
center_kpc = pnb.array.SimArray(center / a, units='kpc')  # Force conversion
r500_kpc = pnb.array.SimArray(r500 / a, units='kpc')  # Force conversion


# caesar halo finding
obj = caesar.load('/scratch/b/babul/fjenning/xray_bubbles/PAPER2_Level1_zoom_252/Groups/Caesar_FT2Myr_halo_252_644.hdf5')
halo_data = obj.halos

gas_positions = s.gas['pos'] # get positions of all gas particles in box

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

#s.gas['u'] /= pnb.array.SimArray(s.properties['a'], units='1')

# Hydrogen mass fraction
if code.lower() == 'simba': # might as well keep the simba stuff
    s.gas['X_H'] = 1-s.gas['metals'][:,metals_idx['Z']]-s.gas['metals'][:,metals_idx['He']]
elif code.lower() == 'simba-c':
    s.gas['X_H'] = s.gas['metals'][:,metals_idx['H']]

s.gas['nh'] = nh(s) # Hydrogen number density
s.gas['ElectronAbundance'].units = '1'
s.gas['T'] = temp(s) # From internal energy to temperature
s.gas['temp'] = s.gas['T']
print(f"Total gas particles: {len(s.gas)}")


## Get hot diffuse gas using provided cuts
print('Isolating IGrM (temperate, hydrogen number density, and decoupled wind particle cuts)')
temp_cut = '5e5 K' # default for now
nh_cut = '0.13 cm**-3' # default for now
temp_filt = pnb.filt.HighPass('T', temp_cut)
density_filt = pnb.filt.LowPass('nh', nh_cut)
wind_filt = ~pnb.filt.HighPass('DelayTime', 0)  # Cut out hot decoupled wind particles
total_filt = temp_filt & density_filt & wind_filt


hotdiffusegas = s.gas[total_filt]
print(f"Particles after filtering: {len(hotdiffusegas)}")

# calc entropy
hotdiffusegas['ne'] = hotdiffusegas['ElectronAbundance'] * hotdiffusegas['nh'].in_units('cm**-3') # eletron number density
hotdiffusegas['K'] = kB * hotdiffusegas['T'] * hotdiffusegas['ne']**(-2,3) # entropy  
hotdiffusegas['pos'].convert_units('kpc')


def gas_profile(snapshot, center, rbins, profile_type=Profile, profile_ndims=3, weight_by='mass'):
    tx = pnb.transformation.inverse_translate(snapshot, center)
    boxsize = snapshot.properties['boxsize']
    
    log_rbins = np.log10(rbins)
    bin_centres = (rbins[:-1] + rbins[1:])/2.
    log_bin_centres = (log_rbins[:-1] + log_rbins[1:])/2.

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
        
        
        
        ## x-axis values (normalized radius)
        
        result['x']['rbins'] = rbins
        result['x']['log_rbins'] = log_rbins
        result['x']['bin_centres'] = bin_centres
        result['x']['log_bin_centres'] = log_bin_centres
        
        
        
        ## y-axis values (entropy)
        
        result['y']['n'] = pnb.array.SimArray(pg['n'], units='') # KEEP, this is the number of particles, we need this
        result['y']['volume'] = pg['volume']
        result['y']['Ne'] = pg['Ne']
        result['y']['mass'] = pg['mass']
        
        result['y']['T_mass'] = prop_weighted_avg(pg, 'temp', 'mass')
        result['y']['K_mass'] = prop_weighted_avg(pg, 'K', 'mass')

                    
        snapshot['pos'] = original_pos
    
    return result


# radial binning, 40 bins

rbins = np.append([0], np.logspace(-2, np.log10(10), 40))
xaxis_units_dict = {'R500':''}
r_units = 'kpc'

rbins = pnb.array.SimArray(rbins, units=xaxis_units_dict['R500'])
r500 = pnb.array.SimArray(r500, units='kpc') 
r500 = r500.in_units(r_units) 

# radial binning, 40 bins
rbins = np.append([0], np.logspace(-2, np.log10(10), 40))
xaxis_units_dict = {'R500':''}
r_units = 'kpc'

rbins = pnb.array.SimArray(rbins, units=xaxis_units_dict['R500'])
r500 = pnb.array.SimArray(r500, units='kpc') 
r500 = r500.in_units(r_units) 
# r500 = pnb.array.SimArray(r500.to(r_units), units=r_units)

## For binning
temp_rbins = rbins * r500
temp_radius = pnb.array.SimArray(temp_rbins[-1], units=r_units)

filter_all = pnb.filt.Sphere(radius=1.1*r500_kpc, cen=center_kpc)
temp_gas = hotdiffusegas[filter_all] #throw this into the gas_profile function then done!

print('temp gas size:', len(temp_gas))
print('temp gas:')
print(temp_gas)
print('center:')
print(center)
print('rbins:')
print(rbins[0])

print(f"temp_gas length: {len(temp_gas)}")
print(f"temp_gas contents: {temp_gas}")

final_prod = gas_profile(temp_gas, center, rbins, profile_type=Profile, profile_ndims=3, weight_by='mass')

# Saving the file

output =  output_path + output_filename
# Save profiles
print()
print('Saving profiles to file')
save_object_with_dill(final_prod, output)
print("saved as "+ output_filename + " in " + output_path)
print('\nDONE\n')