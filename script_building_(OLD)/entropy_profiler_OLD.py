import yt
import unyt
from unyt import unyt_quantity, unyt_array, km, s, Mpc, kg, m, Msun, kpc
from unyt import Gyr, Mpc, km, s
from unyt.physical_constants import G as G_unyt, kboltz, mp as mp_unyt  # G, k_B, m_p

import numpy as np
import matplotlib.pyplot as plt

import caesar
import pynbody as pnb
import pynbody.units as u
import pynbody.array as pna
import argparse
import sys

import dill
import os
from astro_constants import NA_no_units, NA, kB_no_units, kB, mH, mp_no_units, mp, mu_e, mu, G_no_units, G
#from astro_constants import NA_no_units, NA, kB_no_units, kB, mH, mp_no_units, mp, mu_e, mu, G_no_units,


parser = argparse.ArgumentParser(description="Create Entropy profiles for a given halo")
parser.add_argument('--snap', action='store', type=str, required=True, 
                   help='Snapshot filepath')
parser.add_argument('--caesarfile', action='store', type=str, required=False, 
                   help='Caesar file filepath')
parser.add_argument('--halo', action='store', type=int, required=True, 
                   help='Desired Halo ID')
parser.add_argument('--code', action='store', type=str, required=True, 
                   help='Either Simba or Simba-C')
parser.add_argument('--idpath', action='store', type=str, required=False, 
                   help='Path to halo id\'s')
args = parser.parse_args()


snap_file = args.snap
caesar_file = args.caesarfile
target_halo = args.halo
code = args.code
id_filepath = args.idpath


# code = 'Simba'
# target_halo = 3224
# snap_file = '/scratch/b/babul/aspadawe/snapshots/Hyenasc_L1_zooms/halo_3224/snap_halo_3224_151.hdf5'
# id_filepath = "/home/b/babul/stlock/entropy_profiles/finished/ids_formatted.txt" # containing halo data/ids
# caesar_file = "/scratch/b/babul/aspadawe/snapshots/Hyenasc_L1_zooms/halo_3224/Groups/Caesar_halo_3224_151.hdf5"

snap = yt.load(snap_file)#, unit_system='cgs')
s = pnb.load(snap_file)


#------------------------------------------------------------------------------------ METAL INDICES & FUNCTIONS

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

def save_object_with_dill(obj, filename):
    with open(filename, 'wb') as f:  # Overwrites any existing file.
        dill.dump(obj, f, dill.HIGHEST_PROTOCOL)
    
 
data = snap.all_data()

def gas_entropy(field, data):
    return data["gas", "kT"] * data["gas", "El_number_density"] ** (-2.0/3.0)

def gas_entropy2(field, data):
    return data["gas", "kT"] * data["PartType0", "El_number_density"] ** (-2.0/3.0)

snap.add_field(
    ("gas", "gas_entropy"),
    units="keV*cm**2",
    function=gas_entropy,
    sampling_type="cell",
)

snap.add_field(
    ("PartType0", "gas_entropy"),
    units="keV*cm**2",
    function=gas_entropy2,
    sampling_type="particle",
    force_override=True
)

gas_entropy = data['gas', 'entropy']

def get_coordinates_by_id(id_filepath, target_id):
    with open(id_filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 6:
                continue  # skip malformed lines
            
            try:
                obj_id = int(parts[2])
                if obj_id == target_id:
                    # parts[4] contains the coordinate string: (x,y,z)
                    coord_str = parts[4].strip('()')
                    x, y, z = map(float, coord_str.split(','))
                    return (x, y, z)
            except (ValueError, IndexError):
                continue  # skip malformed lines
    
    return None  # if ID not found

# ------------------------------------------------------------------------------------ POSITION & IGrM FILTERING

coords = get_coordinates_by_id(id_filepath, target_halo)
halo_center = unyt_array(coords, 'Mpc')
halo_center/=1000
#nprint(halo_center)

igrm = data.include_above(('PartType0', 'temperature'), 5e5, units='K')
igrm = igrm.include_below(('PartType0', 'H_nuclei_density'), 0.13, units='cm**-3')
igrm = igrm.include_equal(('PartType0', 'DelayTime'), 0)

# Radial Binning

radial_bins = unyt.unyt_array(np.append([0], np.logspace(-1, np.log10(5), 20)), 'Mpc')

igrm_pos_x = igrm['gas', 'particle_position_x']
igrm_pos_y = igrm['gas', 'particle_position_y']
igrm_pos_z = igrm['gas', 'particle_position_z']

igrm_halo_center = ((igrm_pos_x-halo_center[0])**2 + (igrm_pos_y-halo_center[1])**2 + 
                  (igrm_pos_z-halo_center[2])**2)**0.5



# ------------------------------------------------------------------------------------ pos of highest gas density

rho_max_index = np.argmax(data['gas', 'density'])

rho_max_pos = data['gas', 'position'][rho_max_index]
rho_max_pos.in_units('Mpc')

rho_max_part_pos = data['gas', 'particle_position'][rho_max_index]
rho_max_part_pos.in_units('Mpc')

igrm_zoom_radii_rel_rho_max = ((igrm_pos_x-rho_max_pos[0])**2 + 
                  (igrm_pos_y-rho_max_pos[1])**2 + 
                  (igrm_pos_z-rho_max_pos[2])**2)**0.5
# ------------------------------------------------------------------------------------ ENTROPY CALCULATING
igrm_zoom_entropies = igrm['gas', 'entropy']
igrm_zoom_masses = igrm['gas', 'mass']

# change igrm_zoom_radii = .... to use the center you want!
igrm_zoom_radii = igrm_zoom_radii_rel_rho_max
igrm_zoom_radii.in_units('Mpc')

entropy_profile_median = []
entropy_profile_mean = []
entropy_profile_mass_weighted = []
for bin_lo, bin_hi in zip(radial_bins[:-1], radial_bins[1:]):
    radial_filter = np.logical_and(igrm_zoom_radii.in_units('Mpc').value>=bin_lo.in_units('Mpc').value,
                                   igrm_zoom_radii.in_units('Mpc').value<=bin_hi.in_units('Mpc').value)
    entropies = igrm_zoom_entropies[radial_filter]
    masses = igrm_zoom_masses[radial_filter]
    
    entropy_profile_median.append(np.nanmedian(entropies))
    entropy_profile_mean.append(np.nanmean(entropies))
    if sum(masses)>0:
        K_mass = np.average(entropies, weights=masses)
        entropy_profile_mass_weighted.append(K_mass)
    else:
        entropy_profile_mass_weighted.append(np.nan)
    
# ------------------------------------------------------------------------------------ K500
"""
def compute_K500(caesar_object, halo_index):

    halo = caesar_object.halos[halo_index]

    # cosmological parameters from Caesar
    M500 = halo.virial_quantities["m500c"]  # Msun
    R500 = halo.virial_quantities["r500c"]  # kpc

    omega_b = caesar_object.simulation.omega_baryon
    omega_m = caesar_object.simulation.omega_matter
    omega_l = caesar_object.simulation.omega_lambda
    z = caesar_object.simulation.redshift
    h = caesar_object.simulation.hubble_constant

    # Compute the Hubble parameter at some z
    # H0 ~ 100*h km/s/Mpc  → convert to 1/s
    H0_si = (100.0 * h * (km/s) / Mpc).to_value("1/s")  # float in s^-1
    E_z = np.sqrt(omega_m*(1+z)**3 + omega_l)
    Hz_si = H0_si * E_z  # s^-1

    # Critical density (in Msun/kpc^3)
    rho_crit_val = 3.0 * Hz_si**2 / (8.0 * np.pi * G.value)  # float in kg/m^3
    rho_crit_si  = unyt_quantity(rho_crit_val, "kg/m**3")
    rho_crit_z   = rho_crit_si.to("Msun/kpc**3")

    # Electron density at R500:
    f_b         = omega_b / omega_m
    rho_gas_500 = 500.0 * f_b * rho_crit_z  # Msun/kpc^3
    n_e_500     = (rho_gas_500 / (mu_e * mp)).to("cm**-3")

    # T_{500} = (G * M500 * mu * m_p) / (k_B * R500)
    M500_si = M500.to("kg")
    R500_si = R500.to("m")
    T500_si = G * M500_si * mu * mp / (kboltz * R500_si)  # unyt in Kelvin

    # K500 = k_B * T500 * (n_e,500)^{-2/3}, convert to keV cm^2
    K500_si       = kboltz * T500_si * (n_e_500**(-2.0/3.0))
    K500_keV_cm2  = K500_si.to("keV*cm**2")

    return K500_keV_cm2

K500 = compute_K500(obj, halo_id)
print("K500 = ", K500)
"""

# ------------------------------------------------------------------------------------ PLOTTING

radial_bin_centres = (radial_bins[:-1] + radial_bins[1:])/2.
radial_bin_centres_logged = np.log10(radial_bin_centres)
log_radial_bin_centres = (np.log10(radial_bins[:-1]) + np.log10(radial_bins[1:]))/2.
log_radial_bin_centres_unlogged = 10**log_radial_bin_centres

# normalized_r = radial_bin_centres/R500
# normalized_k = np.array(entropy_profile_mean)/K500.value    

plt.title("Halo " + str(target_halo) + " Entropy Profile")
plt.plot(radial_bin_centres, entropy_profile_median, label='median', lw=1)
plt.plot(radial_bin_centres, entropy_profile_mean, label='mass-weighted/mean', lw=1)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('R [Mpc]')
plt.ylabel(r'K [keV cm$^2$]')
plt.legend()
# plt.show()

#plt.savefig("halo_" + str(target_halo)+ "_entropy_plot.png")
print("plot successfully made")

"""
fig = plt.gcf()  # Get current figure

# Save profiles
save_object_with_dill(fig, "halo_" + str(target_halo)+ "_entropy_plot.dill")
print("profile successfully saved")

# load plot with the following:

with open("halo_..._entropy_plot.dill", "rb") as f:
    fig = dill.load(f)

fig.show()
"""