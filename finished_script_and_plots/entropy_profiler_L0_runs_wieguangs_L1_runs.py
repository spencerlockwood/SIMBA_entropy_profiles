# Simulation data

# Observational data
low_mass_range = 12.8
high_mass_range = 13.2

# caesar_group_id = 4

# Plotting/saving
#savefilename = 'Snap' + str(zoom_ID) + '_Halo_contamination.png' # .png name to save plot under
savefilename = '/home/b/babul/stlock/entropy_profiles/outputs/testingtestingtesting_full40_v2.png'
thetitle = "L1 Hyenas Entropy Profiles - Varied Mass Range" # filepath/title of graph


# LOWEST -> HIGHEST BlackHoleFeedbackFactor

# (with 5 and 6) order = [20, 19, 18, 17, 26, 16, 25, 15, 24, 14, 23, 13, 22, 12, 21, 1, 27, 2, 28, 3, 29, 4, 30, 5, 31, 6, 32, 7, 8, 9, 10, 11]

# -------------------------------------------------------------------------------------------------------------------------------------- SIM_DATA
import yt
import unyt
from unyt import unyt_quantity, unyt_array, km, s, Mpc, kg, m, Msun, kpc
from unyt import Gyr, Mpc, km, s
from unyt.physical_constants import G as G_unyt, kboltz, mp as mp_unyt  # G, k_B, m_p

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import caesar
import pynbody as pnb
import pynbody.units as u
import pynbody.array as pna
import seaborn as sns


import dill
import os
from astro_constants import NA_no_units, NA, kB_no_units, kB, mH, mp_no_units, mp, mu_e, mu, G_no_units, G
#from astro_constants import NA_no_units, NA, kB_no_units, kB, mH, mp_no_units, mp, mu_e, mu, G_no_units,

code = 'simba-c'

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

def compute_gas_entropy(field, data):
    return data["gas", "kT"] * data["gas", "El_number_density"] ** (-2.0/3.0)

def compute_gas_entropy2(field, data):
    return data["gas", "kT"] * data["PartType0", "El_number_density"] ** (-2.0/3.0)


# Path to the text file
halo_info_file = "/home/b/babul/stlock/entropy_profiles/finished_script_and_plots/Matched_Elite_halos-at-z0-Level1.txt"

# Skip the first line since it's a comment/header
halo_file_data = np.loadtxt(halo_info_file)


low_mass_range = 12.4
high_mass_range = 14.0

#cmap = cm.plasma  # Choose colormap
cmap = cm.get_cmap('turbo')
norm = mcolors.Normalize(vmin=low_mass_range, vmax=high_mass_range)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

# Start figure
fig, ax = plt.subplots(figsize=(12, 8))



for row in halo_file_data:
    halo_id = int(row[0])               # Origin halo ID (your folder name's number)
    caesar_group_id = int(row[1])       # Caesar group ID
    central_galaxy_id = int(row[2])     # Caesar central galaxy ID
    total_mass = row[3]                  # Total mass of the group
    contam_mass = row[4]                 # Contamination mass

    pos_x, pos_y, pos_z = row[5:8]

    halo_folder = f"halo_{halo_id}"
    snap_file = f'/home/b/babul/wcui/runs/HYENAS/{halo_folder}/snap_halo_{halo_id}_151.hdf5'
    caesarfilepath = f'/home/b/babul/wcui/runs/HYENAS/{halo_folder}/Groups/Caesar_halo_{halo_id}_151.hdf5'

    snap = yt.load(snap_file)#, unit_system='cgs')
    data = snap.all_data()
    # s = pnb.load(snap_file)
    obj = caesar.load(caesarfilepath)

    # convert .txt file position units (kpccm/h) to Mpc
    h = obj.simulation.hubble_constant
    pos_x = pos_x/(1000 * h)
    pos_y = pos_y/(1000 * h)
    pos_z = pos_z/(1000 * h)

    dist_to_ic_center = row[8]

    if ("gas", "gas_entropy") not in snap.field_info:
        snap.add_field(
            ("gas", "gas_entropy"),
            units="keV*cm**2",
            function=compute_gas_entropy,
            sampling_type="cell",
        )
    
    if ("PartType0", "gas_entropy") not in snap.field_info:
            snap.add_field(
            ("PartType0", "gas_entropy"),
            units="keV*cm**2",
            function=compute_gas_entropy2,
            sampling_type="particle",
            force_override=True
        )
            
    gas_entropy = data['gas', 'entropy']

    halo_data = obj.halos

    # halo_center = [pos_x, pos_y, pos_z]
    # print(halo_center) # in units of Mpc

    igrm = data.include_above(('PartType0', 'temperature'), 5e5, units='K')
    igrm = igrm.include_below(('PartType0', 'H_nuclei_density'), 0.13, units='cm**-3')
    igrm = igrm.include_equal(('PartType0', 'DelayTime'), 0)


    # Radial Binning

    # radial_bins = unyt.unyt_array(np.append([0], np.logspace(-3, np.log10(5), 30)), 'Mpc')
    radial_bins = unyt.unyt_array(np.append([0], np.logspace(-4, np.log10(5), 45)), 'Mpc')

    igrm_pos_x = igrm['gas', 'particle_position_x']
    igrm_pos_y = igrm['gas', 'particle_position_y']
    igrm_pos_z = igrm['gas', 'particle_position_z']

    halo_center = halo_data[caesar_group_id].minpotpos
    halo_center = halo_center.in_units('Mpc')
    print(halo_center)

    # isolate igrm

    igrm_halo_center = ((igrm_pos_x-halo_center[0])**2 + (igrm_pos_y-halo_center[1])**2 + (igrm_pos_z-halo_center[2])**2)**0.5

    igrm_zoom_entropies = igrm['gas', 'entropy']
    igrm_zoom_masses = igrm['gas', 'mass']


    entropy_profile_median = []
    entropy_profile_mean = []
    entropy_profile_mass_weighted = []
    for bin_lo, bin_hi in zip(radial_bins[:-1], radial_bins[1:]):
    #     print(bin_lo, bin_hi)
        radial_filter = np.logical_and(igrm_halo_center.in_units('Mpc').value>=bin_lo.in_units('Mpc').value,
                                    igrm_halo_center.in_units('Mpc').value<=bin_hi.in_units('Mpc').value)
    #     print(np.where(radial_filter))
        entropies = igrm_zoom_entropies[radial_filter]
        masses = igrm_zoom_masses[radial_filter]
        # print(len(entropies))
        # print(len(masses))
    #     print()
        
        entropy_profile_median.append(np.nanmedian(entropies))
        entropy_profile_mean.append(np.nanmean(entropies))
        if sum(masses)>0:
            K_mass = np.average(entropies, weights=masses)
            entropy_profile_mass_weighted.append(K_mass)
        else:
            entropy_profile_mass_weighted.append(np.nan)
        
        # print(K_mass)
        # print()



    from unyt import Msun, kpc, cm, keV, unyt_quantity
    from unyt.physical_constants import mp, kboltz, G  # ✅ use unyt's version


    # Constants
    mu = 0.59        # Mean molecular weight
    mu_e = 1.14      # Mean molecular weight per free electron

    # Load CAESAR catalog
    omega_b = obj.simulation.omega_baryon
    omega_m = obj.simulation.omega_matter

    # Update f_b if desired
    f_b = omega_b / omega_m

    # Critical density
    H0_si = (100 * h * 1e3) / (3.086e22)  # H0 in 1/s
    rho_crit_si = 3 * H0_si**2 / (8 * np.pi * 6.6743e-11)  # kg/m^3
    rho_crit = unyt_quantity(rho_crit_si, "kg/m**3").to("Msun/kpc**3")

    # Choose your halo
    halo = obj.halos[caesar_group_id]  # Replace with desired index or search method

    M500 = unyt_quantity(halo.virial_quantities['m500c'].in_units("Msun"), "Msun")
    M200 = unyt_quantity(halo.virial_quantities['m200c'].in_units("Msun"), "Msun")

    R500 = unyt_quantity(halo.virial_quantities['r500c'].in_units("m"), "m") # IN METERS FOR TEMP CALC
    R200 = unyt_quantity(halo.virial_quantities['r200c'].in_units("m"), "m")

    ne_500 = (500 * f_b * rho_crit / (mu_e * mp)).to("cm**-3")
    ne_200 = (200 * f_b * rho_crit / (mu_e * mp)).to("cm**-3")

    T500 = (G * M500 * mu * mp / R500).to("keV") # removing 1/kb here
    T200 = (G * M200 * mu * mp / R200).to("keV") # removing 1/kb here

    K500 = (T500 * ne_500**(-2/3)).to("keV*cm**2") # removing kb here
    K200 = (T200 * ne_200**(-2/3)).to("keV*cm**2") # removing kb here

    R500 = halo_data[caesar_group_id].virial_quantities['r500c'].in_units('Mpc') # redefine r500 and r200 with correct units before binning
    R200 = halo_data[caesar_group_id].virial_quantities['r200c'].in_units('Mpc')
    radial_bin_centres = (radial_bins[:-1] + radial_bins[1:])/2.
    radial_bin_centres_logged = np.log10(radial_bin_centres)
    log_radial_bin_centres = (np.log10(radial_bins[:-1]) + np.log10(radial_bins[1:]))/2.
    log_radial_bin_centres_unlogged = 10**log_radial_bin_centres

    normalized_r = radial_bin_centres/R500.value
    normalized_k_mean = np.array(entropy_profile_mean)/K500.value
    normalized_k_median = np.array(entropy_profile_median)/K500.value


    # axes_.plot(CLoGS_lR,CLoGS_ly,color=CLoGS_fractional_color,lw=6,ls="--", zorder=obs_plot_kwargs['Oppenheimer+21']['zorder']+5,label="CLoGS Groups 1401")

    log_r = np.log10(normalized_r)
    log_k = np.log10(normalized_k_median)

    log_m500c = np.log10(M500)
    color = cmap(norm(log_m500c))

    ax.plot(
        log_r,
        log_k,
        color=color,
        # color = cmap_sim(),
        lw=1,
        ls='--',
        zorder=60,
        alpha=0.7
        #label= f'Halo of M500: ' + str(round(np.log10(halo_data[caesar_group_id].virial_quantities['m500c'].value), 1)) + ', Cont: ' + str(round(halo_data[caesar_group_id].contamination, 3))
        #label=f'simulation {run+1}'

    )


snap_file = '/scratch/b/babul/aspadawe/snapshots/HyenasC/L1/halo_3224_v3/snapshot_151.hdf5'
caesarfilepath = '/scratch/b/babul/stlock/snapshots/caesar_avivs_L1_halo3224.hdf5'

snap = yt.load(snap_file)#, unit_system='cgs')
data = snap.all_data()
# s = pnb.load(snap_file)
obj = caesar.load(caesarfilepath)

# convert .txt file position units (kpccm/h) to Mpc
h = obj.simulation.hubble_constant
pos_x = pos_x/(1000 * h)
pos_y = pos_y/(1000 * h)
pos_z = pos_z/(1000 * h)

dist_to_ic_center = row[8]

if ("gas", "gas_entropy") not in snap.field_info:
    snap.add_field(
        ("gas", "gas_entropy"),
        units="keV*cm**2",
        function=compute_gas_entropy,
        sampling_type="cell",
    )

if ("PartType0", "gas_entropy") not in snap.field_info:
        snap.add_field(
        ("PartType0", "gas_entropy"),
        units="keV*cm**2",
        function=compute_gas_entropy2,
        sampling_type="particle",
        force_override=True
    )
        
gas_entropy = data['gas', 'entropy']

halo_data = obj.halos

# halo_center = [pos_x, pos_y, pos_z]
# print(halo_center) # in units of Mpc

igrm = data.include_above(('PartType0', 'temperature'), 5e5, units='K')
igrm = igrm.include_below(('PartType0', 'H_nuclei_density'), 0.13, units='cm**-3')
igrm = igrm.include_equal(('PartType0', 'DelayTime'), 0)


# Radial Binning

# radial_bins = unyt.unyt_array(np.append([0], np.logspace(-3, np.log10(5), 30)), 'Mpc')
radial_bins = unyt.unyt_array(np.append([0], np.logspace(-4, np.log10(5), 45)), 'Mpc')

igrm_pos_x = igrm['gas', 'particle_position_x']
igrm_pos_y = igrm['gas', 'particle_position_y']
igrm_pos_z = igrm['gas', 'particle_position_z']

halo_center = halo_data[caesar_group_id].minpotpos
halo_center = halo_center.in_units('Mpc')
print(halo_center)

# isolate igrm

igrm_halo_center = ((igrm_pos_x-halo_center[0])**2 + (igrm_pos_y-halo_center[1])**2 + (igrm_pos_z-halo_center[2])**2)**0.5

igrm_zoom_entropies = igrm['gas', 'entropy']
igrm_zoom_masses = igrm['gas', 'mass']


entropy_profile_median = []
entropy_profile_mean = []
entropy_profile_mass_weighted = []
for bin_lo, bin_hi in zip(radial_bins[:-1], radial_bins[1:]):
#     print(bin_lo, bin_hi)
    radial_filter = np.logical_and(igrm_halo_center.in_units('Mpc').value>=bin_lo.in_units('Mpc').value,
                                igrm_halo_center.in_units('Mpc').value<=bin_hi.in_units('Mpc').value)
#     print(np.where(radial_filter))
    entropies = igrm_zoom_entropies[radial_filter]
    masses = igrm_zoom_masses[radial_filter]
    # print(len(entropies))
    # print(len(masses))
#     print()
    
    entropy_profile_median.append(np.nanmedian(entropies))
    entropy_profile_mean.append(np.nanmean(entropies))
    if sum(masses)>0:
        K_mass = np.average(entropies, weights=masses)
        entropy_profile_mass_weighted.append(K_mass)
    else:
        entropy_profile_mass_weighted.append(np.nan)
    
    # print(K_mass)
    # print()



from unyt import Msun, kpc, cm, keV, unyt_quantity
from unyt.physical_constants import mp, kboltz, G  # ✅ use unyt's version


# Constants
mu = 0.59        # Mean molecular weight
mu_e = 1.14      # Mean molecular weight per free electron

# Load CAESAR catalog
omega_b = obj.simulation.omega_baryon
omega_m = obj.simulation.omega_matter

# Update f_b if desired
f_b = omega_b / omega_m

# Critical density
H0_si = (100 * h * 1e3) / (3.086e22)  # H0 in 1/s
rho_crit_si = 3 * H0_si**2 / (8 * np.pi * 6.6743e-11)  # kg/m^3
rho_crit = unyt_quantity(rho_crit_si, "kg/m**3").to("Msun/kpc**3")

# Choose your halo
halo = obj.halos[caesar_group_id]  # Replace with desired index or search method

M500 = unyt_quantity(halo.virial_quantities['m500c'].in_units("Msun"), "Msun")
M200 = unyt_quantity(halo.virial_quantities['m200c'].in_units("Msun"), "Msun")

R500 = unyt_quantity(halo.virial_quantities['r500c'].in_units("m"), "m") # IN METERS FOR TEMP CALC
R200 = unyt_quantity(halo.virial_quantities['r200c'].in_units("m"), "m")

ne_500 = (500 * f_b * rho_crit / (mu_e * mp)).to("cm**-3")
ne_200 = (200 * f_b * rho_crit / (mu_e * mp)).to("cm**-3")

T500 = (G * M500 * mu * mp / R500).to("keV") # removing 1/kb here
T200 = (G * M200 * mu * mp / R200).to("keV") # removing 1/kb here

K500 = (T500 * ne_500**(-2/3)).to("keV*cm**2") # removing kb here
K200 = (T200 * ne_200**(-2/3)).to("keV*cm**2") # removing kb here

R500 = halo_data[caesar_group_id].virial_quantities['r500c'].in_units('Mpc') # redefine r500 and r200 with correct units before binning
R200 = halo_data[caesar_group_id].virial_quantities['r200c'].in_units('Mpc')
radial_bin_centres = (radial_bins[:-1] + radial_bins[1:])/2.
radial_bin_centres_logged = np.log10(radial_bin_centres)
log_radial_bin_centres = (np.log10(radial_bins[:-1]) + np.log10(radial_bins[1:]))/2.
log_radial_bin_centres_unlogged = 10**log_radial_bin_centres

normalized_r = radial_bin_centres/R500.value
normalized_k_mean = np.array(entropy_profile_mean)/K500.value
normalized_k_median = np.array(entropy_profile_median)/K500.value


# axes_.plot(CLoGS_lR,CLoGS_ly,color=CLoGS_fractional_color,lw=6,ls="--", zorder=obs_plot_kwargs['Oppenheimer+21']['zorder']+5,label="CLoGS Groups 1401")

log_r = np.log10(normalized_r)
log_k = np.log10(normalized_k_median)

log_m500c = np.log10(M500)
color = cmap(norm(log_m500c))

ax.plot(
    log_r,
    log_k,
    color='black',
    # color = cmap_sim(),
    lw=1.5,
    ls='--',
    zorder=60
    #label= f'Halo of M500: ' + str(round(np.log10(halo_data[caesar_group_id].virial_quantities['m500c'].value), 1)) + ', Cont: ' + str(round(halo_data[caesar_group_id].contamination, 3))
    #label=f'simulation {run+1}'

)





sns.set_context("paper", font_scale=1.3)
sns.set_style("whitegrid")

plt.minorticks_on()
plt.tick_params(direction='in', which='both', top=False, right=False)
plt.tight_layout()

plt.title(thetitle)
plt.xlabel(r'log($R / R_{500}$)')
plt.ylabel(r'$\log(K_{\mathrm{mass}} / K_{500})$')
# plt.legend(loc='upper left', fontsize=7)

cbar = plt.colorbar(sm, ax=ax)
cbar.set_label(r"$\log_{10}(M_{500c}$)")


plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig(savefilename, dpi=300, bbox_inches='tight')