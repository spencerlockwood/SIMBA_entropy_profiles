# Simulation data
target_halo_id = 4

# Plotting/saving
#savefilename = 'Snap' + str(zoom_ID) + '_Halo_contamination.png' # .png name to save plot under
savefilename = '/home/b/babul/stlock/entropy_profiles/outputs/BAL_runs_energies_testing.png'
thetitle = "Halo 3224 Simba-C L0 Energy - Varied level of gas accreted by BH's" # filepath/title of graph


# LOWEST -> HIGHEST BlackHoleFeedbackFactor

# (with 5 and 6) order = [20, 19, 18, 17, 26, 16, 25, 15, 24, 14, 23, 13, 22, 12, 21, 1, 27, 2, 28, 3, 29, 4, 30, 5, 31, 6, 32, 7, 8, 9, 10, 11]
order = [0, 1]
num_runs = 10


# -------------------------------------------------------------------------------------------------------------------------------------- SIM_DATA
import yt
import unyt
from unyt import unyt_quantity, unyt_array, km, s, Mpc, kg, m, Msun, kpc
from unyt import Gyr, Mpc, km, s
from unyt.physical_constants import G as G_unyt, kboltz, mp as mp_unyt  # G, k_B, m_p

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D  # For custom legend handles

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

#norm = mcolors.Normalize(vmin=low_mass_range, vmax=high_mass_range)
#sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

cmap_total = plt.cm.get_cmap('coolwarm', 9)
cmap_internal = plt.cm.get_cmap('coolwarm', 9)
# ------------------------------------------------------------------------------------------------------------------ Loop for each snap/caesar file!!!
color_n = 7
color_index = 1
# cmap = plt.cm.get_cmap('magma', 9)

fig, ax = plt.subplots(figsize=(8, 6))

for run in range(num_runs):
    if run+1!=6:
        target_halo_id = 4
        if run+1==6:
            target_halo_id=5
        

        snap_file = f'/scratch/b/babul/stlock/snapshots/HyenasC/L0/halo_3224_finished_runs/BAL_variations/run{run+1:02}/snapshot_151.hdf5'
        caesarfilepath = f'/scratch/b/babul/stlock/snapshots/HyenasC/L0/halo_3224_finished_runs/BAL_variations/run{run+1:02}/caesarfile_snap151_test.hdf5'

        # ADD CAESAR FILEPATH HERE TO 

        snap = yt.load(snap_file)#, unit_system='cgs')
        data = snap.all_data()
        s = pnb.load(snap_file)

        # Get data, and some functions

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

        obj = caesar.load(caesarfilepath)
        halo_data = obj.halos

        halo_center = halo_data[target_halo_id].minpotpos
        halo_center = halo_center.in_units('Mpc')
        print(halo_center)


        igrm = data.include_above(('PartType0', 'temperature'), 5e5, units='K')
        igrm = igrm.include_below(('PartType0', 'H_nuclei_density'), 0.13, units='cm**-3')
        igrm = igrm.include_equal(('PartType0', 'DelayTime'), 0)


        # Radial Binning

        # radial_bins = unyt.unyt_array(np.append([0], np.logspace(-3, np.log10(5), 30)), 'Mpc')
        radial_bins = unyt.unyt_array(np.append([0], np.logspace(-4, np.log10(5), 45)), 'Mpc')

        halo_center.in_units('Mpc')

        igrm_pos_x = igrm['gas', 'particle_position_x']
        igrm_pos_y = igrm['gas', 'particle_position_y']
        igrm_pos_z = igrm['gas', 'particle_position_z']

        # isolate igrm

        igrm_halo_center = ((igrm_pos_x-halo_center[0])**2 + (igrm_pos_y-halo_center[1])**2 + (igrm_pos_z-halo_center[2])**2)**0.5

        igrm_zoom_entropies = igrm['gas', 'entropy']
        igrm_zoom_masses = igrm['gas', 'mass']


        entropy_profile_median = []
        entropy_profile_mean = []
        entropy_profile_mass_weighted = []

        internal_energy_profile = []
        kinetic_energy_profile = []
        total_energy_profile = []
        for bin_lo, bin_hi in zip(radial_bins[:-1], radial_bins[1:]):
            radial_filter = np.logical_and(igrm_halo_center.in_units('Mpc').value>=bin_lo.in_units('Mpc').value,
                                        igrm_halo_center.in_units('Mpc').value<=bin_hi.in_units('Mpc').value)
            
            entropies = igrm_zoom_entropies[radial_filter]

            # -------------------------------------------------------------------------------energy profile
            # Get particle data with proper units
            masses = igrm_zoom_masses[radial_filter].to('g')  # Convert to grams
            vx = igrm['gas','velocity_x'][radial_filter].in_units('cm/s')  # Explicitly convert to cm/s
            vy = igrm['gas','velocity_y'][radial_filter].in_units('cm/s')
            vz = igrm['gas','velocity_z'][radial_filter].in_units('cm/s')
            v_squared = (vx**2 + vy**2 + vz**2)  # Now properly in cm²/s²

            temps = igrm['gas','temperature'][radial_filter].to('K')
            # 3. Compute specific energies (keV/g)
            specific_ke = (0.5 * v_squared).to('keV/g')      # Kinetic energy per gram
            specific_u_int = (1.5 * kboltz.to('keV/K') * temps) / mp_unyt.to('g')  # Internal energy per gram
            
            # 4. Mass-weighted averages for the bin
            if len(masses) > 0:
                avg_u_int = np.average(specific_u_int, weights=masses)
                avg_u_tot = np.average(specific_u_int + specific_ke, weights=masses)
            else:
                avg_u_int = np.nan
                avg_u_tot = np.nan
            
            # 5. Append to profiles
            internal_energy_profile.append(avg_u_int)
            total_energy_profile.append(avg_u_tot)
            # -------------------------------------------------------------------------------

        from unyt import Msun, kpc, cm, keV, unyt_quantity
        from unyt.physical_constants import mp, kboltz, G  # ✅ use unyt's version


        # Constants
        mu = 0.59        # Mean molecular weight
        mu_e = 1.14      # Mean molecular weight per free electron

        # Load CAESAR catalog
        h = obj.simulation.hubble_constant
        omega_b = obj.simulation.omega_baryon
        omega_m = obj.simulation.omega_matter

        # Update f_b if desired
        f_b = omega_b / omega_m

        # Critical density
        H0_si = (100 * h * 1e3) / (3.086e22)  # H0 in 1/s
        rho_crit_si = 3 * H0_si**2 / (8 * np.pi * 6.6743e-11)  # kg/m^3
        rho_crit = unyt_quantity(rho_crit_si, "kg/m**3").to("Msun/kpc**3")

        # Choose your halo
        halo = obj.halos[target_halo_id]  # Replace with desired index or search method

        M500 = unyt_quantity(halo.virial_quantities['m500c'].in_units("Msun"), "Msun")

        R500 = unyt_quantity(halo.virial_quantities['r500c'].in_units("m"), "m") # in M for T500 calc

        n_500 = (500 * f_b * rho_crit / (mu * mp)).to("cm**-3")
        # T500 = (G * M500 * mu * mp / R500).to("keV") # removing 1/kb here
        R500 = halo_data[target_halo_id].virial_quantities['r500c'].in_units('Mpc')
        
        
        radial_bin_centres = (radial_bins[:-1] + radial_bins[1:])/2.
        radial_bin_centres_logged = np.log10(radial_bin_centres)
        log_radial_bin_centres = (np.log10(radial_bins[:-1]) + np.log10(radial_bins[1:]))/2.
        log_radial_bin_centres_unlogged = 10**log_radial_bin_centres

        normalized_r = radial_bin_centres/R500.value

        # -------------------------------------------------------------------------------energy profile
        # Calculate T500 (as kT in keV)
        T500_keV = (G * M500 * mu * mp / R500).to("keV")

        # Convert to actual temperature (K)
        T500_K = (T500_keV / kboltz.to('keV/K')).to('K')

        # Calculate n_500
        n_500 = (500 * f_b * rho_crit / (mu * mp)).to("cm**-3")

        # Calculate U500 (energy density at R500)
        U500 = (1.5 * n_500 * kboltz.to('keV/K') * T500_K).to('keV/cm**3')

        # Now normalize your profiles
        U500_specific = (1.5 * kboltz.to('keV/K') * T500_K) / mp_unyt.to('g')  # keV/g
        normalized_u_int = np.array(internal_energy_profile) / U500_specific
        normalized_u_tot = np.array(total_energy_profile) / U500_specific

        non_normalized_u_int = np.array(internal_energy_profile)
        non_normalized_u_tot = np.array(total_energy_profile)


        # Convert to log scale for plotting
        log_u_int = np.log10(normalized_u_int)
        log_u_tot = np.log10(normalized_u_tot)

        # n_log_u_int = np.log10(non_normalized_u_int)
        # n_log_u_tot = np.log10(non_normalized_u_tot)

        # -------------------------------------------------------------------------------energy profile


        log_r = np.log10(normalized_r)
        # Then modify your plotting section to include these new profiles:
        """
            plt.plot(
            log_r,
            log_k,
            color=cmap_sim(run),
            lw=1,
            ls='-',
            zorder=60,
            label=f'simulation {run+1} (Entropy)'
        )
        """

        plt.plot(
            log_r,
            log_u_int,
            color=cmap_internal(run),
            lw=1,
            ls='--',
            zorder=60,
            label=f'Run {run+1} Internal'
        )
        
        plt.plot(
            log_r,
            log_u_tot,
            color=cmap_total(run),
            lw=1,
            ls='-',
            zorder=60,
            label=f'Run {run+1} Total'
        )
        
        color_index+=1
color_index = 1
# ------------------------------------------------------------------------------------------------------------------ end Loop!!!


#====================================================================================================================================================


sns.set_context("paper", font_scale=1.3)
sns.set_style("whitegrid")

"""
note = 'Observational Data Mass Range [M500]: ' + str(low_mass_range) + ' - ' + str(high_mass_range) + ', z = ' + str(round(snap.current_redshift, 3))

plt.text(
    0.98, 0.02,                # x, y (in axes fraction coordinates)
    note,            # your note
    transform=plt.gca().transAxes,  # place relative to axis (not data)
    fontsize=8, 
    ha='right', 
    va='bottom'
)


plt.minorticks_on()
plt.tick_params(direction='in', which='both', top=False, right=False)
plt.tight_layout()

plt.title(thetitle)
plt.xlabel(r'log($R / R_{500}$)')
plt.ylabel(r'KeV/g')
plt.legend(loc='upper left', fontsize=7)

plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.savefig(savefilename, dpi=300, bbox_inches='tight')
"""

# Create custom legend entries
legend_handles = [
    # Style descriptors (Internal/Total)
    Line2D([], [], color='gray', ls='--', lw=2, label='Internal Energy'),
    Line2D([], [], color='gray', ls='-', lw=2, label='Total Energy'),
    
    # Add colored run entries
    *[Line2D([], [], color=cmap_internal(run), lw=2, marker='o', ls='None', label=f'Run {run+1}') 
      for run in range(num_runs)]
]

cmap_bar = plt.cm.get_cmap('coolwarm', 1000)
norm = mpl.colors.Normalize(vmin=0.02, vmax=0.20)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap_bar)

cbar = plt.colorbar(sm, ax=ax, pad=0.02)

# Label
cbar.set_label(
    'BAL_f_accretion value',
    rotation=270,
    labelpad=0,
    fontsize=8
)

# Only show lowest and highest values
cbar.set_ticks([0.02, 0.20])
cbar.set_ticklabels(['0.02', '0.20'])

# Remove tick marks (optional but clean)
cbar.ax.tick_params(size=0)


# Add the legend
plt.legend(handles=legend_handles, loc='upper left', fontsize=8, ncol=1)

# Rest of your formatting
plt.minorticks_on()
plt.tick_params(direction='in', which='both', top=False, right=False)
plt.title(thetitle)
plt.xlabel(r'log($R / R_{500}$)')
plt.ylabel(r'log($U/U_{500}$)')
# plt.ylabel('Energy [Kev/g]')
plt.grid(True, which='both', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig(savefilename, dpi=300, bbox_inches='tight')