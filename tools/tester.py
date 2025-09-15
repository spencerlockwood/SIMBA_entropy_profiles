import numpy as np
import pynbody as pnb
import os
import argparse
import dill
import matplotlib.pyplot as plt
import caesar

from astro_constants import kB, mp  # Boltzmann constant & proton mass

# Argument parser
parser = argparse.ArgumentParser(description="Plot entropy profiles from simulation data.")
parser.add_argument('--snap', action='store', type=str, required=True, 
                    help='Path to snapshot file')
parser.add_argument('--caesar_catalog', action='store', type=str, required=True, 
                    help='Path to CAESAR file')
parser.add_argument('--output', action='store', type=str, required=False, default='entropy_profile.png',
                    help='Output file name for the plot')
parser.add_argument('--output_path', action='store', type=str, required=False, default='.',
                    help='Output file location/path')
args = parser.parse_args()

# Load the snapshot
data = pnb.load(args.snap)

# Load CAESAR catalog
caesar_data = caesar.load(args.caesar_catalog)

# Extract halo center from CAESAR
halo = caesar_data.halos[0]  # Assuming we are analyzing the first halo
halo_center = halo.pos

# Extracting gas properties
rho = data.gas['rho']  # Density
T = data.gas['temp']   # Temperature
mass = data.gas['mass']
r = np.sqrt(np.sum((data.gas['pos'] - halo_center)**2, axis=1))  # Radial distance from halo center

# Sort by radius for profile calculation
sorted_indices = np.argsort(r)
r_sorted = r[sorted_indices]
rho_sorted = rho[sorted_indices]
T_sorted = T[sorted_indices]
mass_sorted = mass[sorted_indices]

# Define radial bins
num_bins = 50
r_bins = np.linspace(np.min(r_sorted), np.max(r_sorted), num_bins+1)

# Compute entropy profile
entropy_profile = []
r_bin_centers = []

for i in range(num_bins):
    in_bin = (r_sorted >= r_bins[i]) & (r_sorted < r_bins[i+1])
    if np.sum(in_bin) > 0:
        mass_weighted_entropy = np.sum((kB * T_sorted[in_bin]) / (rho_sorted[in_bin]**(2/3)) * mass_sorted[in_bin]) / np.sum(mass_sorted[in_bin])
        entropy_profile.append(mass_weighted_entropy)
        r_bin_centers.append((r_bins[i] + r_bins[i+1]) / 2)

# Convert lists to arrays
r_bin_centers = np.array(r_bin_centers)
entropy_profile = np.array(entropy_profile)

# Plot the entropy profile
plt.figure(figsize=(8, 6))
plt.plot(r_bin_centers, entropy_profile, marker='o', linestyle='-', color='b', label='Entropy Profile')
plt.xlabel('Radius [kpc]')
plt.ylabel('Entropy [keV cm^2]')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.grid()
plt.title('Mass-Weighted Entropy Profile')
output_filepath = os.path.join(args.output_path, args.output)
plt.savefig(output_filepath)
plt.show()
