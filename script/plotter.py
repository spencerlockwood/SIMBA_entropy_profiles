import dill
import matplotlib.pyplot as plt
import numpy as np

# Load the profile from the .pkl file
with open("/home/b/babul/stlock/entropy_profiles/outputs/entropy_prof_L1_zoom_252_644_temp_cut=5e5K_nh_cut=0.13_cm-3.pkl", "rb") as f:
    profile = dill.load(f)

# Extract data
r = profile['x']['bin_centres']  # Radial bin centers
entropy = profile['y']['K_mass']  # Entropy

# Convert to log scale if necessary
log_r = np.log10(r)
log_entropy = np.log10(entropy)

# Plot
plt.figure(figsize=(8,6))
plt.plot(r, entropy, marker='o', linestyle='-', label="Entropy Profile")
plt.xscale("log")  # Log scale for radius
plt.yscale("log")  # Log scale for entropy (optional)
plt.xlabel("Radius (kpc)")
plt.ylabel("Entropy (keV cm²)")
plt.title("Entropy vs. Radius")
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.show()
