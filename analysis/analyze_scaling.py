import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# --- Parameters ---
# Theoretical Critical Exponents for 2D Ising Model (Onsager Solution)
Tc_theory = 2.269
gamma_theory = 1.75
nu_theory = 1.0

# Lattice sizes
L_values = [20, 40, 60, 80]

# Setup Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle(f'Finite Size Scaling Analysis (2D Ising Model)', fontsize=16)

colors = ['b', 'g', 'r', 'm']

# Loop through data
for i, L in enumerate(L_values):
    filename = f"ising_L{L}.csv"

    if not os.path.exists(filename):
        print(f"Warning: {filename} not found in {os.getcwd()}")
        continue

    df = pd.read_csv(filename)
    df = df.sort_values('T') # Ensure sorted

    T = df['T']
    Chi = df['Chi']

    # --- Plot 1: Raw Data ---
    ax1.plot(T, Chi, 'o-', markersize=4, label=f'L={L}', color=colors[i])

    # --- Plot 2: Scaled Data (Data Collapse) ---
    # Reduced temperature: t = (T - Tc) / Tc
    t = (T - Tc_theory) / Tc_theory

    # X-axis: t * L^(1/nu)
    x_scaled = t * (L ** (1/nu_theory))

    # Y-axis: Chi / L^(gamma/nu)
    y_scaled = Chi / (L ** (gamma_theory / nu_theory))

    ax2.plot(x_scaled, y_scaled, 'o', markersize=4, alpha=0.6, label=f'L={L}', color=colors[i])

# Styling Plot 1 (Raw)
ax1.set_xlabel('Temperature ($T$)')
ax1.set_ylabel('Susceptibility ($\chi$)')
ax1.set_title('Raw Susceptibility Peaks')
ax1.axvline(Tc_theory, color='k', linestyle='--', alpha=0.3, label='Theoretical $T_c$')
ax1.legend()
ax1.grid(True)

# Styling Plot 2 (Collapsed)
ax2.set_xlabel(r'Scaled Temperature $t \cdot L^{1/\nu}$')
ax2.set_ylabel(r'Scaled Susceptibility $\chi / L^{\gamma/\nu}$')
ax2.set_title('Data Collapse (Verification of Universality)')
ax2.set_xlim(-2, 2) # Zoom in on the critical region
ax2.grid(True)
ax2.legend()

print("Plot generated! Check the popup window.")
plt.show()