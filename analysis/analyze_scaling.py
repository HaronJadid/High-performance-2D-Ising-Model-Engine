import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Configuration
Tc_theory = 2.269
gamma_theory = 1.75
nu_theory = 1.0
L_values = [20, 40, 60, 80]
OUTPUT_DIR = "analysis/plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle(f'Finite Size Scaling Analysis (2D Ising Model)', fontsize=16)

colors = ['b', 'g', 'r', 'm']

for i, L in enumerate(L_values):
    filename = f"ising_L{L}.csv"
    if not os.path.exists(filename):
        print(f"Warning: {filename} not found.")
        continue

    df = pd.read_csv(filename)
    df = df.sort_values('T')

    T = df['T']
    Chi = df['Chi']

    # Raw Data
    ax1.plot(T, Chi, 'o-', markersize=4, label=f'L={L}', color=colors[i])

    # Scaled Data
    t = (T - Tc_theory) / Tc_theory
    x_scaled = t * (L ** (1/nu_theory))
    y_scaled = Chi / (L ** (gamma_theory / nu_theory))
    ax2.plot(x_scaled, y_scaled, 'o', markersize=4, alpha=0.6, label=f'L={L}', color=colors[i])

ax1.set_xlabel('Temperature ($T$)')
ax1.set_ylabel('Susceptibility ($\chi$)')
ax1.set_title('Raw Susceptibility Peaks')
ax1.axvline(Tc_theory, color='k', linestyle='--', alpha=0.3, label='Theoretical $T_c$')
ax1.legend()
ax1.grid(True)

ax2.set_xlabel(r'Scaled Temperature $t \cdot L^{1/\nu}$')
ax2.set_ylabel(r'Scaled Susceptibility $\chi / L^{\gamma/\nu}$')
ax2.set_title('Data Collapse (Verification of Universality)')
ax2.set_xlim(-2, 2)
ax2.grid(True)
ax2.legend()

# Save the plot
output_path = os.path.join(OUTPUT_DIR, "data_collapse.png")
plt.savefig(output_path, dpi=300)
print(f"Plot saved to {output_path}")