import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

# Configuration
INPUT_FILE = "snapshots.txt"
OUTPUT_FILE = "domain_growth.gif"
L = 200  # Must match the C++ simulation

if not os.path.exists(INPUT_FILE):
    print(f"Error: {INPUT_FILE} not found. Run the C++ engine with --viz first.")
    exit(1)

print("Loading data...")
snapshots = []
with open(INPUT_FILE, 'r') as f:
    for line in f:
        # Parse space-separated spins
        spins = np.fromstring(line, dtype=int, sep=' ')
        snapshots.append(spins.reshape((L, L)))

print(f"Loaded {len(snapshots)} frames. Generating animation...")

fig = plt.figure(figsize=(6, 6))
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)

# Initial Frame
im = ax.imshow(snapshots[0], cmap='binary', interpolation='nearest', animated=True)

def update(i):
    im.set_array(snapshots[i])
    return [im]

ani = animation.FuncAnimation(fig, update, frames=len(snapshots), interval=50, blit=True)
ani.save(OUTPUT_FILE, writer='pillow', fps=20)

print(f"Success! Animation saved to {OUTPUT_FILE}")