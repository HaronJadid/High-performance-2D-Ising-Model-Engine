# Discrete 2D Ising Model Monte Carlo Simulation Engine

High-performance C++ engine for simulating the 2D Ising Model using the Metropolis algorithm.

## Structure

- `src/`: Source code (C++).
- `analysis/`: Python scripts for data analysis and plotting.
- `data/`: Directory for storing simulation data.

## Quick Start (Docker)

Run the full simulation and generate the visualization in one command:

```bash
# Build the container
docker build -t ising-hpc .

# Run simulation & generate GIF (Outputs to mapped volume)
docker run -v $(pwd):/app/output ising-hpc cp domain_growth.gif /app/output/
```
## Visualization
To observe domain coarsening (Spontaneous Symmetry Breaking):

- Run `./IsingModelSim --viz` (Generates snapshots.txt)

- Run `python3 analysis/visualize.py`

- Result: `domain_growth.gif`

## Optimization Details
- Bitwise Neighbor Lookup: Replaced modulo operators with conditional logic for PBC.

- Cache Aligned Arrays: std::vector<int> flat layout ensures contiguous memory access.