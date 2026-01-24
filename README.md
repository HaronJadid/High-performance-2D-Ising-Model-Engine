# Discrete 2D Ising Model Monte Carlo Simulation Engine

##  Key Features

* **HPC Architecture:** Flattened 1D memory layout for optimal CPU cache locality.
* **Parallelism:** Thread-safe implementation using OpenMP and Red-Black (Checkerboard) domain decomposition.
* **Optimization:** Precomputed Boltzmann lookup tables to bypass expensive transcendental functions.
* **Scientific Accuracy:** Uses **Mersenne Twister (mt19937)** for high-quality entropy (thread-local seeding).
* **Reproducibility:** Fully Dockerized environment for consistent scientific results.

## Scientific Results

## Build

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

```bash
./IsingModelSim
```
