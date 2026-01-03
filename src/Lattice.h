//
// Created by harJD on 12/30/2025.
//

#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <random> // Crucial for Mersenne Twister
#include <iostream>
#include <fstream>

class Lattice {
private:
    int N; // Dimension of the grid (NxN)
    int n_spins; // Total number of spins (N*N)
    std::vector<int> spins; // Flattened 1D array for Cache Locality

    // Random Number Generation Engine
    std::mt19937 rng;
    std::uniform_int_distribution<int> dist;
    double J = 1.0;
    std::uniform_real_distribution<double> acceptance_dist;
    std::vector<double> prob_table;

public:
    // Constructor
    Lattice(int size);

    // Initialize the lattice with random spins (+1 or -1)
    void initialize();

    // Accessors (Inline for performance)
    int getSpin(int x, int y) const;
    int getN() const { return N; }

    // Flip a spin at a specific location
    void flipSpin(int x, int y);

    // Visualization (for debugging)
    void printLattice() const;

    // Calculate total magnetization
    int getMagnetization() const;

    double computeTotalEnergy() const;

    double getNeighborSum(int x, int y) const;

    void metropolisSweep(double T);

    void precomputeExponentials(double T);

    void metropolisSweepOMP(double T);

    void saveState(std::ofstream& out) const;

};

#endif //LATTICE_H