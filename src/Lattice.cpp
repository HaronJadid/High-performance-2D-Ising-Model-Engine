//
// Created by harJD on 12/30/2025.
//
#include "Lattice.h"
#include <chrono>
#include <cmath> // for std::exp
#include <omp.h>

Lattice::Lattice(int size) : N(size), n_spins(size * size) {
    // Seed the Mersenne Twister with current time
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    rng = std::mt19937(seed);

    // Distribution to pick 0 or 1. We will map 0 -> -1 and 1 -> +1
    dist = std::uniform_int_distribution<int>(0, 1);

    spins.resize(n_spins);

    acceptance_dist = std::uniform_real_distribution<double>(0.0, 1.0);
}

void Lattice::initialize() {
    for (int i = 0; i < n_spins; ++i) {
        // Map 0 to -1, and 1 to 1
        spins[i] = (dist(rng) == 0) ? -1 : 1;
    }
}

int Lattice::getSpin(int x, int y) const {
    // Basic boundary check (can be removed later for raw speed if careful)
    if (x < 0 || x >= N || y < 0 || y >= N) return 0;
    return spins[x + y * N];
}

void Lattice::flipSpin(int x, int y) {
    spins[x + y * N] *= -1;
}

void Lattice::printLattice() const {
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            char c = (getSpin(x, y) == 1) ? '+' : '.';
            std::cout << c << " ";
        }
        std::cout << "\n";
    }
    std::cout << "-------------------\n";
}

int Lattice::getMagnetization() const {
    int M = 0;
    for (int s : spins) {
        M += s;
    }
    return M;
}

double Lattice::getNeighborSum(int x, int y) const {
    double sum = 0.0;

    // Right Neighbor ((x+1) % N)
    sum += getSpin((x + 1) % N, y);

    // Left Neighbor ((x - 1 + N) % N)
    sum += getSpin((x - 1 + N) % N, y);

    // Down Neighbor ((y + 1) % N)
    sum += getSpin(x, (y + 1) % N);

    // Up Neighbor ((y - 1 + N) % N)
    sum += getSpin(x, (y - 1 + N) % N);

    return sum;
}

double Lattice::computeTotalEnergy() const {
    double energy = 0.0;

    // Iterate over all sites
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            // H = -J * s_i * sum(s_neighbors)
            // But be careful! If we sum all 4 neighbors for every site,
            // we count every pair twice (i-j and j-i).
            // So we divide by 2.

            double neighbors = getNeighborSum(x, y);
            energy += -J * getSpin(x, y) * neighbors;
        }
    }

    return energy / 2.0;
}

void Lattice::metropolisSweep(double T) {
    // Ensure table is ready (You could also move this call to main loop)
    // For safety, let's assume precomputeExponentials(T) was called outside
    // or call it here if T changes (but that adds a check).
    // Let's assume the user calls precomputeExponentials in main.

    for (int i = 0; i < n_spins; ++i) {
        int x = std::uniform_int_distribution<int>(0, N - 1)(rng);
        int y = std::uniform_int_distribution<int>(0, N - 1)(rng);

        double neighborSum = getNeighborSum(x, y);
        int currentSpin = getSpin(x, y);

        // Calculate Delta E (Cast to int for array indexing)
        int deltaE = static_cast<int>(2.0 * J * currentSpin * neighborSum);

        if (deltaE <= 0) {
            flipSpin(x, y);
        } else {
            // OPTIMIZATION: Use the lookup table!
            // deltaE will be 4 or 8 here.
            if (acceptance_dist(rng) < prob_table[deltaE]) {
                flipSpin(x, y);
            }
        }
    }
}


void Lattice::precomputeExponentials(double T) {
    // We need to store probabilities for Delta E = 4 and 8.
    // Let's size the vector to 9 so we can access index 4 and 8 directly.
    // It wastes a tiny bit of memory (indices 0-3, 5-7 unused) but gains speed/readability.
    prob_table.resize(9);

    prob_table[4] = std::exp(-4.0 * J / T);
    prob_table[8] = std::exp(-8.0 * J / T);
}

void Lattice::metropolisSweepOMP(double T) {
    for (int parity = 0; parity < 2; ++parity) {

#pragma omp parallel
        {
            // OPTIMIZATION: static thread_local ensures initialization runs ONCE per thread
            static thread_local std::mt19937* local_rng = nullptr;

            // Lazy initialization
            if (!local_rng) {
                int tid = omp_get_thread_num();
                auto now = std::chrono::system_clock::now().time_since_epoch().count();
                local_rng = new std::mt19937(now + tid);
            }

            // Use *local_rng instead of local_rng below...

#pragma omp for schedule(static)
            for (int y = 0; y < N; ++y) {
                int start_x = (y % 2 == parity) ? 0 : 1;
                for (int x = start_x; x < N; x += 2) {

                    double neighborSum = getNeighborSum(x, y);
                    int currentSpin = getSpin(x, y);
                    int deltaE = static_cast<int>(2.0 * J * currentSpin * neighborSum);

                    bool flip = false;
                    if (deltaE <= 0) {
                        flip = true;
                    } else {
                        // DEREFERENCE the pointer here: (*local_rng)
                        if (acceptance_dist(*local_rng) < prob_table[deltaE]) {
                            flip = true;
                        }
                    }

                    if (flip) spins[x + y * N] *= -1;
                }
            }
        }
    }


}
void Lattice::saveState(std::ofstream& out) const {
    for (int s : spins) {
        out << s << " ";
    }
    out << "\n";
}