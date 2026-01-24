#include "Lattice.h"
#include <chrono>
#include <cmath>
#include <omp.h>

Lattice::Lattice(int size) : N(size), n_spins(size * size) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    rng = std::mt19937(seed);
    dist = std::uniform_int_distribution<int>(0, 1);
    spins.resize(n_spins);
    acceptance_dist = std::uniform_real_distribution<double>(0.0, 1.0);
}

void Lattice::initialize() {
    for (int i = 0; i < n_spins; ++i) {
        spins[i] = (dist(rng) == 0) ? -1 : 1;
    }
}

int Lattice::getSpin(int x, int y) const {
    if (x < 0 || x >= N || y < 0 || y >= N) return 0;
    return spins[x + y * N];
}

void Lattice::flipSpin(int x, int y) {
    spins[x + y * N] *= -1;
}

// --- MISSING FUNCTION ADDED HERE ---
int Lattice::getMagnetization() const {
    int M = 0;
    for (int s : spins) {
        M += s;
    }
    return M;
}
// -----------------------------------

// Optimized Neighbor Sum (No Modulo)
double Lattice::getNeighborSum(int x, int y) const {
    double sum = 0.0;
    int row_offset = y * N;

    // Right: (x + 1) % N
    int right = (x + 1 == N) ? 0 : x + 1;
    sum += spins[right + row_offset];

    // Left: (x - 1 + N) % N
    int left = (x == 0) ? N - 1 : x - 1;
    sum += spins[left + row_offset];

    // Up: (y - 1 + N) % N
    int up = (y == 0) ? N - 1 : y - 1;
    sum += spins[x + up * N];

    // Down: (y + 1) % N
    int down = (y + 1 == N) ? 0 : y + 1;
    sum += spins[x + down * N];

    return sum;
}

double Lattice::computeTotalEnergy() const {
    double energy = 0.0;
    for (int y = 0; y < N; ++y) {
        for (int x = 0; x < N; ++x) {
            double neighbors = getNeighborSum(x, y);
            energy += -J * spins[x + y * N] * neighbors;
        }
    }
    return energy / 2.0;
}

void Lattice::precomputeExponentials(double T) {
    prob_table.resize(9);
    prob_table[4] = std::exp(-4.0 * J / T);
    prob_table[8] = std::exp(-8.0 * J / T);
}

// Removed unused 'T' parameter name to fix warning
void Lattice::metropolisSweepOMP(double /*T*/) {
    for (int parity = 0; parity < 2; ++parity) {
#pragma omp parallel
        {
            static thread_local std::mt19937* local_rng = nullptr;
            if (!local_rng) {
                int tid = omp_get_thread_num();
                auto now = std::chrono::system_clock::now().time_since_epoch().count();
                local_rng = new std::mt19937(now + tid);
            }

#pragma omp for schedule(static)
            for (int y = 0; y < N; ++y) {
                int start_x = (y % 2 == parity) ? 0 : 1;
                for (int x = start_x; x < N; x += 2) {
                    double neighborSum = getNeighborSum(x, y);
                    int currentSpin = spins[x + y * N];
                    int deltaE = static_cast<int>(2.0 * J * currentSpin * neighborSum);

                    if (deltaE <= 0) {
                        spins[x + y * N] *= -1;
                    } else {
                        if (acceptance_dist(*local_rng) < prob_table[deltaE]) {
                            spins[x + y * N] *= -1;
                        }
                    }
                }
            }
        }
    }
}

void Lattice::saveState(std::ofstream& out) const {
    for (int i = 0; i < n_spins; ++i) {
        out << spins[i] << (i == n_spins - 1 ? "" : " ");
    }
    out << "\n";
}