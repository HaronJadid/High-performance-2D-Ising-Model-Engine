#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include "src/Lattice.h"

void runSimulation(int L) {
    std::cout << "--- Simulating Lattice L=" << L << " ---" << std::endl;

    Lattice lattice(L);
    lattice.initialize();

    std::string filename = "ising_L" + std::to_string(L) + ".csv";
    std::ofstream file(filename);
    file << "T,M,E,Cv,Chi\n";

    // Focus on the Critical Region
    double T_start = 2.6;
    double T_end = 2.0;
    double T_step = 0.02;

    for (double T = T_start; T >= T_end; T -= T_step) {
        lattice.precomputeExponentials(T);

        // Equilibration: Large systems need more time
        int equilibrate = 2000 + (L * 10);
        for (int i = 0; i < equilibrate; ++i) {
            lattice.metropolisSweepOMP(T);
        }

        // Data Collection
        long long E_sum = 0;
        long long M_sum = 0;
        long long E_sq_sum = 0;
        long long M_sq_sum = 0;

        int collect_steps = 10000; // High precision

        for (int i = 0; i < collect_steps; ++i) {
            lattice.metropolisSweepOMP(T);

            if (i % 5 == 0) {
                double E = lattice.computeTotalEnergy();
                double M = std::abs(lattice.getMagnetization());

                E_sum += E;
                M_sum += M;
                E_sq_sum += (E * E);
                M_sq_sum += (M * M);
            }
        }

        // Statistics
        int samples = collect_steps / 5;
        double E_mean = (double)E_sum / samples;
        double M_mean = (double)M_sum / samples;
        double E_var = ((double)E_sq_sum / samples) - (E_mean * E_mean);
        double M_var = ((double)M_sq_sum / samples) - (M_mean * M_mean);

        // Normalize
        double N_sq = (double)L * L;
        double Cv = E_var / (T * T * N_sq);
        double Chi = M_var / (T * N_sq);

        file << std::fixed << std::setprecision(4)
             << T << "," << (M_mean/N_sq) << "," << (E_mean/N_sq)
             << "," << Cv << "," << Chi << "\n";

        std::cout << "T=" << T << " done. \r" << std::flush;
    }
    std::cout << "\nSaved to " << filename << std::endl;
}

int main() {
    std::vector<int> sizes = {20, 40, 60, 80};
    for (int L : sizes) {
        runSimulation(L);
    }
    return 0;
}