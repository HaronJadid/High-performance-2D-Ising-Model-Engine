#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include "Lattice.h"

void runScalingAnalysis() {
    std::cout << "=== Running Scaling Analysis (Data Collapse) ===" << std::endl;
    std::vector<int> sizes = {20, 40, 60, 80};

    for (int L : sizes) {
        std::cout << "--- Simulating Lattice L=" << L << " ---" << std::endl;
        Lattice lattice(L);
        lattice.initialize();

        std::string filename = "ising_L" + std::to_string(L) + ".csv";
        std::ofstream file(filename);
        file << "T,M,E,Cv,Chi\n";

        double T_start = 2.6, T_end = 2.0, T_step = 0.02;

        for (double T = T_start; T >= T_end; T -= T_step) {
            lattice.precomputeExponentials(T);

            // Equilibration
            for (int i = 0; i < 2000; ++i) lattice.metropolisSweepOMP(T);

            // Collection
            long long M_sum = 0, M_sq_sum = 0;
            int collect_steps = 5000;

            for (int i = 0; i < collect_steps; ++i) {
                lattice.metropolisSweepOMP(T);
                double M = std::abs(lattice.getMagnetization());
                M_sum += M;
                M_sq_sum += (M * M);
            }

            double M_mean = (double)M_sum / collect_steps;
            double M_var = ((double)M_sq_sum / collect_steps) - (M_mean * M_mean);
            double N_sq = (double)L * L;
            double Chi = M_var / (T * N_sq);

            file << T << "," << (M_mean/N_sq) << ",0,0," << Chi << "\n";
            std::cout << "T=" << T << " Chi=" << Chi << "      \r" << std::flush;
        }
        std::cout << "\nSaved " << filename << std::endl;
    }
}

void runVisualization() {
    std::cout << "=== Running Visualization (Domain Growth) ===" << std::endl;
    int L = 200;
    Lattice lattice(L);
    lattice.initialize(); // Infinite T state (Random)

    double T_quench = 1.0; // Deep ferromagnetic phase
    lattice.precomputeExponentials(T_quench);

    std::ofstream out("snapshots.txt");
    std::cout << "Quenching to T=" << T_quench << ". Saving snapshots..." << std::endl;

    // Save Initial State
    lattice.saveState(out);

    int total_steps = 1000;
    int save_interval = 5;

    for (int i = 0; i < total_steps; ++i) {
        lattice.metropolisSweepOMP(T_quench);
        if (i % save_interval == 0) {
            lattice.saveState(out);
        }
    }
    std::cout << "Done. Snapshots saved to 'snapshots.txt'." << std::endl;
    std::cout << "Run 'python3 analysis/visualize.py' to generate GIF." << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc > 1 && std::strcmp(argv[1], "--viz") == 0) {
        runVisualization();
    } else {
        runScalingAnalysis();
        std::cout << "\nHint: Run with '--viz' to generate domain growth animation." << std::endl;
    }
    return 0;
}