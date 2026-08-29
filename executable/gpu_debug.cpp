// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

// Compares the GPU distance histogram kernel against the CPU kernel it is meant to replace,
// and reports the timing of both. Usage: gpu_debug [structure file]

#include <data/Molecule.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/SimpleCPU.h>
#include <settings/All.h>
#include <gpu/WebGPUSimple.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using namespace ausaxs;

namespace {
    struct Comparison {
        double max_absolute = 0;
        double max_relative = 0;
        int max_relative_bin = -1;
        double worst_bin_content = 0;
        double cpu_total = 0, gpu_total = 0;
        double absolute_total = 0; // summed bin-by-bin deviation, to catch contributions moving between bins
    };

    template<bool weighted_bins>
    Comparison compare(
        const typename hist::GenericDistribution1D<weighted_bins>::type& cpu,
        const typename hist::GenericDistribution1D<weighted_bins>::type& gpu
    ) {
        Comparison result;
        for (int i = 0; i < static_cast<int>(cpu.size()); ++i) {
            double c = cpu.get_content(i);
            double g = gpu.get_content(i);
            result.cpu_total += c;
            result.gpu_total += g;

            double absolute = std::abs(c - g);
            result.max_absolute = std::max(result.max_absolute, absolute);
            result.absolute_total += absolute;
            if (c != 0) {
                double relative = absolute/std::abs(c);
                if (relative > result.max_relative) {
                    result.max_relative = relative;
                    result.max_relative_bin = i;
                    result.worst_bin_content = c;
                }
            }
        }
        return result;
    }

    void print(std::string_view name, const Comparison& c) {
        std::cout << "    " << std::left << std::setw(6) << name
                  << " sum " << std::scientific << std::setprecision(3) << c.cpu_total << " vs " << c.gpu_total
                  << " (deviation " << std::abs(c.cpu_total - c.gpu_total)/std::abs(c.cpu_total)
                  << ", summed per-bin deviation " << c.absolute_total/std::abs(c.cpu_total) << ")" << std::endl;
        std::cout << "           worst bin " << c.max_relative_bin << ": " << c.max_relative
                  << " relative, " << c.max_absolute << " of " << c.worst_bin_content << std::endl;
    }

    template<bool weighted_bins>
    void run(const data::Molecule& molecule) {
        constexpr bool vbw = false;
        hist::detail::CompactCoordinates<vbw> atoms(molecule.get_bodies());
        hist::detail::CompactCoordinates<vbw> waters(molecule.get_waters());
        std::cout << (weighted_bins ? "weighted bins" : "unweighted bins")
                  << " (" << atoms.size() << " atoms, " << waters.size() << " waters)" << std::endl;

        auto time = [] (auto&& f) {
            auto start = std::chrono::steady_clock::now();
            auto result = f();
            auto end = std::chrono::steady_clock::now();
            return std::make_pair(std::move(result), std::chrono::duration<double, std::milli>(end - start).count());
        };

        auto [cpu_result, cpu_ms] = time([&] () {
            hist::distance_calculator::SimpleCPU<weighted_bins, vbw> calculator;
            calculator.enqueue_calculate_self(atoms);
            calculator.enqueue_calculate_cross(atoms, waters);
            return calculator.run();
        });

        auto [gpu_result, gpu_ms] = time([&] () {
            hist::distance_calculator::WebGPUSimple<weighted_bins, vbw> calculator;
            calculator.enqueue_calculate_self(atoms);
            calculator.enqueue_calculate_cross(atoms, waters);
            return calculator.run();
        });

        std::cout << "    cpu: " << std::fixed << std::setprecision(1) << cpu_ms << " ms"
                  << "  gpu: " << gpu_ms << " ms"
                  << "  speedup: " << std::setprecision(2) << (gpu_ms == 0 ? 0 : cpu_ms/gpu_ms) << "x" << std::endl;
        print("self", compare<weighted_bins>(cpu_result.self.at(0), gpu_result.self.at(0)));
        print("cross", compare<weighted_bins>(cpu_result.cross.at(0), gpu_result.cross.at(0)));
    }
}

int main(int argc, char const *argv[]) {
    std::string file = argc > 1 ? argv[1] : "tests/files/2epe.pdb";
    settings::general::verbose = false;

    data::Molecule molecule(file);
    molecule.generate_new_hydration();

    run<false>(molecule);
    run<true>(molecule);
    return 0;
}
