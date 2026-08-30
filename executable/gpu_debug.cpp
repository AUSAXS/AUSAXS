// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

// Compares the GPU distance histogram kernels against the CPU kernel they are meant to replace,
// and reports the timing of all of them.
// Usage: gpu_debug [structure file] [threads] [repeats] [backends]
// The backends are given as a comma-separated subset of "sycl,webgpu", defaulting to all of them. The
// CPU kernel is always run, as the GPU results are compared against it.

#include <data/Molecule.h>
#include <hist/detail/CompactCoordinates.h>
#include <hist/distance_calculator/SimpleCPU.h>
#include <settings/All.h>

#ifdef AUSAXS_GPU
    #include <gpu/WebGPUSimple.h>
#endif

#ifdef AUSAXS_GPU_SYCL
    #include <gpu/SYCLSimple.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace ausaxs;

namespace {
    int repeats = 3;

    // which backends to measure. they share the device, and measuring one right after the other gives
    // noticeably noisier numbers than running each on its own, so it has to be possible to pick.
    bool run_sycl = true, run_webgpu = true;

    struct Comparison {
        double max_absolute = 0;
        double max_relative = 0;
        int max_relative_bin = -1;
        double worst_bin_content = 0;
        double cpu_total = 0, gpu_total = 0;
        double absolute_total = 0; // summed bin-by-bin deviation, to catch contributions moving between bins
        double peak_bin = 0;       // largest single bin, to judge the headroom of the fixed-point range
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
            result.peak_bin = std::max(result.peak_bin, std::abs(c));
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
        std::cout << "        " << std::left << std::setw(6) << name
                  << " sum " << std::scientific << std::setprecision(3) << c.cpu_total << " vs " << c.gpu_total
                  << " (deviation " << std::abs(c.cpu_total - c.gpu_total)/std::abs(c.cpu_total)
                  << ", summed per-bin deviation " << c.absolute_total/std::abs(c.cpu_total) << ")" << std::endl;
        std::cout << "               worst bin " << c.max_relative_bin << ": " << c.max_relative
                  << " relative, " << c.max_absolute << " of " << c.worst_bin_content
                  << "; peak bin " << c.peak_bin << std::endl;
    }

    /**
     * @brief Run @a f the requested number of times and return its result along with the shortest time taken.
     *
     * It is run untimed for a while first. Beyond paying for any one-time setup, this is what makes the
     * numbers comparable at all: a few milliseconds of work on an otherwise idle GPU is over before it
     * has left its low-power state, and the same kernel then measures anywhere within a factor of two
     * depending on what happened to run before it.
     */
    template<typename functor_t>
    auto time(functor_t&& f) {
        constexpr auto warmup = std::chrono::milliseconds(500);

        double best = std::numeric_limits<double>::max();
        auto result = f();
        for (auto end = std::chrono::steady_clock::now() + warmup; std::chrono::steady_clock::now() < end;) {
            result = f();
        }
        for (int i = 0; i < repeats; ++i) {
            auto start = std::chrono::steady_clock::now();
            result = f();
            auto end = std::chrono::steady_clock::now();
            best = std::min(best, std::chrono::duration<double, std::milli>(end - start).count());
        }
        return std::make_pair(std::move(result), best);
    }

    template<typename kernel_t, bool weighted_bins, bool vbw>
    void run_kernel(
        std::string_view name,
        const hist::detail::CompactCoordinates<vbw>& atoms,
        const hist::detail::CompactCoordinates<vbw>& waters,
        const typename hist::distance_calculator::SimpleKernel<weighted_bins, vbw>::run_result& cpu_result,
        double pairs, double cpu_ms
    ) {
        auto [result, ms] = time([&] () {
            kernel_t calculator;
            calculator.enqueue_calculate_self(atoms);
            calculator.enqueue_calculate_cross(atoms, waters);
            return calculator.run();
        });

        std::cout << "    " << std::left << std::setw(7) << name << std::right
                  << std::fixed << std::setprecision(1) << std::setw(8) << ms << " ms"
                  << "   speedup vs cpu: " << std::setprecision(2) << (ms == 0 ? 0 : cpu_ms/ms) << "x"
                  << "   (" << std::setprecision(1) << 1e-6*pairs/ms << " Gpair/s)" << std::endl;
        print("self", compare<weighted_bins>(cpu_result.self.at(0), result.self.at(0)));
        print("cross", compare<weighted_bins>(cpu_result.cross.at(0), result.cross.at(0)));
    }

    template<bool weighted_bins>
    void run(const data::Molecule& molecule) {
        constexpr bool vbw = false;
        hist::detail::CompactCoordinates<vbw> atoms(molecule.get_bodies());
        hist::detail::CompactCoordinates<vbw> waters(molecule.get_waters());
        std::cout << (weighted_bins ? "weighted bins" : "unweighted bins")
                  << " (" << atoms.size() << " atoms, " << waters.size() << " waters)" << std::endl;

        auto [cpu_result, cpu_ms] = time([&] () {
            hist::distance_calculator::SimpleCPU<weighted_bins, vbw> calculator;
            calculator.enqueue_calculate_self(atoms);
            calculator.enqueue_calculate_cross(atoms, waters);
            return calculator.run();
        });

        double pairs = 0.5*atoms.size()*(atoms.size() - 1) + static_cast<double>(atoms.size())*waters.size();
        std::cout << "    " << std::left << std::setw(7) << "cpu" << std::right
                  << std::fixed << std::setprecision(1) << std::setw(8) << cpu_ms << " ms"
                  << "                        (" << std::setprecision(1) << 1e-6*pairs/cpu_ms
                  << " Gpair/s over " << settings::general::threads << " threads)" << std::endl;

        #ifdef AUSAXS_GPU_SYCL
            if (run_sycl) {
                run_kernel<hist::distance_calculator::SYCLSimple<weighted_bins, vbw>, weighted_bins, vbw>(
                    "sycl", atoms, waters, cpu_result, pairs, cpu_ms
                );
            }
        #endif
        #ifdef AUSAXS_GPU
            if (run_webgpu) {
                run_kernel<hist::distance_calculator::WebGPUSimple<weighted_bins, vbw>, weighted_bins, vbw>(
                    "webgpu", atoms, waters, cpu_result, pairs, cpu_ms
                );
            }
        #endif
    }

    // the same comparison one level up: the scattering profile as the histogram manager produces it,
    // which is what settings::general::gpu switches for the rest of the library
    void run_manager(const data::Molecule& molecule) {
        auto profile = [&molecule] (bool gpu, settings::general::GPUBackend backend) {
            settings::general::gpu = gpu;
            settings::general::gpu_backend = backend;
            auto result = time([&molecule] () {return molecule.get_total_histogram()->debye_transform();});
            settings::general::gpu = false;
            return result;
        };

        auto report = [] (std::string_view name, const auto& cpu, const auto& gpu, double cpu_ms, double ms) {
            double max_relative = 0;
            int max_relative_index = -1;
            for (int i = 0; i < static_cast<int>(cpu.size()); ++i) {
                double c = cpu.get_counts()[i];
                double g = gpu.get_counts()[i];
                if (c == 0) {continue;}
                if (double relative = std::abs(c - g)/std::abs(c); relative > max_relative) {
                    max_relative = relative;
                    max_relative_index = i;
                }
            }
            std::cout << "    " << std::left << std::setw(7) << name << std::right
                      << std::fixed << std::setprecision(1) << std::setw(8) << ms << " ms"
                      << "   speedup vs cpu: " << std::setprecision(2) << (ms == 0 ? 0 : cpu_ms/ms) << "x"
                      << "   worst q-point " << max_relative_index << ": "
                      << std::scientific << std::setprecision(3) << max_relative << " relative" << std::endl;
        };

        std::cout << "scattering profile through the histogram manager" << std::endl;
        auto [cpu_intensity, cpu_ms] = profile(false, settings::general::GPUBackend::SYCL);
        std::cout << "    " << std::left << std::setw(7) << "cpu" << std::right
                  << std::fixed << std::setprecision(1) << std::setw(8) << cpu_ms << " ms" << std::endl;

        #ifdef AUSAXS_GPU_SYCL
        if (run_sycl) {
            auto [intensity, ms] = profile(true, settings::general::GPUBackend::SYCL);
            report("sycl", cpu_intensity, intensity, cpu_ms, ms);
        }
        #endif
        #ifdef AUSAXS_GPU
        if (run_webgpu) {
            auto [intensity, ms] = profile(true, settings::general::GPUBackend::WEBGPU);
            report("webgpu", cpu_intensity, intensity, cpu_ms, ms);
        }
        #endif
    }
}

int main(int argc, char const *argv[]) {
    std::string file = argc > 1 ? argv[1] : "tests/files/2epe.pdb";
    settings::general::verbose = false;

    // the thread pool is created on first use, so the count can only be set here, before anything runs
    if (argc > 2) {settings::general::threads = std::stoi(argv[2]);}
    if (argc > 3) {repeats = std::stoi(argv[3]);}
    if (argc > 4) {
        std::string backends = argv[4];
        run_sycl = backends.find("sycl") != std::string::npos;
        run_webgpu = backends.find("webgpu") != std::string::npos;
    }

    data::Molecule molecule(file);
    {   // the hydration shell is regenerated before every histogram evaluation in a rigidbody run,
        // so its cost bounds what overlapping it with the GPU calculation could save
        auto start = std::chrono::steady_clock::now();
        molecule.generate_new_hydration();
        auto end = std::chrono::steady_clock::now();
        std::cout << "hydration shell: " << std::fixed << std::setprecision(1)
                  << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;
    }
    std::cout << "best of " << repeats << " runs, " << settings::axes::bin_count << " bins" << std::endl;
    #ifdef AUSAXS_GPU_SYCL
        if (run_sycl) {std::cout << "sycl device: " << gpu::sycl_backend::device_name() << std::endl;}
    #endif

    run<false>(molecule);
    run<true>(molecule);
    run_manager(molecule);
    return 0;
}
