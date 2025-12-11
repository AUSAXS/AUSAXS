#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>

#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogram.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/atoms/AtomFF.h>
#include <form_factor/FormFactorType.h>
#include <settings/All.h>

#include <random>
#include <vector>

using namespace ausaxs;

namespace {
    // Generate a random molecule with N atoms for benchmarking
    data::Molecule generate_random_molecule(size_t n_atoms, unsigned int seed = 42) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<> pos_dist(-100.0, 100.0);
        std::uniform_int_distribution<> element_dist(0, 4);
        
        const std::array<form_factor::form_factor_t, 5> ff_types = {
            form_factor::form_factor_t::C,
            form_factor::form_factor_t::N,
            form_factor::form_factor_t::O,
            form_factor::form_factor_t::H,
            form_factor::form_factor_t::S
        };
        
        std::vector<data::AtomFF> atoms;
        atoms.reserve(n_atoms);
        
        for (size_t i = 0; i < n_atoms; ++i) {
            atoms.emplace_back(
                Vector3<double>(pos_dist(gen), pos_dist(gen), pos_dist(gen)),
                ff_types[element_dist(gen)]
            );
        }
        
        return data::Molecule({data::Body(std::move(atoms))});
    }
}

TEST_CASE("Distance calculation benchmark: large systems") {
    settings::general::verbose = false;
    
    SECTION("20k atoms") {
        auto mol = generate_random_molecule(20000);
        
        BENCHMARK_ADVANCED("HistogramManagerMT (simple)")(Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
            meter.measure([&] { return mol.get_histogram(); });
        };
        
        BENCHMARK_ADVANCED("HistogramManagerMTFFAvg")(Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg);
            meter.measure([&] { return mol.get_histogram(); });
        };
        
        BENCHMARK_ADVANCED("HistogramManagerMTFFExplicit")(Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit);
            meter.measure([&] { return mol.get_histogram(); });
        };
    }
    
    SECTION("50k atoms") {
        auto mol = generate_random_molecule(50000);
        
        BENCHMARK_ADVANCED("HistogramManagerMT (simple)")(Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
            meter.measure([&] { return mol.get_histogram(); });
        };
        
        BENCHMARK_ADVANCED("HistogramManagerMTFFAvg")(Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg);
            meter.measure([&] { return mol.get_histogram(); });
        };
        
        BENCHMARK_ADVANCED("HistogramManagerMTFFExplicit")(Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit);
            meter.measure([&] { return mol.get_histogram(); });
        };
    }
    
    SECTION("100k atoms") {
        auto mol = generate_random_molecule(100000);
        
        BENCHMARK_ADVANCED("HistogramManagerMT (simple)")(Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
            meter.measure([&] { return mol.get_histogram(); });
        };
        
        BENCHMARK_ADVANCED("HistogramManagerMTFFAvg")(Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg);
            meter.measure([&] { return mol.get_histogram(); });
        };
        
        BENCHMARK_ADVANCED("HistogramManagerMTFFExplicit")(Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit);
            meter.measure([&] { return mol.get_histogram(); });
        };
    }
}

TEST_CASE("Distance calculation benchmark: real PDB file", "[.][pdb]") {
    // This test is tagged [.] so it's skipped by default (requires actual PDB file)
    // Run with: ./benchmark_tests "[pdb]"
    settings::general::verbose = false;
    
    // Try to load a real structure if available
    std::string pdb_path = "tests/files/2epe.pdb";
    
    data::Molecule mol(pdb_path);
    INFO("Loaded molecule with " << mol.size_atom() << " atoms");
    
    BENCHMARK_ADVANCED("Real PDB - FFAvg")(Catch::Benchmark::Chronometer meter) {
        mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg);
        meter.measure([&] { return mol.get_histogram(); });
    };
    
    BENCHMARK_ADVANCED("Real PDB - FFExplicit")(Catch::Benchmark::Chronometer meter) {
        mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit);
        meter.measure([&] { return mol.get_histogram(); });
    };
}
