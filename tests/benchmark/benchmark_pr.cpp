#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
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

    // Generate a random molecule with N atoms distributed over M bodies for benchmarking
    data::Molecule generate_random_molecule(size_t n_atoms, size_t n_bodies, unsigned int seed = 42) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<> pos_dist(-100.0, 100.0);
        std::uniform_int_distribution<> element_dist(0, 4);
        std::uniform_int_distribution<> body_dist(0, n_bodies - 1);

        const std::array<form_factor::form_factor_t, 5> ff_types = {
            form_factor::form_factor_t::C,
            form_factor::form_factor_t::N,
            form_factor::form_factor_t::O,
            form_factor::form_factor_t::H,
            form_factor::form_factor_t::S
        };

        std::vector<std::vector<data::AtomFF>> bodies(n_bodies);
        for (size_t i = 0; i < n_atoms; ++i) {
            size_t body_idx = body_dist(gen);
            bodies[body_idx].emplace_back(
                Vector3<double>(pos_dist(gen), pos_dist(gen), pos_dist(gen)),
                ff_types[element_dist(gen)]
            );
        }

        std::vector<data::Body> body_vec;
        body_vec.reserve(n_bodies);
        for (auto& atom_vec : bodies) {
            body_vec.emplace_back(std::move(atom_vec));
        }
        return data::Molecule(std::move(body_vec));
    }
}

TEST_CASE("Distance calculation benchmark: large systems") {
    settings::general::verbose = false;
    auto atoms = GENERATE(1000, 2000, 3000, 4000, 5000, 10000, 15000, 20000, 25000, 50000, 75000, 100000);

    SECTION(std::to_string(atoms) + " atoms") {
        auto mol = generate_random_molecule(atoms);
        
        BENCHMARK_ADVANCED("Simple") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
            meter.measure([&] { return mol.get_histogram(); });
        };
        
        BENCHMARK_ADVANCED("Fraser") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit);
            meter.measure([&] { return mol.get_histogram(); });
        };

        BENCHMARK_ADVANCED("Grid") (Catch::Benchmark::Chronometer meter) {
            mol.set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid);
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
