#include <catch2/catch_test_macros.hpp>

#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <data/Molecule.h>
#include <settings/All.h>

#include <chrono>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>

using namespace ausaxs;

TEST_CASE("Automatic benchmark: A2M_native", "[.][benchmark_selected]") {
    settings::general::verbose = false;

    // Parameters from user request
    const int index = 14;
    const int iterations = 10;
    const int warmup = 3;

    // File for index 14 (A2M_native)
    const std::string pdb_candidates[] = {
        "data/A2M_native/A2M_native_fixed.pdb",
        "data/A2M_native/A2M_native.pdb"
    };

    std::string pdb_path;
    for (auto &p : pdb_candidates) {
        std::ifstream ifs(p);
        if (ifs.good()) { pdb_path = p; break; }
    }
    if (pdb_path.empty()) {
        INFO("Could not find A2M_native PDB under data/A2M_native; tried fixed and original names");
        FAIL("Missing input PDB for benchmark");
    }

    data::Molecule mol(pdb_path);
    INFO("Loaded molecule from " << pdb_path << " with " << mol.size_atom() << " atoms");

    struct Candidate { const char* name; settings::hist::HistogramManagerChoice choice; } candidates[] = {
        {"Simple", settings::hist::HistogramManagerChoice::HistogramManagerMT},
        {"Fraser", settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit},
        {"Grid",   settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid}
    };

    for (auto &c : candidates) {
        SECTION(c.name) {
            mol.set_histogram_manager(c.choice);

            // Warmup runs
            for (int w = 0; w < warmup; ++w) {
                (void)mol.get_histogram();
            }

            // Timed iterations
            using clock = std::chrono::high_resolution_clock;
            auto t0 = clock::now();
            for (int i = 0; i < iterations; ++i) {
                (void)mol.get_histogram();
            }
            auto t1 = clock::now();

            double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            double avg_ms = total_ms / iterations;

            std::ostringstream oss;
            oss << std::fixed << std::setprecision(3);
            oss << "Benchmark '" << c.name << "' on index=" << index
                << " (" << pdb_path << "): iterations=" << iterations
                << ", warmup=" << warmup << ", avg=" << avg_ms << " ms"
                << ", total=" << total_ms << " ms";

            INFO(oss.str());
            REQUIRE(avg_ms > 0.0);
        }
    }
}
