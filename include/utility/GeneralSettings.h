#pragma once

#include <string>
#include <vector>

namespace settings {
    namespace general {
        extern const char* const residue_folder;    // Download location for all ligand files. Must be constexpr.
        extern bool verbose;                        // Whether to print out extra information.
        extern unsigned int threads;                // The number of threads to use for parallelization.
        extern std::string output;                  // The output directory.
        extern bool keep_hydrogens;                 // Whether to keep bound hydrogens when reading a structure.
        extern bool supplementary_plots;            // Whether to generate supplementary plots when possible.

        namespace detail {
            extern unsigned int job_size;           // The number of atoms to process in each job.
        }
    }
}