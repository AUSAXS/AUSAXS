#pragma once

#include <string>

namespace ausaxs::settings::general {
    extern std::string residue_folder;  // Download location for all ligand files.
    extern bool offline;                // Whether to use the offline mode. This will prevent any network requests, but may result in less accurate results.
    extern bool verbose;                // Whether to print out extra information.
    extern bool warnings;               // Whether to print out warnings.
    extern unsigned int threads;        // The number of threads to use for parallelization.
    extern std::string output;          // The output directory.
    extern std::string cache;           // The cache directory.
    extern bool keep_hydrogens;         // Whether to keep bound hydrogens when reading a structure.
    extern bool generate_plots;         // Whether to generate plots when possible.
    extern bool supplementary_plots;    // Whether to generate supplementary plots when possible.

    namespace detail {
        extern unsigned int job_size;   // The number of atoms to process in each job.
    }
}

namespace ausaxs::settings::general {
    enum class QUnit : char {
        A       = (1 << 1),          // Ångström  (system-default, may be converted)
        USER_A  = (1 << 1)+(1 << 3), // Ångström  (user-specified, cannot be converted)
        NM      = (1 << 2),          // Nanometer (system-default, may be converted)
        USER_NM = (1 << 2)+(1 << 3), // Nanometer (user-specified, cannot be converted)
    };
    extern QUnit input_q_unit; // Unit of q values in the input file.

    namespace helper {
        bool is_angstroms(QUnit u);    // Whether the input unit is Ångström.
        bool is_nanometers(QUnit u);   // Whether the input unit is Nanometer.
        bool is_user_defined(QUnit u); // Whether the input unit is user-defined.
    }
}