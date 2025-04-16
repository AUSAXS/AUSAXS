#pragma once

#include <string>

namespace ausaxs::settings {
    struct general {
        static std::string residue_folder;  // Download location for all ligand files.
        static bool offline;                // Whether to use the offline mode. This will prevent any network requests, but may result in less accurate results.
        static bool verbose;                // Whether to print out extra information.
        static bool warnings;               // Whether to print out warnings.
        static unsigned int threads;        // The number of threads to use for parallelization.
        static std::string output;          // The output directory.
        static std::string cache;           // The cache directory.
        static bool keep_hydrogens;         // Whether to keep bound hydrogens when reading a structure.
        static bool generate_plots;         // Whether to generate plots when possible.
        static bool supplementary_plots;    // Whether to generate supplementary plots when possible.

        struct detail {
            static unsigned int job_size;   // The number of atoms to process in each job.
        };

        enum class QUnit : char {
            A       = (1 << 1),          // Ångström  (system-default, may be converted)
            USER_A  = (1 << 1)+(1 << 3), // Ångström  (user-specified, cannot be converted)
            NM      = (1 << 2),          // Nanometer (system-default, may be converted)
            USER_NM = (1 << 2)+(1 << 3), // Nanometer (user-specified, cannot be converted)
        };
        static QUnit input_q_unit; // Unit of q values in the input file.

        struct helper {
            static bool is_angstroms(QUnit u);    // Whether the input unit is Ångström.
            static bool is_nanometers(QUnit u);   // Whether the input unit is Nanometer.
            static bool is_user_defined(QUnit u); // Whether the input unit is user-defined.
        };
    };
}