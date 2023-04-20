#pragma once

#include <vector>
#include <string>

#include <utility/SmartOption.h>

namespace settings {
    namespace general {
        extern settings::detail::SmartOption<std::string> residue_folder;   // Download location for all ligand files. Must be constexpr. 
        extern settings::detail::SmartOption<bool> verbose;                 // Whether to print out extra information.
        extern settings::detail::SmartOption<unsigned int> threads;         // The number of threads to use for parallelization.
        extern settings::detail::SmartOption<std::string> output;           // The output directory.
        extern settings::detail::SmartOption<bool> keep_hydrogens;          // Whether to keep bound hydrogens when reading a structure.
        extern settings::detail::SmartOption<bool> supplementary_plots;     // Whether to generate supplementary plots when possible. 

        namespace detail {
            extern settings::detail::SmartOption<unsigned int> job_size;    // The number of atoms to process in each job.
        }
    }

    namespace detail {
        void parse_option(const std::string& name, const std::vector<std::string>& value);
        bool is_comment_char(char c);
    }
    
    /**
     * @brief Read the settings from a file.
     */
    void read(const std::string& path);

    /**
     * @brief Write the settings to a file.
     */
    void write(const std::string& path);

    /**
     * @brief Check if a settings file exists in the given directory, and read it if so.
     * @return True if a settings file was found and read.
     */
    bool discover(std::string path);
}