#pragma once

#include <string>
#include <vector>

namespace ausaxs::md {
    class TopologyFile {
        public:
            TopologyFile(const std::string& path) : path(path) {}

            /**
             * @brief Fix relative includes in the topology file.
             */
            void fix_relative_includes();

            /**
             * @brief Discover all includes in the topology file.
             */
            void discover_includes();

            /**
             * @brief Replace the position restraints with the given restraints.
             */
            void replace_restraints(const std::vector<std::string>& restraints);

        private:
            struct Include {
                unsigned int line;
                std::string path;
            };

            std::string path;
            std::vector<Include> includes;
    };
}