#pragma once

#include <data/Protein.h>
#include <utility/Utility.h>
#include <utility/Exceptions.h>

#include <fstream>
#include <iomanip>

namespace io {
    /**
     * @brief XYZWriter.
     *        This class can write .xyz trajectory files.
     */
    class XYZWriter {
        public:
            XYZWriter(std::string path) : path(path) {
                utility::create_directory(path);
                file.open(path);
                if (!file.is_open()) {
                    throw except::io_error("io::XYZWriter: Could not open file " + path);
                }
            }

            ~XYZWriter() {
                file.close();
                utility::print_info("Trajectory written to " + path);
            }

            /**
             * @brief Write a frame to the file.
             */
            void write_frame(const Protein& protein) {
                static unsigned int frame = 0;
                auto atoms = protein.atoms();
                file << " " << atoms.size() << std::endl;
                file << " Frame " << frame++ << std::endl;
                for (const auto& atom : atoms) {
                    file << std::setw(10) << atom.serial << " " 
                         << std::setw(10) << std::setprecision(6) << atom.coords.x() << " " 
                         << std::setw(10) << std::setprecision(6) << atom.coords.y() << " " 
                         << std::setw(10) << std::setprecision(6) << atom.coords.z() << std::endl;
                }
            }

        private:
            std::ofstream file;
            std::string path;
    };
}