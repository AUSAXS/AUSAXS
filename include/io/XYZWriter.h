#pragma once

#include <data/DataFwd.h>

#include <fstream>

namespace io {class File;}
namespace io {
    /**
     * @brief This class can write .xyz trajectory files.
     */
    class XYZWriter {
        public:
            XYZWriter(const io::File& path);

            ~XYZWriter();

            /**
             * @brief Write a frame to the file.
             */
            void write_frame(const data::Molecule* protein);

        private:
            std::ofstream file;
            std::string path;
    };
}