// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <io/IOFwd.h>

#include <fstream>

namespace ausaxs::io::detail::xyz {
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