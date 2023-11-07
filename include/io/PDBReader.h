#pragma once

#include <io/Reader.h>
#include <data/detail/DataDetailFwd.h>

namespace io::detail {
    /**
     * @brief This class handles reading of input PDB format data files. 
     */
    class PDBReader : public Reader {
        public:
            /**
             * @brief Constructor.
             * @param file Path to the input PDB format data file. 
             */
            PDBReader(data::detail::AtomCollection* const file);

            ~PDBReader() override;

            /**
             * @brief Read a PDB format data file.
             * 
             * @param path Path to the input PDB format data file. 
             */
            void read(const io::File& path) override;

        private: 
            data::detail::AtomCollection* const file; // The File backing this Reader. 
    };
}