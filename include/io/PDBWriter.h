#pragma once

#include <io/Writer.h>
#include <data/detail/DataDetailFwd.h>

#include <vector>
#include <string>

namespace io::detail {
    /**
     * @brief This class handles writing a File object into a PDB format data file.
     */
    class PDBWriter : public Writer {
        public:
            /**
             * @brief Constructor. 
             * @param file Path to where the backing File object will be saved. 
             */
            PDBWriter(data::detail::AtomCollection* file);

            ~PDBWriter() override;

            /**
             * @brief Write the backing File to disk. If the size of the protein is greater than 100 000, multiple files will be written.
             * 
             * @param path The output path. 
             */
            void write(const io::File& path) override;

        private: 
            data::detail::AtomCollection* file; // The File backing this Reader. 

            /**
             * @brief Create a string representation of this File.
             * @return The string representation. Each index is a file to be written. 
             */
            std::vector<std::string> as_pdb() const;
    };
}