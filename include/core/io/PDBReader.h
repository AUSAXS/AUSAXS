#pragma once

#include <io/Reader.h>
#include <data/detail/DataDetailFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::io::detail {
    /**
     * @brief This class handles reading of input PDB format data files. 
     */
    class PDBReader : public Reader {
        public:
            PDBReader(observer_ptr<data::detail::AtomCollection> const file);
            ~PDBReader() override;

            void read(const io::File& path) override;

        private: 
            observer_ptr<data::detail::AtomCollection> const file; // The File backing this Reader. 
    };
}