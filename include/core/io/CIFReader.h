#pragma once

#include <io/Reader.h>
#include <data/detail/DataDetailFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::io::detail {
    /**
     * @brief This class handles reading of input PDB format data files. 
     */
    class CIFReader : public Reader {
        public:
            CIFReader(observer_ptr<data::detail::AtomCollection> const file);
            ~CIFReader() override;

            void read(const io::File& path) override;

        private: 
            observer_ptr<data::detail::AtomCollection> file; // The File backing this Reader. 
    };
}