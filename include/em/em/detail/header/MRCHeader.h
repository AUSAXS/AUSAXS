// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <em/detail/header/MapHeader.h>
#include <em/detail/header/data/MRCData.h>
#include <io/IOFwd.h>

namespace ausaxs::em::detail::header {
    /**
     * @brief Wrapper class for MRCData.
     */
    class MRCHeader : public MapHeader<MRCData> {
        public:
            MRCHeader();

            MRCHeader(MRCData&& data);

            ~MRCHeader() override;

            /**
             * @brief Create a string representation of this object.
             */
            std::string to_string() const override;

            /**
             * @brief Get the size of the header.
             */
            unsigned int get_header_size() const override;

            /**
             * @brief Get the data type for the map data. 
             * 
             * @return std::string representation of the data type. 
             */
            em::detail::header::DataType get_data_type() const override;

            /**
             * @brief Get the axes of this map.
             */
            Axis3D get_axes() const noexcept override;

            /**
             * @brief Get the index ordering of the data.
             * 
             * @return [x, y, z] where x, y, and z are the indices of the axes in the order they appear in the map.
             */
            std::tuple<unsigned int, unsigned int, unsigned int> get_axis_order() const noexcept override;

            /**
             * @brief Rotate the map contents. This does not affect the operation of this program.
             *        The arguments must be some permutation of {1, 2, 3}.
             */
            void rotate(int x, int y, int z) noexcept;

            /**
             * @brief Check if a given path is a MRC format file.
             */
            static bool is_mrc(const io::ExistingFile& file);
        
        private:
            MRCData& cast_data() const noexcept;
    };
}