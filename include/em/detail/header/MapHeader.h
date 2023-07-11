#pragma once

#include <em/detail/header/data/DummyData.h>
#include <utility/UtilityFwd.h>

#include <memory>
#include <string>
#include <iosfwd>

namespace em::detail::header {
    struct HeaderData;
    class MapHeader {
        public:
            MapHeader(std::unique_ptr<HeaderData> data);

            virtual ~MapHeader();

            /**
             * @brief Create a string representation of this header.
             */
            virtual std::string to_string() const noexcept = 0;

            /**
             * @brief Get the data type for the map data. 
             */
            virtual em::detail::header::DataType get_data_type() const = 0;

            /**
             * @brief Get the size of the header.
             */
            virtual unsigned int get_header_size() const = 0;

            /**
             * @brief Get the axes of this map.
             */
            virtual Axis3D get_axes() const noexcept = 0;

            /**
             * @brief Get the index ordering of the data.
             * 
             * @return [x, y, z] where x, y, and z are the indices of the axes in the order they appear in the map.
             */
            virtual std::tuple<unsigned int, unsigned int, unsigned int> get_axis_order() const noexcept = 0;

            /**
             * @brief Get the byte size of each voxel.
             */
            unsigned int get_byte_size() const;

            /**
             * @brief Get the header data.
             */
            HeaderData* get_data() const noexcept;

            /**
             * @brief Set the header data.
             */
            void set_data(std::unique_ptr<HeaderData> data);

        private: 
            std::unique_ptr<HeaderData> data = nullptr;
    };
    std::ostream& operator<<(std::ostream& os, const MapHeader& h);
}
