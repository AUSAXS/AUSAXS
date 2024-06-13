#include <em/detail/header/MapHeader.h>
#include <em/detail/header/data/RECData.h>
#include <io/IOFwd.h>

namespace em::detail::header {
    /**
     * @brief Wrapper class for RECData.
     */
    class RECHeader : public MapHeader<RECData> {
        public:
            RECHeader();
            ~RECHeader() override;

            /**
             * @brief Create a string representation of this object.
             */
            std::string to_string() const override;

            /**
             * @brief Get the data type for the map data. 
             * 
             * @return std::string representation of the data type. 
             */
            em::detail::header::DataType get_data_type() const override;

            /**
             * @brief Get the size of the header.
             */
            unsigned int get_header_size() const override;

            /**
             * @brief Get the index ordering of the data.
             * 
             * @return [x, y, z] where x, y, and z are the indices of the axes in the order they appear in the map.
             */
            std::tuple<unsigned int, unsigned int, unsigned int> get_axis_order() const noexcept override;

            /**
             * @brief Get the axes of this map.
             */
            Axis3D get_axes() const noexcept override;

            /**
             * @brief Check if a given path is a REC format file.
             */
            static bool is_rec(const io::ExistingFile& file);

        private:
            bool flags_enabled() const noexcept;

            bool flag_signed_bytes() const noexcept;
            bool flag_four_bit_vals() const noexcept;

            RECData& cast_data() const noexcept;
    };
}