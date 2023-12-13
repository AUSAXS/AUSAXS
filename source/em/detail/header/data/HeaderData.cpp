#include <em/detail/header/data/HeaderData.h>

#include <cstdint>

namespace em::detail::header {
    std::unordered_map<em::detail::header::DataType, unsigned int> byte_sizes = {
        {em::detail::header::DataType::int8,    static_cast<unsigned int>(sizeof(int8_t))},
        {em::detail::header::DataType::int16,   static_cast<unsigned int>(sizeof(int16_t))},
        {em::detail::header::DataType::uint8,   static_cast<unsigned int>(sizeof(uint8_t))},
        {em::detail::header::DataType::uint16,  static_cast<unsigned int>(sizeof(uint16_t))},
        {em::detail::header::DataType::float16, 2u},
        {em::detail::header::DataType::float32, static_cast<unsigned int>(sizeof(float))},
    };
}