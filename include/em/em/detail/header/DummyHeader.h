// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <em/detail/header/MapHeader.h>

namespace ausaxs::em::detail::header {
    /**
     * @brief Dummy header class. 
     */
    class DummyHeader : public MapHeader<DummyData> {
        public:
            DummyHeader();
            ~DummyHeader() override;

            std::string to_string() const override;

            unsigned int get_header_size() const override;

            em::detail::header::DataType get_data_type() const override;

            Axis3D get_axes() const noexcept override;
    };
}