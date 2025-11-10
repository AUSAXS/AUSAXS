#pragma once

#include <form_factor/xray/FormFactor.h>

namespace ausaxs::form_factor {
    class NormalizedFormFactor : public FormFactor {
        NormalizedFormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : FormFactor(a, b, c) {
            set_normalization(1);
        }
    };
}