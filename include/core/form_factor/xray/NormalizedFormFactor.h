#pragma once

#include <form_factor/xray/FormFactor.h>

namespace ausaxs::form_factor {
    class NormalizedFormFactor : public UnnormalizedFormFactor {
        NormalizedFormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : UnnormalizedFormFactor(a, b, c) {
            set_normalization(1);
        }
    };
}