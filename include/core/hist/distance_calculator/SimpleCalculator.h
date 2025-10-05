#pragma once

#include <hist/distance_calculator/SimpleKernel.h>
#include <hist/distance_calculator/SimpleCPU.h>
#include <hist/distance_calculator/SimpleGPU.h>

#include <memory>

namespace ausaxs::hist::distance_calculator {
    template<bool weighted_bins>
    struct SimpleCalculator {
        static std::unique_ptr<SimpleKernel<weighted_bins>> create() {
            #if defined(__ACPP__)
                return std::make_unique<SimpleCPU<weighted_bins>>();
            #elif defined(__WEBGPU__)
                return std::make_unique<WebGPUSimple<weighted_bins>>();
            #else
                return std::make_unique<SimpleCPU<weighted_bins>>();
            #endif
        }
    };
}