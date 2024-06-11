#pragma once

#include <dataset/SimpleDataset.h>

namespace fitter {
    struct FitPlots {
        SimpleDataset intensity;              // The full intensity line
        SimpleDataset intensity_interpolated; // The intensity line interpolated at the data points. 
        SimpleDataset data;                   // The data itself

        bool operator==(const FitPlots& other) const = default;
    };
}