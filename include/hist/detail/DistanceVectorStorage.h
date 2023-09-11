#pragma once

#include <vector>

namespace hist::detail {
    class DistanceVectorsBase {
        public: 
            DistanceVectorsBase() = default;
            ~DistanceVectorsBase() = default;

            virtual std::vector<double> debye_transform() const = 0;
    };

    class DistanceVectors : public DistanceVectorsBase {

    };

    class DistanceVectorsCollection : public DistanceVectorsBase {

    };
}