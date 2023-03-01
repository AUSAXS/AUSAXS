#pragma once

#include <crystal/miller/ReducedMillers.h>
#include <math/Vector3.h>

#include <math.h>
#include <iostream>

namespace crystal {
    class FibonacciMillers : public ReducedMillers {
        public:
            FibonacciMillers(unsigned int h, unsigned int k, unsigned int l);

            std::vector<Miller> generate() const override;

        private:
            int h, k, l;

            std::vector<Vector3<double>> generate_fibonacci_sphere(int n) const;
    };
}