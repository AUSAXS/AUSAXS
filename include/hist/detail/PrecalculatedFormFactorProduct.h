#pragma once

#include <vector>

template<typename T> class Container2D;
namespace hist::detail {
    class FormFactor;
    class PrecalculatedFormFactorProduct {
        public:
            PrecalculatedFormFactorProduct() = default;
            PrecalculatedFormFactorProduct(const FormFactor& ff1, const FormFactor& ff2, const std::vector<double>& q);

            double evaluate(unsigned int index) const;

            static Container2D<PrecalculatedFormFactorProduct> generate_table();

        private:
            std::vector<double> precalculated_ff_q;
    };
}