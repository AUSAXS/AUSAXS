#pragma once

#include <vector>

template<typename T> class Container2D;
namespace form_factor {
    class FormFactor;
    class PrecalculatedFormFactorProduct {
        public:
            PrecalculatedFormFactorProduct() = default;
            PrecalculatedFormFactorProduct(const FormFactor& ff1, const FormFactor& ff2, const std::vector<double>& q);

            double evaluate(unsigned int index) const;

            /**
             * @brief Generate a table of precalculated form factor products.
             *        Each entry is a precalculated form factor product for the default q axis. 
             *        The table is indexed by the form factor types, with the last column being the excluded volume. 
             *        It thus has the dimensions (#no_of_form_factors+1, #no_of_form_factors+1)
             *        The table is symmetric, so only the upper triangle is filled (i.e. i <= j)
             */
            static Container2D<PrecalculatedFormFactorProduct> generate_table();

        private:
            std::vector<double> precalculated_ff_q;
    };
}