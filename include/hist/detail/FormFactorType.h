#pragma once

#include <string_view>

namespace hist::detail::ff {
    // The form factor type of an atom. This is intended to be used as an index for best performance.
    enum class form_factor_t {
        HYDROGEN,
        CARBON,
        NITROGEN,
        OXYGEN,
        OTHER,
        COUNT
    };

    /**
     * @brief Get the number of unique form factors.
     */
    constexpr unsigned int get_count();

    /**
     * @brief Get the form factor type of an element from the Atom class. 
     */
    form_factor_t get_type(std::string_view element);

    class EvalFormFactor {
        public:
            constexpr EvalFormFactor(form_factor_t ff1, form_factor_t ff2);

            double operator()(double q) const;

        private: 
            double factor = 1;
    };

    constexpr EvalFormFactor get_form_factor(int ff_type1, int ff_type2);

    struct FormFactor {
        constexpr FormFactor(std::array<double, 5> a, std::array<double, 5> b, double c) : a(a), b(b), c(c) {}

        double evaluate(double q) const {
            double sum = 0;
            for (unsigned int i = 0; i < 5; ++i) {
                sum += a[i]*std::exp(-b[i]*q*q);
            }
            return sum + c;
        }

        std::array<double, 5> a;
        std::array<double, 5> b;
        double c;
    };

    FormFactor C({0.489918, 0.262003, 0.196767, 0.049879, 0.003971}, {20.6593, 10.5109, 0.893666, 51.6512, 0.215195}, 0.000000);
    FormFactor O({0.489918, 0.262003, 0.196767, 0.049879, 0.003971}, {20.6593, 10.5109, 0.893666, 51.6512, 0.215195}, 0.000000);
}