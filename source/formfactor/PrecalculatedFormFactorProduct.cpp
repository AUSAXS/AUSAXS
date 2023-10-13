#include <formfactor/PrecalculatedFormFactorProduct.h>
#include <formfactor/FormFactor.h>
#include <settings/HistogramSettings.h>
#include <utility/Axis.h>
#include <utility/Container2D.h>

#include <cmath>

using namespace form_factor;

PrecalculatedFormFactorProduct::PrecalculatedFormFactorProduct(const FormFactor& ff1, const FormFactor& ff2, const std::vector<double>& q) {
    std::vector<double> res(q.size());
    for (unsigned int i = 0; i < q.size(); ++i) {
        res[i] = ff1.evaluate(q[i])*ff2.evaluate(q[i]);
    }
    precalculated_ff_q = std::move(res);
}

double PrecalculatedFormFactorProduct::evaluate(unsigned int index) const {
    return precalculated_ff_q[index];
}

Container2D<PrecalculatedFormFactorProduct> PrecalculatedFormFactorProduct::generate_table() {
    Container2D<PrecalculatedFormFactorProduct> table(FormFactor::get_count(), FormFactor::get_count());
    std::vector<double> q = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins).as_vector();
    for (unsigned int i = 0; i < FormFactor::get_count(); ++i) {
        for (unsigned int j = 0; j < i; ++j) {
            table.index(i, j) = std::move(
                PrecalculatedFormFactorProduct(
                    storage::get_form_factor(static_cast<form_factor_t>(i)), 
                    storage::get_form_factor(static_cast<form_factor_t>(j)), 
                    q
                )
            );
            table.index(j, i) = table.index(i, j);
        }
        table.index(i, i) = std::move(
            PrecalculatedFormFactorProduct(
                storage::get_form_factor(static_cast<form_factor_t>(i)), 
                storage::get_form_factor(static_cast<form_factor_t>(i)), 
                q
            )
        );
    }
    return table;
}
