#include <minimizer/Scan.h>
#include <utility/Exceptions.h>

using namespace mini;

Dataset Scan::landscape(unsigned int evals) const {
    if (parameters.size() == 1) {
        Limit& bounds = parameters[0].bounds.value();
        for (double val = bounds.min; val < bounds.max; val += bounds.span()/bins) {
            function(&val);
        }
        return get_evaluated_points();
    } 
    
    else if (parameters.size() == 2) {
        throw except::unexpected("Error in Scan::landscape: Not implemented.");
    }
}

void Scan::looper(std::vector<double>& p, unsigned int index) const {
    throw except::unexpected("Error in Scan::looper: Not implemented.");
    // Limit bounds = parameters[index].bounds.value();
    // for (double val = bounds.min; val < bounds.max; val += bounds.span()/bins) {
    //     p[index] = val;
    //     if (index < parameters.size()) {
    //         looper(p, index+1);
    //     } else {
    //         if (function(p.data()) < current_best) {
    //             // update best pars
    //         }
    //     }
    // }
}

Dataset Scan::get_evaluated_points() const {
    if (evaluations.empty()) {throw except::bad_order("Error in Scan::get_evaluated_points: Cannot get evaluated points before a minimization call has been made.");}

    unsigned int N = evaluations.size();
    std::vector<double> x(N), y(N);
    for (unsigned int i = 0; i < N; i++) {
        x[i] = evaluations[i].vals[0];
        y[i] = evaluations[i].fval;
    }
    return Dataset(x, y, "x", "f(x)");
}

void Scan::add_parameter(const Parameter& param) {
    if (!param.has_bounds()) {throw except::invalid_argument("Error in Scan::add_parameter: The parameter must be supplied with limits for this minimizer.");}
    if (!parameters.empty()) {throw except::invalid_operation("Error in Scan::add_parameter: This minimizer only supports 1D problems.");}
    if (param.has_guess()) {utility::print_warning("Warning in Scan::add_parameter: Guess value will be ignored.");}
    parameters.push_back(param);
}

Result Scan::minimize_override() {
    Dataset data = landscape(bins);
    auto min = data.find_minimum();
    FittedParameter p(parameters[0], min.x, 0);
    return Result(p, min.y, fevals);
}
