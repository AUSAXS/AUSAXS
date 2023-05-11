#include <Symbols.h>
#include <mini/Scan.h>
#include <mini/Golden.h>
#include <mini/detail/Parameter.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

using namespace mini;

Scan::Scan(double(&func)(std::vector<double>), unsigned int evals) : Minimizer(func) {
    set_max_evals(evals);
}

Scan::Scan(std::function<double(std::vector<double>)> func, unsigned int evals) : Minimizer(func) {
    set_max_evals(evals);
}

Scan::Scan(double(&func)(std::vector<double>), const Parameter& param, unsigned int evals) : Minimizer(func) {
    set_max_evals(evals);
    add_parameter(param);
}

Scan::Scan(std::function<double(std::vector<double>)> func, const Parameter& param, unsigned int evals) : Minimizer(func) {
    set_max_evals(evals);
    add_parameter(param);
}

mini::Landscape Scan::landscape(unsigned int evals) {
    // check if the minimizer has already been called
    if (!evaluations.evals.empty()) {
        // if so, we can just reuse its result
        return evaluations;
    }

    if (parameters.size() == 1) {
        const Limit& bounds = parameters[0].bounds.value();
        for (double val = bounds.min; val < bounds.max; val += bounds.span()/evals) {
            function({val});
        }
        return get_evaluated_points();
    } 
    
    else if (parameters.size() == 2) {
        throw except::unexpected("Scan::landscape: Not implemented.");
    } 
    
    else {
        throw except::unexpected("Scan::landscape: Not implemented.");
    }
}

void Scan::add_parameter(const Parameter& param) {
    if (!param.has_bounds()) {throw except::invalid_argument("Scan::add_parameter: The parameter must be supplied with limits for this minimizer.");}
    if (!parameters.empty()) {throw except::invalid_operation("Scan::add_parameter: This minimizer only supports 1D problems.");}
    parameters.push_back(param);
}

Result Scan::minimize_override() {
    SimpleDataset data = landscape(max_evals).as_dataset();
    auto min = data.find_minimum();

    // find local minimum
    auto width = data.span_x().span()/data.size(); // find width of each step
    auto prev_bounds = parameters[0].bounds;
    parameters[0].bounds = Limit(std::max(min.x - width, prev_bounds->min), std::min(min.x + width, prev_bounds->max)); // update bounds
    parameters[0].guess = {}; // remove guess to avoid warning

    // record_evaluations(false);
    mini::Golden golden(function, parameters[0]);
    return golden.minimize();
}
