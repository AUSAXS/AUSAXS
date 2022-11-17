#include <mini/Scan.h>
#include <mini/Golden.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

using namespace mini;

Scan::Scan(double(&func)(std::vector<double>), unsigned int evals) : Minimizer(func), bins(evals) {}

Scan::Scan(std::function<double(std::vector<double>)> func, unsigned int evals) : Minimizer(func), bins(evals) {}

Scan::Scan(double(&func)(std::vector<double>), const Parameter& param, unsigned int evals) : Minimizer(func), bins(evals) {
    add_parameter(param);
}

Scan::Scan(std::function<double(std::vector<double>)> func, const Parameter& param, unsigned int evals) : Minimizer(func), bins(evals) {
    add_parameter(param);
}

void Scan::set_evals(unsigned int evals) noexcept {
    bins = evals;
}

Dataset2D Scan::landscape(unsigned int evals) {
    // check if the minimizer has already been called
    if (!evaluations.empty()) {
        // if so, we can just reuse its result
        Dataset2D data;
        std::for_each(evaluations.begin(), evaluations.end(), [&data] (const Evaluation& eval) {data.push_back(eval.vals[0], eval.fval);});
        return data;
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

void Scan::looper(std::vector<double>&, unsigned int) const {
    throw except::unexpected("Scan::looper: Not implemented.");
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

Dataset2D Scan::get_evaluated_points() const {
    if (evaluations.empty()) {throw except::bad_order("Scan::get_evaluated_points: Cannot get evaluated points before a minimization call has been made.");}

    unsigned int N = evaluations.size();
    std::vector<double> x(N), y(N);
    for (unsigned int i = 0; i < N; i++) {
        x[i] = evaluations[i].vals[0];
        y[i] = evaluations[i].fval;
    }
    return Dataset2D(x, y, "x", "f(x)");
}

void Scan::add_parameter(const Parameter& param) {
    if (!param.has_bounds()) {throw except::invalid_argument("Scan::add_parameter: The parameter must be supplied with limits for this minimizer.");}
    if (!parameters.empty()) {throw except::invalid_operation("Scan::add_parameter: This minimizer only supports 1D problems.");}
    if (param.has_guess()) {utility::print_warning("Warning in Scan::add_parameter: Guess value will be ignored.");}
    parameters.push_back(param);
}

Result Scan::minimize_override() {
    Dataset2D data = landscape(bins);
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
