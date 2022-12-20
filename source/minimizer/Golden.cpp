#include <Symbols.h>
#include <mini/Golden.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

using namespace mini;

Golden::Golden(double(&func)(std::vector<double>)) : Minimizer(func) {}

Golden::Golden(std::function<double(std::vector<double>)> func) : Minimizer(func) {}

Golden::Golden(double(&func)(std::vector<double>), const Parameter& param) : Minimizer(func) {
    add_parameter(param);
}

Golden::Golden(std::function<double(std::vector<double>)> func, const Parameter& param) : Minimizer(func) {
    add_parameter(param);
}

Limit Golden::search(Limit bounds) const {
    // Code adapted from the python implementation from Wikipedia: https://en.wikipedia.org/wiki/Golden-section_search 
    double a = bounds.min, b = bounds.max;
    double temp = a + b;
    
    // sort such that a < b
    a = std::min(a, b);
    b = temp - a;

    double diff = b - a;
    if (diff < tol) [[unlikely]] {
        return Limit(a, b);
    }

    // expected number of steps to reach tolerance
    unsigned int n = std::ceil(std::log(tol/diff)/std::log(invphi));

    double c = a + invphi2*diff;
    double d = a + invphi*diff;
    double fc = function({c});
    double fd = function({d});

    for (unsigned int k = 0; k < n-1; k++) {
        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            diff = invphi*diff;
            c = a + invphi2*diff;
            fc = function({c});
        } else {
            a = c;
            c = d; 
            fc = fd;
            diff = invphi*diff;
            d = a + invphi*diff;
            d = a + invphi*diff;
            fd = function({d});
        }
    }

    if (fc < fd) {
        return Limit(a, d);
    } else {
        return Limit(c, b);
    }
}

Result Golden::minimize_override() {
    Limit optimal_interval = search(parameters[0].bounds.value());
    FittedParameter p(parameters[0], optimal_interval.center(), optimal_interval-optimal_interval.center());
    return Result(p, function({p.value}), fevals);
}

void Golden::add_parameter(const Parameter& param) {
    if (!param.has_bounds()) {throw except::invalid_argument("Golden::add_parameter: The parameter must be supplied with limits for this minimizer.");}
    if (!parameters.empty()) {throw except::invalid_operation("Golden::add_parameter: This minimizer only supports 1D problems.");}
    if (param.has_guess()) {debug_print("Warning in Golden::add_parameter: Guess value will be ignored.");}
    parameters.push_back(param);
}