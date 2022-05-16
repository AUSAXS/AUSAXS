#include <minimizer/Golden.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

using namespace mini;

Golden::Golden(std::string par, Limit limits) {
    add_parameter(par, limits);
}

Limit Golden::search(double a, double b) const {
    // Code adapted from the python implementation from Wikipedia: https://en.wikipedia.org/wiki/Golden-section_search 
    double temp = a + b;
    
    // sort such that a < b
    a = std::min(a, b);
    b = temp - a;

    double diff = b - a;
    if (__builtin_expect(diff < tol, false)) {
        return Limit(a, b);
    }

    // expected number of steps to reach tolerance
    int n = std::ceil(std::log(tol/diff)/std::log(invphi));

    double c = a + invphi2*diff;
    double d = a + invphi*diff;
    double fc = function(&c);
    double fd = function(&d);

    for (unsigned int k = 0; k < n-1; k++) {
        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            diff = invphi*diff;
            c = a + invphi2*diff;
            fc = function(&c);
        } else {
            a = c;
            c = d; 
            fc = fd;
            diff = invphi*diff;
            d = a + invphi*diff;
            d = a + invphi*diff;
            fd = function(&d);
        }
    }

    if (fc < fd) {
        return Limit(a, d);
    } else {
        return Limit(c, b);
    }
}

Minimizer::Result Golden::minimize() const {
    Limit optimal_interval = search(limit.min, limit.max);
    double opt = optimal_interval.center();
    double err = std::max(opt-optimal_interval.min, optimal_interval.max-opt); // use largest error
    double opt_val = function(&opt); // function value at minimum

    return Minimizer::Result({{param_names[0], opt}}, {{param_names[0], err}}, opt_val);
}

void Golden::add_parameter(std::string par, double guess) {
    throw except::disabled("Error in Golden::add_parameter: The parameter must be supplied with limits for this minimizer.");
}

void Golden::add_parameter(std::string par, double guess, Limit limits) {
    if (!param_names.empty()) {throw except::invalid_operation("Error in Golden::add_parameter: This minimizer only supports 1D problems.");}

    utility::print_warning("Warning in Golden::add_parameter: Guess value will be ignored.");
    add_parameter(par, limits);
}

void Golden::add_parameter(std::string par, Limit limits) {
    if (!param_names.empty()) {throw except::invalid_operation("Error in Golden::add_parameter: This minimizer only supports 1D problems.");}

    param_names.push_back(par);
    this->limit = limits;
}