#include <minimizer/Golden.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

using namespace mini;

Golden::Golden(double(&func)(double*), std::string par, Limit bounds) : Minimizer(func, 1) {
    add_parameter(par, bounds);
}

Golden::Golden(std::function<double(double*)> func, std::string par, Limit bounds) : Minimizer(func, 1) {
    add_parameter(par, bounds);
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
    unsigned int n = std::ceil(std::log(tol/diff)/std::log(invphi));

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

Dataset Golden::landscape(unsigned int evals) const {
    std::vector<double> x(evals), y(evals);

    double step = bounds.span()/evals;
    for (unsigned int i = 0; i < evals; i++) {
        double val = bounds.min + i*step;
        x[i] = val;
        y[i] = function(&val);
    }
    return Dataset(x, y);
}

Dataset Golden::get_evaluated_points() const {
    if (evaluations.empty()) {throw except::bad_order("Error in Golden::get_evaluated_points: Cannot get evaluated points before a minimization call has been made.");}

    unsigned int N = evaluations.size();
    std::vector<double> x(N), y(N);
    for (unsigned int i = 0; i < N; i++) {
        x[i] = evaluations[i].vals[0];
        y[i] = evaluations[i].fval;
    }
    return Dataset(x, y, "x", "f(x)");
}

Minimizer::Result Golden::minimize() const {
    Limit optimal_interval = search(bounds.min, bounds.max);
    double opt = optimal_interval.center();
    double err = std::max(opt-optimal_interval.min, optimal_interval.max-opt); // use largest error
    double opt_val = function(&opt); // function value at minimum

    return Minimizer::Result({{param_names[0], opt}}, {{param_names[0], err}}, opt_val);
}

void Golden::add_parameter(std::string, double) {
    throw except::disabled("Error in Golden::add_parameter: The parameter must be supplied with limits for this minimizer.");
}

void Golden::add_parameter(std::string par, double, Limit bounds) {
    if (!param_names.empty()) {throw except::invalid_operation("Error in Golden::add_parameter: This minimizer only supports 1D problems.");}

    utility::print_warning("Warning in Golden::add_parameter: Guess value will be ignored.");
    add_parameter(par, bounds);
}

void Golden::add_parameter(std::string par, Limit bounds) {
    if (!param_names.empty()) {throw except::invalid_operation("Error in Golden::add_parameter: This minimizer only supports 1D problems.");}

    param_names.push_back(par);
    this->bounds = bounds;
}