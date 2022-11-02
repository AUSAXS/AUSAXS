#include <minimizer/dlibMinimizer.h>

using namespace mini;

dlibMinimizer::dlibMinimizer(std::function<double(double)> function, Parameter param) {
    // this->function = function;
    // wrapper = [this](column_vector x) { return this->function(x(0));};
}

dlibMinimizer::dlibMinimizer(std::function<double(const double*)> function, Parameter param) {
    add_parameter(param);
    set_function(function);
    dlib_fwrapper = [this](column_vector x) {return this->function(&x(0));};
}

dlibMinimizer::~dlibMinimizer() = default;

Dataset2D dlibMinimizer::landscape(unsigned int evals) {
    throw except::unexpected("dlibMinimizer::landscape: Not implemented yet.");
}

Dataset2D dlibMinimizer::get_evaluated_points() const {
    if (evaluations.empty()) {throw except::bad_order("dlibMinimizer::get_evaluated_points: Cannot get evaluated points before a minimization call has been made.");}

    unsigned int N = evaluations.size();
    std::vector<double> x(N), y(N);
    for (unsigned int i = 0; i < N; i++) {
        x[i] = evaluations[i].vals[0];
        y[i] = evaluations[i].fval;
    }
    return Dataset2D(x, y, "x", "f(x)");    
}

Result dlibMinimizer::minimize_override() {
    if (!is_parameter_set()) {throw except::bad_order("dlibMinimizer::minimize: No parameters were supplied.");}
    if (!is_function_set()) {throw except::bad_order("dlibMinimizer::minimize: No function was set.");}

    // prepare guess value
    column_vector x;
    if (parameters[0].has_guess()) {
        x = {*parameters[0].guess};
    } else if (parameters[0].has_bounds()) {
        x = {parameters[0].bounds->center()};
    } else {
        throw except::invalid_argument("dlibMinimizer::minimize: Either a guess or bounds must be supplied.");
    }

    double fmin;
    if (parameters[0].has_bounds()) {
        fmin = dlib::find_min_box_constrained(
            dlib::bfgs_search_strategy(), 
            dlib::objective_delta_stop_strategy(1e-7), 
            dlib_fwrapper, 
            dlib::derivative(dlib_fwrapper), 
            x, 
            parameters[0].bounds->min,
            parameters[0].bounds->max
        );

    } else {
        fmin = dlib::find_min_using_approximate_derivatives(
            dlib::bfgs_search_strategy(),
            dlib::objective_delta_stop_strategy(1e-7),
            dlib_fwrapper,
            x,
            -1
        );
    }

    Result res;
    res.fval = fmin;
    res.fevals = fevals;
    res.status = 0;
    for (unsigned int i = 0; i < parameters.size(); i++) {
        FittedParameter param;
        param.name = parameters[i].name;
        param.value = x;
        param.error = {0, 0};
        res.add_parameter(param);
    }
    return res;
}