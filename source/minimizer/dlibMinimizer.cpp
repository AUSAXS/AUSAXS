#include <mini/dlibMinimizer.h>

using namespace mini;

dlibMinimizer::dlibMinimizer(std::function<double(double)> function, Parameter param) {
    auto f = [=] (std::vector<double> x) {
        return function(x[0]);
    };
    add_parameter(param);
    set_function(f);
    dlib_fwrapper = [this](column_vector x) {return this->function({x(0)});};
}

dlibMinimizer::dlibMinimizer(std::function<double(std::vector<double>)> function, Parameter param) {
    add_parameter(param);
    set_function(function);
    dlib_fwrapper = [this](column_vector x) {return this->function({x(0)});};
}

dlibMinimizer::~dlibMinimizer() = default;

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