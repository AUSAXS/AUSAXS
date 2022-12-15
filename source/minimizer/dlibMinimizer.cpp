#include <mini/dlibMinimizer.h>
#include <utility/Utility.h>

using namespace mini;

dlibMinimizer::dlibMinimizer() {
    dlib_fwrapper = [this](column_vector x) {return this->function(std::vector<double>(x.begin(), x.end()));};
}

dlibMinimizer::dlibMinimizer(std::function<double(double)> function, Parameter param) : dlibMinimizer() {
    auto f = [=] (std::vector<double> x) {
        return function(x[0]);
    };
    add_parameter(param);
    set_function(f);
}

dlibMinimizer::dlibMinimizer(std::function<double(std::vector<double>)> function, std::vector<Parameter> param) : dlibMinimizer() {
    for (auto& p : param) {
        add_parameter(p);
    }
    set_function(function);
}

dlibMinimizer::~dlibMinimizer() = default;

Result dlibMinimizer::minimize_override() {
    if (!is_parameter_set()) {throw except::bad_order("dlibMinimizer::minimize: No parameters were supplied.");}
    if (!is_function_set()) {throw except::bad_order("dlibMinimizer::minimize: No function was set.");}

    // prepare guess value
    bool bounds = true;
    column_vector x(parameters.size());
    column_vector min(parameters.size());
    column_vector max(parameters.size());
    for (unsigned int i = 0; i < parameters.size(); i++) {
        if (parameters[i].has_guess()) {
            x(i) = parameters[i].guess.value();
        } else if (parameters[i].has_bounds()) {
            x(i) = parameters[i].bounds->center();
        } else {
            throw except::invalid_argument("dlibMinimizer::minimize: Either a guess or bounds must be supplied.");
        }

        if (!parameters[i].has_bounds()) {
            bounds = false;
            if (i != 0) {
                utility::print_warning("dlibMinimizer::minimize_override: Bounds supplied for some parameters, but not all. Disabling bounds.");
            }
        } else {
            min(i) = parameters[i].bounds->min;
            max(i) = parameters[i].bounds->max;
        }
    }

    double fmin;
    if (bounds) {
        fmin = dlib::find_min_box_constrained(
            dlib::bfgs_search_strategy(), 
            dlib::objective_delta_stop_strategy(1e-7), 
            dlib_fwrapper, 
            dlib::derivative(dlib_fwrapper), 
            x, 
            min,
            max
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
        param.value = x(i);
        param.error = {0, 0};
        res.add_parameter(param);
    }
    return res;
}