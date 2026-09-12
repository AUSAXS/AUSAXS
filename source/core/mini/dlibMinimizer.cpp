// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#if defined(DLIB_AVAILABLE)
    #include <mini/detail/FittedParameter.h>
    #include <mini/detail/Parameter.h>
    #include <mini/dlibMinimizer.h>
    #include <utility/Console.h>

    #include <dlib/global_optimization.h>
    #include <dlib/optimization.h>

    using namespace ausaxs;
    using namespace ausaxs::mini;
    namespace ausaxs::mini {
        struct column_vector : dlib::matrix<double, 0, 1> {
            using dlib::matrix<double, 0, 1>::matrix;
        };
    }

    template<mini::algorithm algo>
    dlibMinimizer<algo>::dlibMinimizer() = default;

    template<mini::algorithm algo>
    dlibMinimizer<algo>::dlibMinimizer(std::function<double(double)> function, const Parameter& param) {
        auto f = [_function = std::move(function)] (std::vector<double> x) {
            return _function(x[0]);
        };
        add_parameter(param);
        set_function(std::move(f));
    }

    template<mini::algorithm algo>
    dlibMinimizer<algo>::dlibMinimizer(std::function<double(std::vector<double>)> function, const std::vector<Parameter>& param) {
        for (const auto& p : param) {
            add_parameter(p);
        }
        set_function(std::move(function));
    }

    template<mini::algorithm algo>
    dlibMinimizer<algo>::~dlibMinimizer() = default;

    template<mini::algorithm algo>
    Result dlibMinimizer<algo>::minimize_override() {
        if (!is_parameter_set()) {throw ausaxs::except::invalid_argument("dlibMinimizer::minimize: No parameters were supplied.");}
        if (!is_function_set()) {throw ausaxs::except::invalid_argument("dlibMinimizer::minimize: No function was set.");}

        // prepare guess value
        bool bounds = true;
        column_vector x(parameters.size());
        column_vector min(parameters.size());
        column_vector max(parameters.size());
        for (int i = 0; i < static_cast<int>(parameters.size()); i++) {
            const auto& param = parameters[i];
            if (param.guess.has_value()) {
                x(i) = *param.guess;
            } else if (param.bounds.has_value()) {
                x(i) = param.bounds->center();
            } else {
                throw ausaxs::except::invalid_argument("dlibMinimizer::minimize: Either a guess or bounds must be supplied.");
            }

            if (param.bounds.has_value()) {
                min(i) = param.bounds->min;
                max(i) = param.bounds->max;
            } else {
                bounds = false;
                if (i != 0) {
                    console::print_warning("dlibMinimizer::minimize_override: Bounds supplied for some parameters, but not all. Disabling bounds.");
                }
            }
        }

        double fmin;
        auto fwrapper = [this](dlib::matrix<double, 0, 1> x) {return this->function(std::vector<double>(x.begin(), x.end()));};
        if (bounds) {
            if constexpr (algo == mini::algorithm::DLIB_GLOBAL) {
                auto eval = dlib::find_min_global(
                    fwrapper, 
                    min, 
                    max,
                    dlib::max_function_calls(max_evals)
                );
                x = eval.x;
                fmin = eval.y;
            } else if constexpr (algo == mini::algorithm::BFGS) {
                fmin = dlib::find_min_box_constrained(
                    dlib::bfgs_search_strategy(), 
                    dlib::objective_delta_stop_strategy(1e-7), 
                    fwrapper, 
                    dlib::derivative(fwrapper), 
                    x, 
                    min,
                    max
                );
            }
        } else {
            console::print_warning("dlibMinimizer::minimize_override: No bounds supplied. Using unconstrained minimization.");
            fmin = dlib::find_min_using_approximate_derivatives(
                dlib::bfgs_search_strategy(),
                dlib::objective_delta_stop_strategy(1e-7),
                fwrapper,
                x,
                -1
            );
        }

        Result res;
        res.fval = fmin;
        res.fevals = fevals;
        res.status = 0;
        for (int i = 0; i < static_cast<int>(parameters.size()); i++) {
            FittedParameter param;
            param.name = parameters[i].name;
            param.value = x(i);
            param.error = {0, 0};
            res.add_parameter(param);
        }
        return res;
    }

    template class mini::dlibMinimizer<mini::algorithm::DLIB_GLOBAL>;
    template class mini::dlibMinimizer<mini::algorithm::BFGS>;
#endif