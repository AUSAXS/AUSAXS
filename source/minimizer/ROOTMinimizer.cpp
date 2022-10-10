#include <minimizer/ROOTMinimizer.h>
#include <utility/Exceptions.h>
#include <utility/Settings.h>

using namespace mini;

ROOTMinimizer::ROOTMinimizer(std::string package, std::string algorithm) {
    create_minimizer(package, algorithm);
}

ROOTMinimizer::ROOTMinimizer(std::string package, std::string algorithm, double(&function)(const double*), Parameter param) {
        create_minimizer(package, algorithm);
        set_function(function);
        if (!param.empty()) {add_parameter(param);}
}

ROOTMinimizer::ROOTMinimizer(std::string package, std::string algorithm, std::function<double(const double*)> function, Parameter param) {
        create_minimizer(package, algorithm);
        set_function(function);
        if (!param.empty()) {add_parameter(param);}
}

ROOTMinimizer::ROOTMinimizer(std::string package, std::string algorithm, double(&function)(const double*), std::vector<Parameter> params) {
        create_minimizer(package, algorithm);
        set_function(function);
        std::for_each(params.begin(), params.end(), [this] (const Parameter& param) {add_parameter(param);});
}

ROOTMinimizer::ROOTMinimizer(std::string package, std::string algorithm, std::function<double(const double*)> function, std::vector<Parameter> params) {
        create_minimizer(package, algorithm);
        set_function(function);
        std::for_each(params.begin(), params.end(), [this] (const Parameter& param) {add_parameter(param);});
}

ROOTMinimizer::~ROOTMinimizer() {
    mini->~Minimizer();
}

void ROOTMinimizer::create_minimizer(std::string package, std::string algorithm) {
    mini = ROOT::Math::Factory::CreateMinimizer(package, algorithm);
    if (mini == nullptr) {
        throw except::invalid_argument("Error in ROOTMinimizer::create_minimizer: \"" + package + "\" and \"" + algorithm + "\" is not a valid algorithm.");
    }
}

Dataset2D ROOTMinimizer::get_evaluated_points() const {
    if (evaluations.empty()) {throw except::bad_order("Error in ROOTMinimizer::get_evaluated_points: Cannot get evaluated points before a minimization call has been made.");}

    unsigned int N = evaluations.size();
    std::vector<double> x(N), y(N);
    for (unsigned int i = 0; i < N; i++) {
        x[i] = evaluations[i].vals[0];
        y[i] = evaluations[i].fval;
    }
    return Dataset2D(x, y, "x", "f(x)");
}

Result ROOTMinimizer::minimize_override() {
    if (!is_parameter_set()) {throw except::bad_order("Error in ROOTMinimizer::minimize: No parameters were supplied.");}
    if (!is_function_set()) {throw except::bad_order("Error in ROOTMinimizer::minimize: No function was set.");}
    prepare_minimizer();
    if (!setting::fit::verbose) {
        mini->SetPrintLevel(-1);
    }
    mini->Minimize();
    if (setting::fit::verbose) {mini->PrintResults();}

    Result res;
    res.fval = mini->MinValue();
    res.fevals = fevals;
    res.status = mini->Status();
    unsigned int vars = mini->NDim();
    const double* result = mini->X();
    const double* error = mini->Errors();

    for (unsigned int i = 0; i < vars; i++) {
        FittedParameter param;
        param.name = mini->VariableName(i);
        param.value = result[i];
        if (mini->ProvidesError()) {
            param.error = {-error[i], error[i]};
        } else {
            param.error = {0, 0};
        }
        res.add_parameter(param);
    }
    return res;
}

Dataset2D ROOTMinimizer::landscape(unsigned int) {
    throw except::unexpected("Error in ROOTMinimizer::landscape: Not implemented yet.");
}

void ROOTMinimizer::prepare_minimizer() {
    if (mini == nullptr) {throw except::unexpected("Error in ROOTMinimizer::prepare_minimizer: Minimizer has not been initialized.");}
    functor = ROOT::Math::Functor(function, parameters.size());
    mini->SetFunction(functor);
    mini->SetTolerance(tol);

    unsigned int var_counter = 0;
    for (const Parameter& param : parameters) {
        if (param.has_bounds()) {
            double guess = param.has_guess() ? param.guess.value() : param.bounds.value().center();
            mini->SetLimitedVariable(var_counter++, param.name, guess, 0.1, param.bounds.value().min, param.bounds.value().max);
        } else if (param.has_guess()) {
            mini->SetVariable(var_counter++, param.name, param.guess.value(), 0.1);
        } else {
            throw except::unexpected("Error in ROOTMinimizer::prepare_minimizer: This shouldn't be able to happen.");
        }
    }
}