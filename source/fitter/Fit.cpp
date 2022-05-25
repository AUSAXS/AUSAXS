#include <fitter/Fit.h>
#include <fitter/Fitter.h>
#include <utility/Exceptions.h>

#include <iomanip>
#include <sstream>

Fit::Fit(Fitter& fitter, const mini::Result& res, double chi2) : Fit(res, chi2, 0) {
    dof = fitter.dof() - res.dim();

    add_fit(fitter);
    figures = fitter.plot();
    residuals = fitter.plot_residuals();
}

Fit::Fit(const mini::Result& res, double chi2, double dof) : chi2(chi2), dof(dof) {
    params = res.parameters;
    calls = res.fevals;
    converged = res.status == 0;
}

void Fit::add_fit(Fitter& fitter) {
    add_fit(fitter.get_fit());
}

void Fit::add_fit(std::shared_ptr<Fit> fit) {
    for (const auto e : fit->params) {
        params.push_back(e);
    }
}

const mini::FittedParameter& Fit::get_parameter(unsigned int index) const {
    if (params.size() <= index) {throw except::out_of_bounds("Error in Fit::get_parameter: Index \"" + std::to_string(index) + "\" is out of range (" + std::to_string(params.size()) + ").");}
    return params[index];
} 

mini::FittedParameter& Fit::get_parameter(unsigned int index) {
    return const_cast<mini::FittedParameter&>(std::as_const(*this).get_parameter(index));
} 

const mini::FittedParameter& Fit::get_parameter(std::string name) const {
    auto pos = std::find_if(params.begin(), params.end(), [&name] (const mini::FittedParameter& param) {return param.name == name;});
    if (pos == params.end()) {throw except::unknown_argument("Error in Fit::get_parameter: No parameter named \"" + name + "\" was found.");}
    return *pos;
}

mini::FittedParameter& Fit::get_parameter(std::string name) {
    return const_cast<mini::FittedParameter&>(std::as_const(*this).get_parameter(name));
}

template<typename T>
struct print_element {
    print_element(T t, int width) : t(t), width(width) {}

    friend std::ostream& operator<<(std::ostream& os, const print_element<T> e) {
        std::stringstream ss; ss << e.t;
        std::string val = ss.str();
        if (val.size() > e.width) {val = val.substr(0, e.width);}

        os << std::left << std::setw(e.width) << e.t; return os;
    }

    T t;
    unsigned int width;
};

std::string Fit::to_string() const {
    std::stringstream ss;
    ss << "+----------------------------------------------------------+"
       << "\n|                       FIT REPORT                         |"
       << "\n+----------------------------------------------------------+"
       << "\n| Converged: " << (converged ? "yes" : "no ") << "                              Fevals: " << print_element(calls, 4) << " |"
       << "\n| chi2: " << print_element(chi2, 10) << "   dof: " << print_element(dof, 6) << "   chi2/dof: " << print_element(chi2/dof, 12) << " |"
       << "\n+----------------------------------------------------------+"
       << "\n| PAR      | VAL          | UNC          |                 |";
    for (const auto& e : params) {
        ss << "\n| " << print_element(e.name, 8) << " | " << print_element(e.value, 12) << " | " << print_element(e.mean_error(), 12)  << " |                 |";
    }
    ss << "\n+----------------------------------------------------------+";

    return ss.str();
}