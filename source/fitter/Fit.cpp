#include <fitter/Fit.h>
#include <iomanip>
#include <sstream>

Fit::Fit(Fitter& fitter, const ROOT::Math::Minimizer* const minimizer, double chi2) : chi2(chi2) {
    unsigned int vars = minimizer->NDim();
    const double* result = minimizer->X();
    const double* errs = minimizer->Errors();
    for (unsigned int i = 0; i < vars; i++) {
        params.insert({minimizer->VariableName(i), result[i]});
        errors.insert({minimizer->VariableName(i), errs[i]});
    }
    add_fit(fitter);

    figures = fitter.plot();
    residuals = fitter.plot_residuals();

    converged = minimizer->Status() == 0;
    calls = minimizer->NCalls();
    dof = fitter.dof() - vars;
}

Fit::Fit(const ROOT::Math::Minimizer* const minimizer, double chi2) : chi2(chi2) {
    unsigned int vars = minimizer->NDim();
    const double* result = minimizer->X();
    const double* errs = minimizer->Errors();
    for (unsigned int i = 0; i < vars; i++) {
        params.insert({minimizer->VariableName(i), result[i]});
        errors.insert({minimizer->VariableName(i), errs[i]});
    }

    if (chi2 == -1) {chi2 = minimizer->MinValue();}
}

Fit::Fit(std::map<std::string, double>& params, std::map<std::string, double>& errs, const double chi2, const int dof, const int calls, const bool converged) : 
    params(params), errors(errs), chi2(chi2), dof(dof), calls(calls), converged(converged) {}

void Fit::add_fit(Fitter& fitter) {
    add_fit(fitter.get_fit());
}

void Fit::add_fit(std::shared_ptr<Fit> fit) {
    for (const auto e : fit->params) {
        params.insert({e.first, e.second});
        errors.insert({e.first, fit->errors.at(e.first)});
    }
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
        ss << "\n| " << print_element(e.first, 8) << " | " << print_element(e.second, 12) << " | " << print_element(errors.at(e.first), 12)  << " |                 |";
    }
    ss << "\n+----------------------------------------------------------+";

    return ss.str();
}