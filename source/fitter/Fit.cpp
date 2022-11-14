#include <fitter/Fit.h>
#include <fitter/Fitter.h>
#include <utility/Exceptions.h>

#include <iomanip>
#include <sstream>

Fit::Fit(Fitter& fitter, const mini::Result& res, double chi2) noexcept : Fit(res, chi2, 0) {
    dof = fitter.dof();
    add_fit(fitter);
    add_plots(fitter);
}

Fit::Fit(const mini::Result& res, double chi2, double dof) noexcept : Result(res), dof(dof) {}

void Fit::add_plots(Fitter& fitter) {
    figures = fitter.plot();
    residuals = fitter.plot_residuals();
}

void Fit::add_fit(Fitter& fitter) noexcept {
    add_fit(fitter.get_fit());
}

void Fit::add_fit(std::shared_ptr<Fit> fit) noexcept {
    for (const auto& e : fit->parameters) {
        parameters.push_back(e);
        dof--;
    }
}

template<typename T>
struct print_element {
    print_element(T t, int width) : t(t), width(width) {}

    friend std::ostream& operator<<(std::ostream& os, const print_element<T> e) noexcept {
        std::stringstream ss; ss << e.t;
        std::string val = ss.str();
        if (val.size() > e.width) {val = val.substr(0, e.width);}

        os << std::left << std::setw(e.width) << e.t; return os;
    }

    T t;
    unsigned int width;
};

std::string Fit::to_string() const noexcept {
    std::stringstream ss;
    ss << "+----------------------------------------------------------+"
       << "\n|                       FIT REPORT                         |"
       << "\n+----------------------------------------------------------+"
       << "\n| Converged: " << (status == 0 ? "yes" : "no ") << "                              Fevals: " << print_element(fevals, 4) << " |"
       << "\n| chi2: " << print_element(fval, 10) << "   dof: " << print_element(dof, 6) << "    chi2/dof: " << print_element(fval/dof, 12) << " |"
       << "\n+----------------------------------------------------------+"
       << "\n| PAR      | VAL          | UNC          |                 |";
    for (const auto& e : parameters) {
        ss << "\n| " << print_element(e.name, 8) << " | " << print_element(e.value, 12) << " | " << print_element(e.mean_error(), 12)  << " |                 |";
    }
    ss << "\n+----------------------------------------------------------+";

    return ss.str();
}

std::string EMFit::to_string() const noexcept {
    std::stringstream ss;
    ss << Fit::to_string();
    ss << "\n| Cutoff corresponds to PyMOL level " << print_element(level, 12) << "           |";
    ss << "\n+----------------------------------------------------------+";
}