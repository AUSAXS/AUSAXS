#include <fitter/Fit.h>
#include <fitter/Fitter.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>

#include <sstream>

using namespace fitter;

Fit::Fit(Fitter& fitter, const mini::Result& res, double chi2) noexcept : Fit(res, chi2, fitter.dof()) {
    add_fit(fitter);
    add_plots(fitter);
}

Fit::Fit(const mini::Result& res, double chi2, double dof) noexcept : Result(res), dof(dof) {
    fval = chi2;
}

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

bool Fit::operator==(const Fit& other) const = default;

std::string Fit::to_string() const noexcept {
    std::stringstream ss;
    ss <<   "+----------------------------------------------------------+"
       << "\n|                       FIT REPORT                         |"
       << "\n+----------------------------------------------------------+"
       << "\n| Converged: " << (status == 0 ? "yes" : "no ") << "                              Fevals: " << utility::print_element(fevals, 4) << " |"
       << "\n| chi2: " << utility::print_element(fval, 10) << "   dof: " << utility::print_element(dof, 6) << "    chi2/dof: " << utility::print_element(fval/dof, 12) << " |"
       << "\n+----------------------------------------------------------+"
       << "\n| PAR      | VAL          | UNC          |                 |";
    for (const auto& e : parameters) {
        ss << "\n| " << utility::print_element(e.name, 8) << " | " << utility::print_element(e.value, 12) << " | " << utility::print_element(e.mean_error(), 12)  << " |                 |";
    }
    ss << "\n+----------------------------------------------------------+";

    return ss.str();
}

std::string EMFit::to_string() const noexcept {
    std::stringstream ss;
    ss << Fit::to_string();
    ss << "\n| Cutoff corresponds to PyMOL level " << utility::print_element(level, 12) << "           |";
    ss << "\n+----------------------------------------------------------+";
    return ss.str();
}