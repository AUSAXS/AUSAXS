/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/FitResult.h>
#include <fitter/Fitter.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>

#include <sstream>

using namespace fitter;

FitResult::FitResult(const mini::Result& res, double chi2, unsigned int dof) noexcept : Result(res), dof(dof) {
    fval = chi2;
}

void FitResult::add_fit(observer_ptr<FitResult> fit, bool front) noexcept {
    if (front) {
        parameters.insert(parameters.begin(), fit->parameters.begin(), fit->parameters.end());
    } else {
        parameters.insert(parameters.end(), fit->parameters.begin(), fit->parameters.end());
    }
    dof -= fit->parameters.size();
}

void FitResult::set_data_curves(Dataset&& curves_) {
    curves = std::move(curves_);
    if (curves.size() != 5) {throw except::invalid_argument("FitResult::set_data_curves: Invalid number of columns. Expected | q | I | I_err | I_fit | residuals |.");}
    if (curves.is_named() && curves.get_col_names() != std::vector<std::string>{"q", "I", "I_err", "I_fit", "residuals"}) {
        throw except::invalid_argument("FitResult::set_data_curves: Invalid column names. Expected | q | I | I_err | I_fit | residuals |.");
    }
    else {curves.set_col_names({"q", "I", "I_err", "I_fit", "residuals"});}
}

void FitResult::set_data_curves(std::vector<double>&& q, std::vector<double>&& data, std::vector<double>&& data_err, std::vector<double>&& model, std::vector<double>&& residuals) {
    curves = Dataset({q, data, data_err, model, residuals}, {"q", "I", "I_err", "I_fit", "residuals"});
}

std::string FitResult::to_string() const noexcept {
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