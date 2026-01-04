// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/api_pyausaxs.h>
#include <api/ObjectStorage.h>
#include <fitter/SmartFitter.h>

#include <string>

using namespace ausaxs;

struct _fit_get_fit_info_obj {
    explicit _fit_get_fit_info_obj(unsigned int n_pars) : 
        pars(n_pars), pars_ptr(n_pars), pvals(n_pars), perr_n(n_pars), perr_p(n_pars)
    {}
    std::vector<std::string> pars;
    std::vector<const char*> pars_ptr;
    std::vector<double> pvals, perr_n, perr_p;
};
int fit_get_fit_info(
    int fit_id,
    const char*** pars, double** pvals, double** perr_min, double** perr_max, int* n_pars,
    double* chi_squared, int* dof,
    int* status
) {return execute_with_catch([&]() {
    auto fit_result = api::ObjectStorage::get_object<fitter::FitResult>(fit_id);
    if (!fit_result) {ErrorMessage::last_error = "Invalid fit result id: \"" + std::to_string(fit_id) + "\""; return -1;}

    _fit_get_fit_info_obj data(fit_result->parameters.size());
    for (unsigned int i = 0; i < fit_result->parameters.size(); ++i) {
        const auto& par = fit_result->parameters[i];
        data.pars[i] = par.name;
        data.pvals[i] = par.value;
        data.perr_n[i] = par.error.min;
        data.perr_p[i] = par.error.max;
        data.pars_ptr[i] = data.pars[i].c_str();
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_fit_get_fit_info_obj>(data_id);

    *dof = fit_result->dof;
    *chi_squared = fit_result->fval;
    *n_pars = static_cast<int>(fit_result->parameters.size());
    *pars = ref->pars_ptr.data();
    *pvals = ref->pvals.data();
    *perr_min = ref->perr_n.data();
    *perr_max = ref->perr_p.data();
    return data_id;
}, status);}

struct _fit_get_fit_curves_obj {
    explicit _fit_get_fit_curves_obj(unsigned int size) :
        q(size), I_data(size), I_err(size), I_model(size)
    {}
    std::size_t size() const {return q.size();}
    std::vector<double> q, I_data, I_err, I_model;
};
int fit_get_fit_curves(
    int fit_id,
    double** q, double** I_data, double** I_err, double** I_model, int* n_points,
    int* status
) {return execute_with_catch([&]() {
    auto fit_result = api::ObjectStorage::get_object<fitter::FitResult>(fit_id);
    if (!fit_result) {ErrorMessage::last_error = "Invalid fit result id: \"" + std::to_string(fit_id) + "\""; return -1;}
    _fit_get_fit_curves_obj data(fit_result->curves.size_rows());
    for (unsigned int i = 0; i < data.size(); ++i) {
        data.q[i]       = fit_result->curves.col(0)[i];
        data.I_data[i]  = fit_result->curves.col(1)[i];
        data.I_err[i]   = fit_result->curves.col(2)[i];
        data.I_model[i] = fit_result->curves.col(3)[i];
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_fit_get_fit_curves_obj>(data_id);
    *q = ref->q.data();
    *I_data = ref->I_data.data();
    *I_err = ref->I_err.data();
    *I_model = ref->I_model.data();
    *n_points = static_cast<int>(ref->size());
    return data_id;
}, status);}