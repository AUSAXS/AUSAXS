// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_iterative_fit.h>
#include <api/ObjectStorage.h>
#include <data/Molecule.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <fitter/SmartFitter.h>

#include <string>

using namespace ausaxs;
using namespace ausaxs::data;

struct _iterative_fit_state_obj {
    explicit _iterative_fit_state_obj(Molecule* protein) : protein(protein) {}
    Molecule* protein;
    std::vector<double> q, I;
    std::unique_ptr<hist::ICompositeDistanceHistogram> hist; 
    fitter::SmartFitter::EnabledFitParameters enabled_pars = fitter::SmartFitter::EnabledFitParameters::initialize_from_settings();
};
int iterative_fit_init(
    int molecule_id, 
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return -1;}
    molecule->reset_histogram_manager();
    auto obj = _iterative_fit_state_obj(molecule);
    obj.enabled_pars.validate_model(molecule->get_histogram().get());
    obj.hist = molecule->get_histogram();
    return api::ObjectStorage::register_object(std::move(obj));
}, status);}

int iterative_fit_init_userq(
    int molecule_id, 
    double* q, int n_points,
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return -1;}
    molecule->reset_histogram_manager();
    auto obj = _iterative_fit_state_obj(molecule);
    obj.q = std::vector<double>(q, q + n_points);
    obj.hist = molecule->get_histogram();
    obj.enabled_pars.validate_model(molecule->get_histogram().get());
    return api::ObjectStorage::register_object(std::move(obj));
}, status);}

void iterative_fit_evaluate(
    int iterative_fit_id, 
    double* pars, int n_pars, 
    double** return_I, int* n_points,
    int* status
) {return execute_with_catch([&]() {
    auto iterative_fit_state = api::ObjectStorage::get_object<_iterative_fit_state_obj>(iterative_fit_id);
    if (!iterative_fit_state) {ErrorMessage::last_error = "Invalid iterative fit id: \"" + std::to_string(iterative_fit_id) + "\""; return;}
    if (iterative_fit_state->q.empty()) {
        iterative_fit_state->q = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).as_vector();
    }
    auto& enabled_pars = iterative_fit_state->enabled_pars;
    if (n_pars != static_cast<int>(enabled_pars.get_enabled_pars_count())) {
        throw std::runtime_error(
            "Number of provided parameters (" + std::to_string(n_pars) + ") " 
            "does not match number of enabled fit parameters (" + std::to_string(enabled_pars.get_enabled_pars_count()) + ")"
        );
    }
    enabled_pars.apply_pars(
        std::vector<double>(pars, pars+n_pars),
        iterative_fit_state->hist.get()
    );

    iterative_fit_state->I = iterative_fit_state->hist->debye_transform(iterative_fit_state->q).y();
    *return_I = iterative_fit_state->I.data();
    *n_points = static_cast<int>(iterative_fit_state->I.size());
}, status);}

void iterative_fit_evaluate_userq(
    int iterative_fit_id, 
    double* pars, int n_pars, 
    double* q, double* I, int n_points,
    int* status
) {return execute_with_catch([&]() {
    auto iterative_fit_state = api::ObjectStorage::get_object<_iterative_fit_state_obj>(iterative_fit_id);
    if (!iterative_fit_state) {ErrorMessage::last_error = "Invalid iterative fit id: \"" + std::to_string(iterative_fit_id) + "\""; return;}
    auto& enabled_pars = iterative_fit_state->enabled_pars;
    if (n_pars != static_cast<int>(enabled_pars.get_enabled_pars_count())) {
        throw std::runtime_error(
            "Number of provided parameters (" + std::to_string(n_pars) + ") " 
            "does not match number of enabled fit parameters (" + std::to_string(enabled_pars.get_enabled_pars_count()) + ")"
        );
    }
    enabled_pars.apply_pars(
        std::vector<double>(pars, pars+n_pars),
        iterative_fit_state->hist.get()
    );

    std::vector<double> q_vals(q, q + n_points);
    std::vector<double> I_vals = iterative_fit_state->hist->debye_transform(q_vals).y();
    std::copy(I_vals.begin(), I_vals.end(), I);
}, status);}