// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_molecule.h>
#include <api/ObjectStorage.h>
#include <io/Reader.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/intensity_calculator/ExactDebyeCalculator.h>
#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/histogram_manager/IHistogramManager.h>
#include <hist/distribution/Distribution1D.h>
#include <hist/detail/SimpleExvModel.h>
#include <fitter/SmartFitter.h>
#include <settings/All.h>

#include <string>

using namespace ausaxs;
using namespace ausaxs::data;

int molecule_from_file(const char* filename, int* status) {return execute_with_catch([&]() {
    auto molecule = data::Molecule(std::string(filename));
    auto molecule_id = api::ObjectStorage::register_object(std::move(molecule));
    return molecule_id;
}, status);}

int molecule_from_pdb_id(int pdb_id, int* status) {return execute_with_catch([&]() {
    auto pdb = api::ObjectStorage::get_object<io::pdb::PDBStructure>(pdb_id);
    if (!pdb) {ErrorMessage::last_error = "Invalid pdb id: \"" + std::to_string(pdb_id) + "\""; return -1;}
    if (settings::molecule::implicit_hydrogens) {pdb->add_implicit_hydrogens();}
    auto data = pdb->reduced_representation();
    auto molecule = data.waters.empty() 
        ? Molecule({Body{std::move(data.atoms)}})
        : Molecule({Body{std::move(data.atoms), std::move(data.waters)}})
    ;
    auto molecule_id = api::ObjectStorage::register_object(std::move(molecule));
    return molecule_id;
}, status);}

int molecule_from_arrays(double* xx, double* yy, double* zz, double* ww, int n_atoms, int* status) {return execute_with_catch([&]() {
    std::vector<data::Atom> atoms(n_atoms);
    std::vector<double> x(xx, xx + n_atoms);
    std::vector<double> y(yy, yy + n_atoms);
    std::vector<double> z(zz, zz + n_atoms);
    std::vector<double> w(ww, ww + n_atoms);
    for (int i = 0; i < n_atoms; ++i) {
        atoms[i] = data::Atom({x[i], y[i], z[i]}, w[i]);
    }
    auto molecule = Molecule({Body{atoms}});
    auto molecule_id = api::ObjectStorage::register_object(std::move(molecule));
    return molecule_id;
}, status);}

struct _molecule_get_data_obj {
    _molecule_get_data_obj(unsigned int n_atoms, unsigned int n_waters) :
        ax(n_atoms), ay(n_atoms), az(n_atoms), aw(n_atoms),
        wx(n_waters), wy(n_waters), wz(n_waters), ww(n_waters), 
        aform_factors(n_atoms), aform_factors_ptr(n_atoms)
    {}
    std::vector<double> ax, ay, az, aw, wx, wy, wz, ww;
    std::vector<std::string> aform_factors;
    std::vector<const char*> aform_factors_ptr;
};
int molecule_get_data(
    int molecule_id,
    double** ax_out, double** ay_out, double** az_out, double** aw_out, const char*** aform_factors_out,
    double** wx_out, double** wy_out, double** wz_out, double** ww_out,
    int* na, int* nw, int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return -1;}
    auto atoms = molecule->get_atoms();
    auto waters = molecule->get_waters();
    _molecule_get_data_obj data(atoms.size(), waters.size());
    for (unsigned int i = 0; i < atoms.size(); ++i) {
        const auto& atom = atoms[i];
        data.ax[i] = atom.coordinates().x();
        data.ay[i] = atom.coordinates().y();
        data.az[i] = atom.coordinates().z();
        data.aw[i] = atom.weight();
        data.aform_factors[i] = form_factor::to_string(atom.form_factor_type());
        data.aform_factors_ptr[i] = data.aform_factors[i].c_str();
    }
    for (unsigned int i = 0; i < waters.size(); ++i) {
        const auto& water = waters[i];
        data.wx[i] = water.coords.x();
        data.wy[i] = water.coords.y();
        data.wz[i] = water.coords.z();
        data.ww[i] = water.weight();
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_molecule_get_data_obj>(data_id);
    *ax_out = ref->ax.data();
    *ay_out = ref->ay.data();
    *az_out = ref->az.data();
    *aw_out = ref->aw.data();
    *wx_out = ref->wx.data();
    *wy_out = ref->wy.data();
    *wz_out = ref->wz.data();
    *ww_out = ref->ww.data();
    *aform_factors_out = ref->aform_factors_ptr.data();
    *na = static_cast<int>(atoms.size());
    *nw = static_cast<int>(waters.size());
    return data_id;
}, status);}

void molecule_hydrate(
    int molecule_id,
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return;}
    molecule->generate_new_hydration();
}, status);}

struct _molecule_distance_histogram_obj {
    explicit _molecule_distance_histogram_obj(unsigned int n_bins) : aa(n_bins), aw(n_bins), ww(n_bins), axis(n_bins) {}
    std::vector<double> aa, aw, ww, axis;
};
int molecule_distance_histogram(
    int molecule_id,
    double** aa, double** aw, double** ww, double** axis, int* n_bins, 
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return -1;}
    molecule->reset_histogram_manager();
    auto hist = molecule->get_histogram();
    _molecule_distance_histogram_obj data(settings::axes::bin_count);
    {   // copy to avoid issues with mismatching sizes
        //? is this necessary? I think these should always have the same size (or at least may be truncated to the same length)
        auto& aa_dist = hist->get_aa_counts();
        auto& aw_dist = hist->get_aw_counts();
        auto& ww_dist = hist->get_ww_counts();
        auto& axis = hist->get_d_axis();
        std::copy(aa_dist.get_content().begin(), aa_dist.get_content().end(), data.aa.begin());
        std::copy(aw_dist.get_content().begin(), aw_dist.get_content().end(), data.aw.begin());
        std::copy(ww_dist.get_content().begin(), ww_dist.get_content().end(), data.ww.begin());
        std::copy(axis.begin(), axis.end(), data.axis.begin());
    }
    assert(data.aa.size() == settings::axes::bin_count);
    assert(data.aw.size() == settings::axes::bin_count);
    assert(data.ww.size() == settings::axes::bin_count);
    assert(data.axis.size() == settings::axes::bin_count);

    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_molecule_distance_histogram_obj>(data_id);
    *aa = ref->aa.data();
    *aw = ref->aw.data();
    *ww = ref->ww.data();
    *axis = ref->axis.data();
    *n_bins = static_cast<int>(settings::axes::bin_count);
    return data_id;
}, status);}

struct _molecule_debye_obj {
    explicit _molecule_debye_obj(unsigned int size) :
        q(size), I(size)
    {}
    std::vector<double> q, I;
};
int molecule_debye(
    int molecule_id, 
    double** q, double** I, int* n_points,
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return -1;}
    molecule->reset_histogram_manager();
    auto hist = molecule->get_histogram();
    auto debye_I = hist->debye_transform();
    _molecule_debye_obj data(debye_I.size());
    for (unsigned int i = 0; i < debye_I.size(); ++i) {
        data.q[i] = constants::axes::q_vals[i];
        data.I[i] = debye_I[i];
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_molecule_debye_obj>(data_id);
    *q = ref->q.data();
    *I = ref->I.data();
    *n_points = static_cast<int>(debye_I.size());
    return data_id;
}, status);}

void molecule_debye_userq(
    int molecule_id, 
    double* q, double* I, int n_points,
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return;}
    molecule->reset_histogram_manager();
    std::vector<double> q_vals(q, q + n_points);
    auto hist = molecule->get_histogram();
    auto debye_I = hist->debye_transform(q_vals);
    if (static_cast<int>(debye_I.size()) != n_points) {*status = 3; return;}
    for (int i = 0; i < n_points; ++i) {
        I[i] = debye_I.y(i);
    }
}, status);}

int molecule_debye_raw(
    int molecule_id,
    double** q, double** I, int* n_points,
    int* status
) {return execute_with_catch([&]() {
    hist::detail::SimpleExvModel::disable(); // disable exv contributions to HistogramManager

    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return -1;}
    auto hist = hist::factory::construct_histogram_manager(molecule, settings::hist::HistogramManagerChoice::HistogramManagerMT)->calculate();
    auto debye_I = hist->debye_transform();
    _molecule_debye_obj data(debye_I.size());
    for (unsigned int i = 0; i < debye_I.size(); ++i) {
        data.q[i] = constants::axes::q_vals[i];
        data.I[i] = debye_I[i]*std::exp(data.q[i]*data.q[i]); // remove form factor added by debye transform
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_molecule_debye_obj>(data_id);
    *q = ref->q.data();
    *I = ref->I.data();
    *n_points = static_cast<int>(debye_I.size());

    hist::detail::SimpleExvModel::enable(); // re-enable exv contributions to ensure consistency elsewhere
    return data_id;
}, status);}

void molecule_debye_raw_userq(
    int molecule_id, 
    double* q, double* I, int n_points,
    int* status
) {return execute_with_catch([&]() {
    hist::detail::SimpleExvModel::disable(); // disable exv contributions to HistogramManager

    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return;}
    std::vector<double> q_vals(q, q + n_points);
    auto hist = hist::factory::construct_histogram_manager(molecule, settings::hist::HistogramManagerChoice::HistogramManagerMT)->calculate();
    auto debye_I = hist->debye_transform(q_vals);
    if (static_cast<int>(debye_I.size()) != n_points) {*status = 3; return;}
    for (int i = 0; i < n_points; ++i) {
        I[i] = debye_I.y(i)*std::exp(q_vals[i]*q_vals[i]); // remove form factor added by debye transform
    }
    hist::detail::SimpleExvModel::enable(); // re-enable exv contributions to ensure consistency elsewhere
}, status);}

int molecule_debye_exact(
    int molecule_id,
    double** q, double** I, int* n_points,
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return -1;}
    auto qv = constants::axes::q_axis.sub_axis(settings::axes::qmin, settings::axes::qmax).as_vector();
    auto Iq = hist::exact_debye_transform(*molecule, qv);
    _molecule_debye_obj data(Iq.size());
    for (unsigned int i = 0; i < Iq.size(); ++i) {
        data.q[i] = qv[i];
        data.I[i] = Iq[i]*std::exp(qv[i]*qv[i]); // remove form factor added by exact_debye
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_molecule_debye_obj>(data_id);
    *q = ref->q.data();
    *I = ref->I.data();
    *n_points = static_cast<int>(Iq.size());
    return data_id;
}, status);}

void molecule_debye_exact_userq(
    int molecule_id, 
    double* q, double* I, int n_points,
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return;}
    std::vector<double> q_vals(q, q + n_points);
    auto Iq = hist::exact_debye_transform(*molecule, q_vals);
    if (static_cast<int>(Iq.size()) != n_points) {*status = 3; return;}
    for (int i = 0; i < n_points; ++i) {
        I[i] = Iq[i]*std::exp(q_vals[i]*q_vals[i]); // remove form factor added by exact_debye
    }
}, status);}

int molecule_debye_fit(
    int molecule_id, int data_id,
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return -1;}
    molecule->reset_histogram_manager();
    auto dataset = api::ObjectStorage::get_object<SimpleDataset>(data_id);
    if (!dataset) {ErrorMessage::last_error = "Invalid dataset id: \"" + std::to_string(data_id) + "\""; return -1;}
    auto fitter = fitter::SmartFitter(*dataset, molecule->get_histogram());
    int fit_result_id = api::ObjectStorage::register_object(fitter.fit());
    return fit_result_id;
}, status);}

void molecule_clear_hydration(
    int molecule_id,
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return;}
    molecule->clear_hydration();
}, status);}

void molecule_Rg(
    int molecule_id,
    double* Rg,
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return;}
    *Rg = molecule->get_Rg();
}, status);}