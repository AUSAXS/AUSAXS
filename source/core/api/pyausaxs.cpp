// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs.h>
#include <api/ObjectStorage.h>
#include <api/ErrorMessage.h>
#include <io/detail/PDBReader.h>
#include <io/Reader.h>
#include <dataset/SimpleDataset.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/distribution/Distribution1D.h>
#include <fitter/SmartFitter.h>
#include <settings/MoleculeSettings.h>
#include <settings/SettingsIO.h>

#include <string>
#include <type_traits>

using namespace ausaxs;
using namespace ausaxs::data;

template<typename Func>
auto execute_with_catch(Func&& f, int* status) -> decltype(f()) {
    try {
        *status = 1;
        if constexpr (std::is_void_v<decltype(f())>) {
            f();
            *status = 0;
            return;
        } else {
            auto v = f();
            *status = 0;
            return v;
        }
    } catch (const std::exception& e) {
        ErrorMessage::last_error = e.what();
        *status = 1;
    } catch (...) {
        ErrorMessage::last_error = "An unknown error occurred.";
        *status = 1;
    }
    if constexpr (!std::is_void_v<decltype(f())>) {return {};}
}

int read_pdb(
    const char* filename,
    int* status
) {return execute_with_catch([&]() {
    auto pdb = io::detail::pdb::read(std::string(filename));
    auto pdb_id = api::ObjectStorage::register_object(std::move(pdb));
    return pdb_id;
}, status);}

int pdb_read(
    const char* filename,
    int* status
) {return execute_with_catch([&]() {
    auto pdb = io::Reader::read(std::string(filename));
    auto pdb_id = api::ObjectStorage::register_object(std::move(pdb));
    return pdb_id;
}, status);}

struct _pdb_get_data_obj {
    _pdb_get_data_obj(unsigned int size) :
        serial(size), resSeq(size), name(size), altLoc(size), resName(size), iCode(size), element(size), charge(size), 
        name_ptr(size), altLoc_ptr(size), resName_ptr(size), iCode_ptr(size), element_ptr(size), charge_ptr(size),
        chainID(size), x(size), y(size), z(size), occupancy(size), tempFactor(size) 
    {}
    std::vector<int> serial, resSeq;
    std::vector<std::string> name, altLoc, resName, iCode, element, charge;
    std::vector<const char*> name_ptr, altLoc_ptr, resName_ptr, iCode_ptr, element_ptr, charge_ptr;
    std::vector<char> chainID;
    std::vector<double> x, y, z, occupancy, tempFactor;
};
int pdb_get_data(
    int object_id,
    int** serial_out, const char*** name_out, const char*** altLoc_out, const char*** resName_out, const char** chainID_out, int** resSeq_out, 
    const char*** iCode_out, double** x_out, double** y_out, double** z_out, double** occupancy_out, double** tempFactor_out, const char*** element_out, 
    const char*** charge_out, int* n_atoms_out, int* status
) {return execute_with_catch([&]() {
    auto pdb = api::ObjectStorage::get_object<io::pdb::PDBStructure>(object_id);
    if (!pdb) {ErrorMessage::last_error = "Invalid pdb id: \"" + std::to_string(object_id) + "\""; return -1;}
    const auto& atoms = pdb->atoms;
    _pdb_get_data_obj data(atoms.size());
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
        const auto& atom = atoms[i];
        data.serial[i] = atom.serial;
        data.name[i] = atom.name;
        data.altLoc[i] = atom.altLoc;
        data.resName[i] = atom.resName;
        data.chainID[i] = atom.chainID;
        data.resSeq[i] = atom.resSeq;
        data.iCode[i] = atom.iCode;
        data.x[i] = atom.coords.x();
        data.y[i] = atom.coords.y();
        data.z[i] = atom.coords.z();
        data.occupancy[i] = atom.occupancy;
        data.tempFactor[i] = atom.tempFactor;
        data.element[i] = constants::symbols::to_string(atom.element);
        data.charge[i] = atom.charge;

        data.name_ptr[i] = data.name[i].c_str();
        data.altLoc_ptr[i] = data.altLoc[i].c_str();
        data.resName_ptr[i] = data.resName[i].c_str();
        data.iCode_ptr[i] = data.iCode[i].c_str();
        data.element_ptr[i] = data.element[i].c_str();
        data.charge_ptr[i] = data.charge[i].c_str();
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_pdb_get_data_obj>(data_id);
    *serial_out = ref->serial.data();
    *name_out = ref->name_ptr.data();
    *altLoc_out = ref->altLoc_ptr.data();
    *resName_out = ref->resName_ptr.data();
    *chainID_out = ref->chainID.data();
    *resSeq_out = ref->resSeq.data();
    *iCode_out = ref->iCode_ptr.data();
    *x_out = ref->x.data();
    *y_out = ref->y.data();
    *z_out = ref->z.data();
    *occupancy_out = ref->occupancy.data();
    *tempFactor_out = ref->tempFactor.data();
    *element_out = ref->element_ptr.data();
    *charge_out = ref->charge_ptr.data();
    *n_atoms_out = static_cast<int>(atoms.size());
    *status = 0;
    return data_id;
}, status);}

int data_read(
    const char* filename,
    int* status
) {return execute_with_catch([&]() {
    auto dataset = SimpleDataset(std::string(filename));
    auto data_id = api::ObjectStorage::register_object(std::move(dataset));
    return data_id;
}, status);}

struct _data_get_data_obj {
    _data_get_data_obj(unsigned int size) :
        q(size), I(size), Ierr(size)
    {}
    std::vector<double> q, I, Ierr;
};
int data_get_data(
    int object_id,
    double** q, double** I, double** Ierr, int* n_points,
    int* status
) {return execute_with_catch([&]() {
    auto dataset = api::ObjectStorage::get_object<SimpleDataset>(object_id);
    if (!dataset) {ErrorMessage::last_error = "Invalid dataset id: \"" + std::to_string(object_id) + "\""; return -1;}
    _data_get_data_obj data(dataset->size());
    for (unsigned int i = 0; i < dataset->size(); ++i) {
        data.q[i] = dataset->x(i);
        data.I[i] = dataset->y(i);
        data.Ierr[i] = dataset->yerr(i);
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_data_get_data_obj>(data_id);
    *q = ref->q.data();
    *I = ref->I.data();
    *Ierr = ref->Ierr.data();
    *n_points = static_cast<int>(dataset->size());
    return data_id;
}, status);}

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
    const char* hydration_model,
    int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return;}
    settings::detail::parse_option("hydration_model", {std::string(hydration_model)});
    molecule->generate_new_hydration();
}, status);}

struct _molecule_distance_histogram_obj {
    _molecule_distance_histogram_obj(unsigned int n_bins) : aa(n_bins), aw(n_bins), ww(n_bins) {}
    std::vector<double> aa, aw, ww;
};
int molecule_distance_histogram(
    int molecule_id,
    double** aa, double** aw, double** ww,
    int* n_bins, double* delta_r, int* status
) {return execute_with_catch([&]() {
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {ErrorMessage::last_error = "Invalid molecule id: \"" + std::to_string(molecule_id) + "\""; return -1;}
    auto hist = molecule->get_histogram();
    _molecule_distance_histogram_obj data(constants::axes::q_axis.bins);
    data.aa = hist->get_aa_counts();
    data.aw = hist->get_aw_counts();
    data.ww = hist->get_ww_counts();
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_molecule_distance_histogram_obj>(data_id);
    *aa = ref->aa.data();
    *aw = ref->aw.data();
    *ww = ref->ww.data();
    *n_bins = static_cast<int>(constants::axes::q_axis.bins);
    *delta_r = constants::axes::d_axis.step();
    return data_id;
}, status);}

struct _molecule_debye_obj {
    _molecule_debye_obj(unsigned int size) :
        q(size), I(size)
    {}
    std::vector<double> q, I;
};
int molecule_debye(
    int molecule_id, 
    const char* exv_model, double** q, double** I, int* n_points,
    int* status
) {return execute_with_catch([&]() {
    settings::detail::parse_option("exv_model", {std::string(exv_model)});
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    molecule->reset_histogram_manager();
    if (!molecule) {*status = 2; return -1;}
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
    const char* exv_model, double* q, double* I, int n_points,
    int* status
) {return execute_with_catch([&]() {
    settings::detail::parse_option("exv_model", {std::string(exv_model)});
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    molecule->reset_histogram_manager();
    if (!molecule) {*status = 2; return;}
    std::vector<double> q_vals(q, q + n_points);
    auto hist = molecule->get_histogram();
    auto debye_I = hist->debye_transform(q_vals);
    if (static_cast<int>(debye_I.size()) != n_points) {*status = 3; return;}
    for (int i = 0; i < n_points; ++i) {
        I[i] = debye_I.y(i);
    }
}, status);}

int molecule_debye_fit(
    int molecule_id, int data_id,
    const char* exv_model,
    int* status
) {return execute_with_catch([&]() {
    settings::detail::parse_option("exv_model", {std::string(exv_model)});
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

int pdb_debye_fit(
    int pdb_id, int data_id,
    const char* exv_model,
    int* status
) {return execute_with_catch([&]() {
    settings::detail::parse_option("exv_model", {std::string(exv_model)});
    auto pdb = api::ObjectStorage::get_object<io::pdb::PDBStructure>(pdb_id);
    if (!pdb) {ErrorMessage::last_error = "Invalid pdb id: \"" + std::to_string(pdb_id) + "\""; return -1;}
    if (settings::molecule::implicit_hydrogens) {pdb->add_implicit_hydrogens();}
    auto data = pdb->reduced_representation();
    auto molecule = data.waters.empty() 
        ? Molecule({Body{std::move(data.atoms)}})
        : Molecule({Body{std::move(data.atoms), std::move(data.waters)}})
    ;
    molecule.reset_histogram_manager();
    auto dataset = api::ObjectStorage::get_object<SimpleDataset>(data_id);
    if (!dataset) {ErrorMessage::last_error = "Invalid dataset id: \"" + std::to_string(data_id) + "\""; return -1;}
    auto fitter = fitter::SmartFitter(*dataset, molecule.get_histogram());
    int fit_result_id = api::ObjectStorage::register_object(fitter.fit());
    return fit_result_id;
}, status);}

struct _fit_get_fit_info_obj {
    _fit_get_fit_info_obj(unsigned int n_pars) : 
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
    _fit_get_fit_curves_obj(unsigned int size) :
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

// #include <em/ImageStack.h>
// int map_read(
//     const char* filename,
//     int* status
// ) {return execute_with_catch([&]() {
//     return map_id;
// }, status);}

// void map_get_slice(
//     int map_id,
//     double z_position,
//     double** slice_data,
//     int* width, int* height,
//     int* status
// );

// int map_fit(
//     int map_id, int data_id,
//     int* status
// );