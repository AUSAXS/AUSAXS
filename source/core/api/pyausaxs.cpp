// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs.h>
#include <api/ObjectStorage.h>
#include <io/detail/PDBReader.h>
#include <io/Reader.h>
#include <dataset/SimpleDataset.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <settings/MoleculeSettings.h>

#include <string>

using namespace ausaxs;
using namespace ausaxs::data;

int read_pdb(
    const char* filename,
    int* status
) {
    *status = 1;
    auto pdb = io::detail::pdb::read(ausaxs::io::ExistingFile(std::string(filename)));
    auto pdb_id = api::ObjectStorage::register_object(std::move(pdb));
    *status = 0;
    return pdb_id;
}

void deallocate(int object_id, int* status) {
    *status = 1;
    api::ObjectStorage::unregister_object(object_id);
    *status = 0;
}

int pdb_read(
    const char* filename,
    int* status
) {
    *status = 1;
    auto pdb = io::Reader::read(std::string(filename));
    auto pdb_id = api::ObjectStorage::register_object(std::move(pdb));
    *status = 0;
    return pdb_id;
}

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
) {
    *status = 1;
    auto pdb = api::ObjectStorage::get_object<io::pdb::PDBStructure>(object_id);
    if (!pdb) {*status = 2; return -1;}
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
}

int data_read(
    const char* filename,
    int* status
) {
    *status = 1;
    auto dataset = SimpleDataset(std::string(filename));
    auto data_id = api::ObjectStorage::register_object(std::move(dataset));
    *status = 0;
    return data_id;
}

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
) {
    *status = 1;
    auto dataset = api::ObjectStorage::get_object<SimpleDataset>(object_id);
    if (!dataset) {*status = 2; return -1;}
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
    *status = 0;
    return data_id;
}

int molecule_from_file(const char* filename, int* status) {
    *status = 1;
    auto molecule = data::Molecule(std::string(filename));
    auto molecule_id = api::ObjectStorage::register_object(std::move(molecule));
    *status = 0;
    return molecule_id;
}

int molecule_from_pdb_id(int pdb_id, int* status) {
    *status = 1;
    auto pdb = api::ObjectStorage::get_object<io::pdb::PDBStructure>(pdb_id);
    if (!pdb) {*status = 2; return -1;}
    if (settings::molecule::implicit_hydrogens) {pdb->add_implicit_hydrogens();}
    auto data = pdb->reduced_representation();
    auto molecule = data.waters.empty() 
        ? Molecule({Body{std::move(data.atoms)}})
        : Molecule({Body{std::move(data.atoms), std::move(data.waters)}})
    ;
    auto molecule_id = api::ObjectStorage::register_object(std::move(molecule));
    *status = 0;
    return molecule_id;
}

int molecule_from_arrays(double* xx, double* yy, double* zz, double* ww, int n_atoms, int* status) {
    *status = 1;
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
    *status = 0;
    return molecule_id;
}

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
) {
    *status = 1;
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {*status = 2; return -1;}
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
    *status = 0;
    return data_id;
}

void molecule_hydrate(
    int molecule_id,
    int* status
) {
    *status = 1;
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {*status = 2; return;}
    molecule->generate_new_hydration();
    *status = 0;
}

struct _molecule_distance_histogram_obj {
    _molecule_distance_histogram_obj(unsigned int n_bins) :
        aa(n_bins), aw(n_bins), ww(n_bins), ax(n_bins), xx(n_bins), wx(n_bins)
    {}
    std::vector<double> aa, aw, ww, ax, xx, wx;
};
int molecule_distance_histogram(
    int molecule_id,
    double** aa, double** aw, double** ww,
    double** ax, double** xx, double** wx,
    int* n_bins, double* delta_r, bool* exv_hists, int* status
) {
    *status = 1;
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {*status = 2; return -1;}
    auto hist = molecule->get_histogram();
    _molecule_distance_histogram_obj data(constants::axes::q_axis.bins);
    if (auto cast = dynamic_cast<hist::ICompositeDistanceHistogramExv*>(hist.get())) {
        data.aa = cast->get_profile_aa().get_counts();
        data.aw = cast->get_profile_aw().get_counts();
        data.ww = cast->get_profile_ww().get_counts();
        data.ax = cast->get_profile_ax().get_counts();
        data.xx = cast->get_profile_xx().get_counts();
        data.wx = cast->get_profile_wx().get_counts();
        *exv_hists = true;
    } else {
        data.aa = hist->get_profile_aa().get_counts();
        data.aw = hist->get_profile_aw().get_counts();
        data.ww = hist->get_profile_ww().get_counts();
        *exv_hists = false;
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_molecule_distance_histogram_obj>(data_id);
    *aa = ref->aa.data();
    *aw = ref->aw.data();
    *ww = ref->ww.data();
    *ax = ref->ax.data();
    *xx = ref->xx.data();
    *wx = ref->wx.data();
    *n_bins = static_cast<int>(constants::axes::q_axis.bins);
    *delta_r = constants::axes::q_axis.step();
    *status = 0;
    return data_id;
}

int molecule_debye(
    int molecule_id,
    double** q, double** I,
    int* n_points, int* status
) {
    *status = 1;
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {*status = 2; return -1;}
    auto hist = molecule->get_histogram();
    auto debye_I = hist->debye_transform();
    _data_get_data_obj data(debye_I.size());
    for (unsigned int i = 0; i < debye_I.size(); ++i) {
        data.q[i] = constants::axes::q_vals[i];
        data.I[i] = debye_I[i];
    }
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_data_get_data_obj>(data_id);
    *q = ref->q.data();
    *I = ref->I.data();
    *n_points = static_cast<int>(debye_I.size());
    *status = 0;
    return data_id;
}

void molecule_debye_userq(
    int molecule_id, double* q, int n_points,
    double* I, int* status
) {
    *status = 1;
    auto molecule = api::ObjectStorage::get_object<Molecule>(molecule_id);
    if (!molecule) {*status = 2; return;}
    std::vector<double> q_vals(q, q + n_points);
    auto hist = molecule->get_histogram();
    auto debye_I = hist->debye_transform(q_vals);
    if (static_cast<int>(debye_I.size()) != n_points) {*status = 3; return;}
    for (int i = 0; i < n_points; ++i) {
        I[i] = debye_I.y(i);
    }
    *status = 0;
}