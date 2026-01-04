// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/api_pyausaxs.h>
#include <api/ObjectStorage.h>
#include <io/Reader.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <fitter/SmartFitter.h>
#include <settings/All.h>

#include <string>

using namespace ausaxs;
using namespace ausaxs::data;

int pdb_read(
    const char* filename,
    int* status
) {return execute_with_catch([&]() {
    auto pdb = io::Reader::read(std::string(filename));
    auto pdb_id = api::ObjectStorage::register_object(std::move(pdb));
    return pdb_id;
}, status);}

struct _pdb_get_data_obj {
    explicit _pdb_get_data_obj(unsigned int size) :
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

int pdb_debye_fit(
    int pdb_id, int data_id,
    int* status
) {return execute_with_catch([&]() {
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