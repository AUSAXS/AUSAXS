// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <api/pyausaxs/api_pdb.h>
#include <api/ObjectStorage.h>
#include <io/Reader.h>
#include <io/pdb/PDBStructure.h>
#include <data/Molecule.h>
#include <data/Body.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <rigidbody/sequencer/detail/SymmetryFit.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <fitter/SmartFitter.h>
#include <settings/All.h>

#include <utility/Exceptions.h>
#include <string>
#include <unordered_map>
#include <vector>

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

namespace {
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
}
int pdb_get_data(
    int object_id,
    int** serial_out, const char*** name_out, const char*** altLoc_out, const char*** resName_out, const char** chainID_out, int** resSeq_out, 
    const char*** iCode_out, double** x_out, double** y_out, double** z_out, double** occupancy_out, double** tempFactor_out, const char*** element_out, 
    const char*** charge_out, int* n_atoms_out, int* status
) {return execute_with_catch([&]() {
    auto pdb = api::ObjectStorage::get_object<io::pdb::PDBStructure>(object_id);
    if (!pdb) {throw except::invalid_argument("Invalid pdb id: \"" + std::to_string(object_id) + "\"");}
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

namespace {
struct _pdb_decompose_obj {
    std::vector<double> x, y, z;
    std::vector<int> copy_index;
};
}

int pdb_decompose_symmetry(
    int pdb_id, const char* symmetry_name,
    double** x_out, double** y_out, double** z_out, int** copy_index_out, int* n_atoms_out,
    double* rmsd_out, int* status
) {return execute_with_catch([&]() {
    auto pdb = api::ObjectStorage::get_object<io::pdb::PDBStructure>(pdb_id);
    if (!pdb) {throw except::invalid_argument("Invalid pdb id: \"" + std::to_string(pdb_id) + "\"");}

    // group atoms into chains, preserving first-seen chain order (chain 0 = reference)
    std::unordered_map<char, int> chain_index;
    std::vector<std::vector<Vector3<double>>> chains;
    for (const auto& atom : pdb->atoms) {
        auto [it, inserted] = chain_index.try_emplace(atom.chainID, static_cast<int>(chains.size()));
        if (inserted) {chains.emplace_back();}
        chains[it->second].push_back(atom.coords);
    }
    if (chains.size() < 2) {throw except::invalid_argument("pdb_decompose_symmetry: at least two chains are required.");}
    for (const auto& c : chains) {
        if (c.size() != chains[0].size()) {
            throw except::invalid_argument("pdb_decompose_symmetry: chains have differing atom counts; they must be copies of the same molecule.");
        }
    }

    auto base = symmetry::create(std::string(symmetry_name));
    if (chains.size() != base->repetitions() + 1) {
        throw except::invalid_argument(
            "pdb_decompose_symmetry: symmetry \"" + std::string(symmetry_name) + "\" needs "
            + std::to_string(base->repetitions() + 1) + " chains, but " + std::to_string(chains.size()) + " were found."
        );
    }

    // reference centre = centroid of the first chain
    Vector3<double> cm{0, 0, 0};
    for (const auto& p : chains[0]) {cm += p;}
    cm /= static_cast<double>(chains[0].size());

    // fit (order-independent; accept_rmsd 0 forces the search to keep the globally best ordering)
    auto fit = rigidbody::sequencer::detail::fit_symmetry_best_order(*base, cm, chains, 0.0);

    // build the decomposed structure: reference chain (copy 0) + the fitted symmetry copies
    auto reconstructed = rigidbody::sequencer::detail::reconstruct_copies(*fit.symmetry, cm, chains[0]);
    _pdb_decompose_obj data;
    int per = static_cast<int>(chains[0].size());
    int reps = static_cast<int>(fit.symmetry->repetitions());
    data.x.reserve(per*(reps + 1)); data.y.reserve(per*(reps + 1)); data.z.reserve(per*(reps + 1)); data.copy_index.reserve(per*(reps + 1));
    for (int k = 0; k <= reps; ++k) {
        for (const auto& q : reconstructed[k]) {
            data.x.push_back(q.x());
            data.y.push_back(q.y());
            data.z.push_back(q.z());
            data.copy_index.push_back(k);
        }
    }

    *rmsd_out = fit.rmsd;
    int data_id = api::ObjectStorage::register_object(std::move(data));
    auto ref = api::ObjectStorage::get_object<_pdb_decompose_obj>(data_id);
    *x_out = ref->x.data();
    *y_out = ref->y.data();
    *z_out = ref->z.data();
    *copy_index_out = ref->copy_index.data();
    *n_atoms_out = static_cast<int>(ref->x.size());
    *status = 0;
    return data_id;
}, status);}

int pdb_debye_fit(
    int pdb_id, int data_id,
    int* status
) {return execute_with_catch([&]() {
    auto pdb = api::ObjectStorage::get_object<io::pdb::PDBStructure>(pdb_id);
    if (!pdb) {throw except::invalid_argument("Invalid pdb id: \"" + std::to_string(pdb_id) + "\"");}
    if (settings::molecule::implicit_hydrogens) {pdb->add_implicit_hydrogens();}
    auto data = pdb->reduced_representation();
    auto molecule = data.waters.empty() ? Molecule({Body{std::move(data.atoms)}}) : Molecule({Body{std::move(data.atoms), std::move(data.waters)}});
    molecule.reset_histogram_manager();
    auto dataset = api::ObjectStorage::get_object<SimpleDataset>(data_id);
    if (!dataset) {throw except::invalid_argument("Invalid dataset id: \"" + std::to_string(data_id) + "\"");}
    auto fitter = fitter::SmartFitter(*dataset, molecule.get_histogram());
    int fit_result_id = api::ObjectStorage::register_object(fitter.fit());
    return fit_result_id;
}, status);}