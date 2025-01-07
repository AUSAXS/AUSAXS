/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/manager/SimpleProteinManager.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <data/Molecule.h>
#include <utility/Exceptions.h>
#include <data/Body.h>

using namespace ausaxs;

void em::managers::SimpleProteinManager::update_protein(double cutoff) {
    auto atoms = generate_atoms(cutoff);
    std::vector<data::AtomFF> converted(atoms.size());
    std::transform(atoms.begin(), atoms.end(), converted.begin(), [] (const data::EMAtomFF& atom) {return atom.get_atom_ff();});
    protein = std::make_unique<data::Molecule>(std::vector{data::Body{converted}});
    protein->set_histogram_manager(std::make_unique<hist::HistogramManagerMT<true>>(protein.get()));
}