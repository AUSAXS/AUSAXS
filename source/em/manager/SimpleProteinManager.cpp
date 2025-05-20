/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/manager/SimpleProteinManager.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <data/Molecule.h>
#include <utility/Exceptions.h>
#include <utility/Logging.h>
#include <data/Body.h>

using namespace ausaxs;

void em::managers::SimpleProteinManager::update_protein(double cutoff) {
    logging::log("SimpleProteinManager::update_protein: cutoff = " + std::to_string(cutoff));
    auto atoms = generate_atoms(cutoff);

    // this sorting step is principially not necessary, but required for consistency with SmartProteinManager
    // otherwise we may see tiny differences in the generated hydration shells which breaks the tests
    std::sort(
        atoms.begin(), 
        atoms.end(), 
        [] (const data::EMAtom& atom1, const data::EMAtom& atom2) {return atom1.charge_density() < atom2.charge_density();}
    );

    std::vector<data::AtomFF> converted(atoms.size());
    std::transform(atoms.begin(), atoms.end(), converted.begin(), [] (const data::EMAtom& atom) {return atom.get_atom_ff();});
    protein = std::make_unique<data::Molecule>(std::vector{data::Body{converted}});
    protein->set_histogram_manager(settings::hist::HistogramManagerChoice::HistogramManagerMT);
}