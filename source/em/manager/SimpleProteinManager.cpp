/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <em/manager/SimpleProteinManager.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <data/Molecule.h>
#include <utility/Exceptions.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>

using namespace ausaxs;

void em::managers::SimpleProteinManager::update_protein(double cutoff) {
    protein = std::make_unique<data::Molecule>(generate_atoms(cutoff));
    protein->set_histogram_manager(std::make_unique<hist::HistogramManagerMT<true>>(protein.get()));
}