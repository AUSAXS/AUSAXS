#include <em/manager/SimpleProteinManager.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <data/Molecule.h>
#include <utility/Exceptions.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <data/Body.h>

void em::managers::SimpleProteinManager::update_protein(double cutoff) {
    protein = std::make_unique<data::Molecule>(generate_atoms(cutoff));
    // throw except::disabled("em::managers::SimpleProteinManager::update_protein: This function is disabled and will be removed.");
    // protein->set_histogram_manager<hist::HistogramManagerMT>();
}