#include <em/manager/SimpleProteinManager.h>
#include <hist/HistogramManagerMT.h>
#include <data/Protein.h>
#include <utility/Exceptions.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <data/Body.h>

void em::managers::SimpleProteinManager::update_protein(double cutoff) {
    protein = std::make_unique<Protein>(generate_atoms(cutoff));
    // throw except::disabled("em::managers::SimpleProteinManager::update_protein: This function is disabled and will be removed.");
    // protein->set_histogram_manager<hist::HistogramManagerMT>();
}