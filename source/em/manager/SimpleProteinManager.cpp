#include <em/manager/SimpleProteinManager.h>
#include <hist/HistogramManagerMT.h>
#include <data/Protein.h>
#include <memory>

void em::managers::SimpleProteinManager::update_protein(double cutoff) {
    protein = std::make_shared<Protein>(generate_atoms(cutoff));
    protein->set_histogram_manager<hist::HistogramManagerMT>();
}