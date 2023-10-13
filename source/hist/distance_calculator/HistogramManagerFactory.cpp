#include <hist/distance_calculator/HistogramManagerFactory.h>
#include <settings/HistogramSettings.h>
#include <hist/distance_calculator/HistogramManager.h>
#include <hist/distance_calculator/HistogramManagerMT.h>
#include <hist/distance_calculator/HistogramManagerMTFF.h>
#include <hist/distance_calculator/PartialHistogramManager.h>
#include <hist/distance_calculator/PartialHistogramManagerMT.h>
#include <hist/distance_calculator/DebugManager.h>
#include <data/Protein.h>

using namespace hist::factory;

std::unique_ptr<hist::HistogramManager> hist::factory::construct_histogram_manager(Protein* protein) {
    return hist::factory::construct_histogram_manager(protein, settings::hist::histogram_manager);
}

std::unique_ptr<hist::HistogramManager> hist::factory::construct_histogram_manager(Protein* protein, const settings::hist::HistogramManagerChoice& choice) {
    switch (choice) {
        case settings::hist::HistogramManagerChoice::HistogramManager: 
            return std::make_unique<HistogramManager>(protein);
        case settings::hist::HistogramManagerChoice::HistogramManagerMT:
            return std::make_unique<HistogramManagerMT>(protein);
        case settings::hist::HistogramManagerChoice::HistogramManagerMTFF:
            return std::make_unique<HistogramManagerMTFF>(protein);
        case settings::hist::HistogramManagerChoice::PartialHistogramManager:
            return std::make_unique<PartialHistogramManager>(protein);
        case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
            return std::make_unique<PartialHistogramManagerMT>(protein);
        case settings::hist::HistogramManagerChoice::DebugManager:
            return std::make_unique<DebugManager>(protein);
        default:
            throw except::unknown_argument("hist::factory::construct_histogram_manager: Unkown HistogramManagerChoice. Did you forget to add it to the switch statement?");
    }            
}