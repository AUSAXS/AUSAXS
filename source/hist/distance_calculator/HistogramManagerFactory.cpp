#include <hist/distance_calculator/HistogramManagerFactory.h>
#include <hist/distance_calculator/HistogramManager.h>
// #include <hist/distance_calculator/HistogramManagerMT.h>
// #include <hist/distance_calculator/HistogramManagerMTFFAvg.h>
// #include <hist/distance_calculator/HistogramManagerMTFFExplicit.h>
// #include <hist/distance_calculator/HistogramManagerMTFFGrid.h>
// #include <hist/distance_calculator/PartialHistogramManager.h>
// #include <hist/distance_calculator/PartialHistogramManagerMT.h>
// #include <hist/distance_calculator/DebugManager.h>
#include <settings/HistogramSettings.h>
#include <data/Molecule.h>
#include <utility/Exceptions.h>

using namespace hist::factory;

std::unique_ptr<hist::IHistogramManager> hist::factory::construct_histogram_manager(view_ptr<const data::Molecule> protein, bool use_weighted_distribution) {
    return hist::factory::construct_histogram_manager(protein, settings::hist::histogram_manager, use_weighted_distribution);
}

std::unique_ptr<hist::IHistogramManager> hist::factory::construct_histogram_manager(view_ptr<const data::Molecule> protein, const settings::hist::HistogramManagerChoice& choice, bool use_weighted_distribution) {
    if (use_weighted_distribution) {
        switch (choice) {
            case settings::hist::HistogramManagerChoice::HistogramManager: 
                return std::make_unique<HistogramManager<true>>(protein);
            // case settings::hist::HistogramManagerChoice::HistogramManagerMT:
            //     return std::make_unique<HistogramManagerMT<true>>(protein);
            // case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
            //     return std::make_unique<HistogramManagerMTFFAvg<true>>(protein);
            // case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit:
            //     return std::make_unique<HistogramManagerMTFFExplicit<true>>(protein);
            // case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid: 
            //     return std::make_unique<HistogramManagerMTFFGrid<true>>(protein);
            // case settings::hist::HistogramManagerChoice::PartialHistogramManager:
            //     return std::make_unique<PartialHistogramManager<true>>(protein);
            // case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
            //     return std::make_unique<PartialHistogramManagerMT<true>>(protein);
            // case settings::hist::HistogramManagerChoice::DebugManager:
            //     return std::make_unique<DebugManager<true>>(protein);
            default:
                throw except::unknown_argument("hist::factory::construct_histogram_manager: Unkown HistogramManagerChoice. Did you forget to add it to the switch statement?");
        }
    } else {
        switch (choice) {
            case settings::hist::HistogramManagerChoice::HistogramManager: 
                return std::make_unique<HistogramManager<false>>(protein);
            // case settings::hist::HistogramManagerChoice::HistogramManagerMT:
            //     return std::make_unique<HistogramManagerMT<false>>(protein);
            // case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
            //     return std::make_unique<HistogramManagerMTFFAvg<false>>(protein);
            // case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit:
            //     return std::make_unique<HistogramManagerMTFFExplicit<false>>(protein);
            // case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid: 
            //     return std::make_unique<HistogramManagerMTFFGrid<false>>(protein);
            // case settings::hist::HistogramManagerChoice::PartialHistogramManager:
            //     return std::make_unique<PartialHistogramManager<false>>(protein);
            // case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
            //     return std::make_unique<PartialHistogramManagerMT<false>>(protein);
            // case settings::hist::HistogramManagerChoice::DebugManager:
            //     return std::make_unique<DebugManager<false>>(protein);
            default:
                throw except::unknown_argument("hist::factory::construct_histogram_manager: Unkown HistogramManagerChoice. Did you forget to add it to the switch statement?");
        }
    }
}