/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <hist/histogram_manager/HistogramManagerFactory.h>
#include <hist/histogram_manager/HistogramManager.h>
#include <hist/histogram_manager/HistogramManagerMT.h>
#include <hist/histogram_manager/HistogramManagerMTFFAvg.h>
#include <hist/histogram_manager/HistogramManagerMTFFExplicit.h>
#include <hist/histogram_manager/HistogramManagerMTFFGrid.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridSurface.h>
#include <hist/histogram_manager/HistogramManagerMTFFGridScalableExv.h>
#include <hist/histogram_manager/SymmetryManagerMT.h>
#include <hist/histogram_manager/PartialHistogramManager.h>
#include <hist/histogram_manager/PartialHistogramManagerMT.h>
#include <hist/histogram_manager/PartialSymmetryManagerMT.h>
#include <data/Molecule.h>
#include <utility/Exceptions.h>
#include <utility/Console.h>
#include <utility/Logging.h>
#include <settings/HistogramSettings.h>
#include <settings/GeneralSettings.h>
#include <settings/ExvSettings.h>
#include <settings/FitSettings.h>
#include <settings/Flags.h>

using namespace ausaxs;
using namespace ausaxs::hist::factory;

std::unique_ptr<hist::IHistogramManager> hist::factory::construct_histogram_manager(
    observer_ptr<const data::Molecule> protein, bool use_weighted_distribution
) {
    auto choice = settings::hist::get_histogram_manager();
    bool has_syms = protein->symmetry().has_symmetries();
    if (has_syms) {
        switch (choice) {
            case settings::hist::HistogramManagerChoice::HistogramManager:
            case settings::hist::HistogramManagerChoice::HistogramManagerMT:
                choice = settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT;
                break;
            default: 
                console::print_warning(
                    "construct_histogram_manager: Molecule contains symmetries, but the chosen excluded volume method does not support them. "
                    "Symmetries will be ignored. "
                );
                break;
        }
    }
    return construct_histogram_manager(protein, choice, use_weighted_distribution);
}

std::unique_ptr<hist::IHistogramManager> hist::factory::construct_histogram_manager(
    observer_ptr<const data::Molecule> protein, settings::hist::HistogramManagerChoice choice, bool use_weighted_distribution
) {
    if (use_weighted_distribution) {
        switch (choice) {
            case settings::hist::HistogramManagerChoice::HistogramManager: 
                return std::make_unique<HistogramManager<true>>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMT:
                return std::make_unique<HistogramManagerMT<true>>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
                return std::make_unique<HistogramManagerMTFFAvg<true>>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit:
                return std::make_unique<HistogramManagerMTFFExplicit<true>>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid: 
                return std::make_unique<HistogramManagerMTFFGrid>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface: 
                return std::make_unique<HistogramManagerMTFFGridSurface>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridScalableExv: 
                return std::make_unique<HistogramManagerMTFFGridScalableExv>(protein);

            case settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT:
                return std::make_unique<SymmetryManagerMT<true>>(protein);

            case settings::hist::HistogramManagerChoice::PartialHistogramManager:
                return std::make_unique<PartialHistogramManager<true>>(protein);

            case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
                return std::make_unique<PartialHistogramManagerMT<true>>(protein);

            case settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT:
                return std::make_unique<PartialSymmetryManagerMT<true>>(protein);

            // case settings::hist::HistogramManagerChoice::DebugManager:
            //     return std::make_unique<DebugManager<true>>(protein);

            case settings::hist::HistogramManagerChoice::FoXSManager:
            case settings::hist::HistogramManagerChoice::PepsiManager:
            case settings::hist::HistogramManagerChoice::CrysolManager:
                // FoXSManager, PepsiManager, and CrysolManager are all extensions of the HistogramManagerMTFFExplicit method
                return std::make_unique<HistogramManagerMTFFExplicit<true>>(protein);

            default:
                throw except::unknown_argument("hist::factory::construct_histogram_manager: Unkown HistogramManagerChoice. Did you forget to add it to the switch statement?");
        }
    } else {
        switch (choice) {
            case settings::hist::HistogramManagerChoice::HistogramManager: 
                return std::make_unique<HistogramManager<false>>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMT:
                return std::make_unique<HistogramManagerMT<false>>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFAvg:
                return std::make_unique<HistogramManagerMTFFAvg<false>>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFExplicit:
                return std::make_unique<HistogramManagerMTFFExplicit<false>>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGrid:
                return std::make_unique<HistogramManagerMTFFGrid>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridSurface:
                return std::make_unique<HistogramManagerMTFFGridSurface>(protein);

            case settings::hist::HistogramManagerChoice::HistogramManagerMTFFGridScalableExv:
                return std::make_unique<HistogramManagerMTFFGridScalableExv>(protein);

            case settings::hist::HistogramManagerChoice::HistogramSymmetryManagerMT:
                return std::make_unique<SymmetryManagerMT<false>>(protein);

            case settings::hist::HistogramManagerChoice::PartialHistogramManager:
                return std::make_unique<PartialHistogramManager<false>>(protein);

            case settings::hist::HistogramManagerChoice::PartialHistogramManagerMT:
                return std::make_unique<PartialHistogramManagerMT<false>>(protein);

            case settings::hist::HistogramManagerChoice::PartialHistogramSymmetryManagerMT:
                return std::make_unique<PartialSymmetryManagerMT<false>>(protein);

            case settings::hist::HistogramManagerChoice::FoXSManager:
            case settings::hist::HistogramManagerChoice::PepsiManager:
            case settings::hist::HistogramManagerChoice::CrysolManager:
                // FoXSManager, PepsiManager, and CrysolManager are all extensions of the HistogramManagerMTFFExplicit method
                return std::make_unique<HistogramManagerMTFFExplicit<false>>(protein);
            // case settings::hist::HistogramManagerChoice::DebugManager:
            //     return std::make_unique<DebugManager<false>>(protein);

            default:
                throw except::unknown_argument("hist::factory::construct_histogram_manager: Unkown HistogramManagerChoice. Did you forget to add it to the switch statement?");
        }
    }
}