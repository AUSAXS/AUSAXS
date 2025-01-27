#pragma once

#include <data/DataFwd.h>
#include <utility/observer_ptr.h>
#include <hist/HistFwd.h>

#include <memory>

namespace ausaxs::settings::hist {enum class HistogramManagerChoice;}
namespace ausaxs::hist {
    namespace factory {
        /**
         * @brief Create a HistogramManager object.
         */
        std::unique_ptr<IHistogramManager> construct_histogram_manager(observer_ptr<const data::Molecule> protein, bool use_weighted_distribution = false);

        /**
         * @brief Create a HistogramManager object.
         */
        std::unique_ptr<IHistogramManager> construct_histogram_manager(observer_ptr<const data::Molecule> protein, const settings::hist::HistogramManagerChoice& choice, bool use_weighted_distribution = false);
    }
}