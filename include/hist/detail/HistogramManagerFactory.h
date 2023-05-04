#pragma once

#include <settings/HistogramSettings.h>

#include <memory>

class Protein;
namespace hist {
    class HistogramManager;
    namespace factory {
        /**
         * @brief Create a HistogramManager object.
         */
        std::unique_ptr<HistogramManager> construct_histogram_manager(Protein* protein, settings::hist::HistogramManagerChoice choice = settings::hist::histogram_manager);
    }
}