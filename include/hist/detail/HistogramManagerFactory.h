#pragma once

#include <memory>

class Protein;
namespace settings::hist {enum class HistogramManagerChoice;}
namespace hist {
    class HistogramManager;
    namespace factory {
        /**
         * @brief Create a HistogramManager object.
         */
        std::unique_ptr<HistogramManager> construct_histogram_manager(Protein* protein);

        /**
         * @brief Create a HistogramManager object.
         */
        std::unique_ptr<HistogramManager> construct_histogram_manager(Protein* protein, const settings::hist::HistogramManagerChoice& choice);
    }
}