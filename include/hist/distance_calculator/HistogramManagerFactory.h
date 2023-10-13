#pragma once

#include <data/DataFwd.h>

#include <memory>

namespace settings::hist {enum class HistogramManagerChoice;}
namespace hist {
    class HistogramManager;
    namespace factory {
        /**
         * @brief Create a HistogramManager object.
         */
        std::unique_ptr<HistogramManager> construct_histogram_manager(data::Molecule* protein);

        /**
         * @brief Create a HistogramManager object.
         */
        std::unique_ptr<HistogramManager> construct_histogram_manager(data::Molecule* protein, const settings::hist::HistogramManagerChoice& choice);
    }
}