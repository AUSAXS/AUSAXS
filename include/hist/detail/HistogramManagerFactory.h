#pragma once

#include <hist/HistogramManager.h>
#include <hist/HistogramManagerMT.h>
#include <hist/PartialHistogramManager.h>
#include <hist/PartialHistogramManagerMT.h>

namespace hist {
    /**
     * @brief A factory class for creating HistogramManager objects. 
     */
    class HistogramManagerFactory {
        public:
            /**
             * @brief Create a HistogramManager object.
             * 
             * @tparam T: Manager to create. Only the PartialHistogramManagerMT is intended for production. Options:
             *         - HistogramManager: A simple manager that recalculates the entire histogram every time.
             *         - HistogramManagerMT: A multithreaded implementation of the simple manager.
             *         - PartialHistogramManager: A smart manager that only recalculates the parts of the histogram that are needed.
             *         - PartialHistogramManagerMT: A multithreaded implementation of the partial manager.
             */
            template<detail::HistogramManagerType T>
            static std::unique_ptr<HistogramManager> create(Protein* protein) {
                return std::make_unique<T>(protein);
            }
    };
}