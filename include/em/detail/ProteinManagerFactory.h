#pragma once

#include <em/ProteinManager.h>
#include <em/SimpleProteinManager.h>

namespace em {
    class ProteinManagerFactory {
        public:
            /**
             * @brief Create a ProteinManager from an ImageStack.
             * 
             * @tparam T: Manager to create. Only the ProteinManager is intended for production. Options:
             *          - SimpleProteinManager: A simple manager that reconstructs the entire protein for each cutoff.
             *          - ProteinManager: A smart manager that only reconstructs the parts of the protein that are needed.
             */
            template<detail::ProteinManagerType T>
            static std::unique_ptr<ProteinManager> create(const ImageStackBase* images) {
                return std::make_unique<T>(images);
            }
    };
}