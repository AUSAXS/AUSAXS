#pragma once

#include <data/DataFwd.h>
#include <hist/detail/HistDetailFwd.h>
#include <utility/observer_ptr.h>

namespace ausaxs::hist::detail {
    /**
     * @brief This class exists as a non-templated way to disable the effective charge excluded volume approximation used in the simple histogram managers.
     */
    class SimpleExvModel {
        public:
            SimpleExvModel(observer_ptr<const data::Molecule> protein);
            virtual ~SimpleExvModel();

            /**
             * @brief Enable the effective charge excluded volume model, thus subtracting the average excluded volume charge from each atom.
             */
            static void enable();

            /**
             * @brief Disable the effective charge excluded volume model.
             */
            static void disable();

        protected:
            /**
             * @brief Account for the excluded volume in the data.
             *		  Note: this should not be done for models with explicit excluded volume terms.
             *
             * This is done by subtracting the average excluded volume charge from each atom.
             *
             * @param data_a The atomic data to apply the excluded volume transformation to.
             */
            virtual void apply_simple_excluded_volume(hist::detail::CompactCoordinates& data_a);
        
        private:
            observer_ptr<const data::Molecule> protein;
    };
}