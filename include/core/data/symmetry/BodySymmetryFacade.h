#pragma once

#include <data/DataFwd.h>
#include <data/symmetry/Symmetry.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <utility/observer_ptr.h>
#include <io/IOFwd.h>

namespace ausaxs::data::detail {
    template<typename BODY, bool CONST = std::is_const_v<BODY>>
    class BodySymmetryFacade {
        public:
            BodySymmetryFacade(observer_ptr<BODY> body) : body(body) {}

            /**
             * @brief Add a symmetry to this body.
             */
            void add(const symmetry::Symmetry& symmetry) requires (!CONST);
            void add(symmetry::Symmetry&& symmetry) requires (!CONST); //< @copydoc add_symmetry()
            void add(symmetry::type symmetry) requires (!CONST); //< @copydoc add_symmetry()

            /**
             * @brief Get the symmetries of this body.
             */
            std::vector<symmetry::Symmetry>& get() requires (!CONST);
            const std::vector<symmetry::Symmetry>& get() const; //< @copydoc get_symmetries()

            /**
             * @brief Get the symmetry at the specified index.
             */
            symmetry::Symmetry& get(unsigned int index) requires (!CONST);
            const symmetry::Symmetry& get(unsigned int index) const; //< @copydoc get_symmetry()

            /**
             * @brief Write the Body and all its symmetries to a file.
             */
            void save(const io::File& path) const;

            /**
             * @brief Apply all symmetries to the body, returning a new larger body with all the symmetries applied.
             */
            Body get_explicit_structure() const;

        private:
            observer_ptr<BODY> body;
    };
}