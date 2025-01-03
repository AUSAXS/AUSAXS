#pragma once

#include <data/DataFwd.h>
#include <data/Symmetry.h>
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
            void add(const detail::Symmetry& symmetry) requires (!CONST);
            void add(detail::Symmetry&& symmetry) requires (!CONST); //< @copydoc add_symmetry()

            /**
             * @brief Get the symmetries of this body.
             */
            std::vector<detail::Symmetry>& get() requires (!CONST);
            const std::vector<detail::Symmetry>& get() const; //< @copydoc get_symmetries()

            /**
             * @brief Get the symmetry at the specified index.
             */
            detail::Symmetry& get(unsigned int index) requires (!CONST);
            const detail::Symmetry& get(unsigned int index) const; //< @copydoc get_symmetry()

            /**
             * @brief Write the Body and all its symmetries to a file.
             */
            void save(const io::File& path) const;

        private:
            observer_ptr<BODY> body;
    };
}

template<typename BODY, bool CONST>
void ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::add(const detail::Symmetry& symmetry) requires (!CONST) {
    body->symmetries.emplace_back(symmetry);
    body->changed_internal_state();
    body->changed_external_state();
}

template<typename BODY, bool CONST>
void ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::add(detail::Symmetry&& symmetry) requires (!CONST) {
    body->symmetries.emplace_back(std::move(symmetry));
    body->changed_internal_state();
    body->changed_external_state();
}

template<typename BODY, bool CONST>
std::vector<ausaxs::data::detail::Symmetry>& ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::get() requires (!CONST) {
    return body->symmetries;
}

template<typename BODY, bool CONST>
const std::vector<ausaxs::data::detail::Symmetry>& ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::get() const {
    return body->symmetries;
}

template<typename BODY, bool CONST>
ausaxs::data::detail::Symmetry& ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::get(unsigned int index) requires (!CONST) {
    return body->symmetries[index];
}

template<typename BODY, bool CONST>
const ausaxs::data::detail::Symmetry& ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::get(unsigned int index) const {
    return body->symmetries[index];
}

template<typename BODY, bool CONST>
void ausaxs::data::detail::BodySymmetryFacade<BODY, CONST>::save(const io::File& path) const {
}