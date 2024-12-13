#pragma once

#include <data/DataFwd.h>
#include <data/Symmetry.h>
#include <io/IOFwd.h>

namespace ausaxs::data::detail {
    template<typename BODY>
    class BodySymmetryFacade {
        public:
            BodySymmetryFacade(BODY& body) : body(body) {}

            /**
             * @brief Add a symmetry to this body.
             */
            template<bool CONST = std::is_const_v<BODY>> requires (!CONST)
            void add_symmetry(const detail::Symmetry& symmetry);

            template<bool CONST = std::is_const_v<BODY>> requires (!CONST)
            void add_symmetry(detail::Symmetry&& symmetry); //< @copydoc add_symmetry()

            /**
             * @brief Get the symmetries of this body.
             */
            template<bool CONST = std::is_const_v<BODY>> requires (!CONST)
            std::vector<detail::Symmetry>& get_symmetries();
            const std::vector<detail::Symmetry>& get_symmetries() const; //< @copydoc get_symmetries()

            /**
             * @brief Get the symmetry at the specified index.
             */
            template<bool CONST = std::is_const_v<BODY>> requires (!CONST)
            detail::Symmetry& get_symmetry(unsigned int index);
            const detail::Symmetry& get_symmetry(unsigned int index) const; //< @copydoc get_symmetry()

            /**
             * @brief Write the Body and all its symmetries to a file.
             */
            void save(const io::File& path) const;

        private:
            BODY& body;
    };
}

// std::vector<data::detail::Symmetry>& Body::get_symmetries() {return symmetries;}
// const std::vector<data::detail::Symmetry>& Body::get_symmetries() const {return symmetries;}

// void Body::add_symmetry(const data::detail::Symmetry& symmetry) {
//     symmetries.push_back(symmetry);
//     changed_internal_state();
//     changed_external_state();
// }

// void Body::add_symmetry(data::detail::Symmetry&& symmetry) {
//     symmetries.push_back(std::move(symmetry));
//     changed_internal_state();
//     changed_external_state();
// }

