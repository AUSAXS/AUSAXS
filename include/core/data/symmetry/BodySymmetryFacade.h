// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <data/detail/SimpleBody.h>
#include <data/symmetry/Symmetry.h>
#include <data/symmetry/SymmetryStorage.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <utility/observer_ptr.h>
#include <io/IOFwd.h>

//? GCC and probably also Clang do not like how this class is used as a temporary front for access into the symmetries
//? of a body. Specifically, the dangling reference warning is triggered by the chain 'body->symmetry()->get()'
//? since it believes a reference to a temporary object is being returned. These pragmas silence this warning.
#if defined(__clang__)
    #if __has_warning("-Wdangling")
        #pragma clang diagnostic push
        #pragma clang diagnostic ignored "-Wdangling"
        #define AUSAXS_DANGLING_WARNING
    #endif
#elif defined(__GNUC__) && ( __GNUC__ >= 13)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdangling-reference"
    #define AUSAXS_DANGLING_WARNING
#endif

namespace ausaxs::symmetry::detail {
    template<typename BODY, bool NONCONST = !std::is_const_v<BODY>>
    class BodySymmetryFacade {
        public:
            BodySymmetryFacade(observer_ptr<BODY> body) : body(body) {}

            /**
             * @brief Add a symmetry to this body.
             */
            void add(symmetry::Symmetry&& symmetry) requires (NONCONST);
            void add(symmetry::type symmetry) requires (NONCONST); //< @copydoc add_symmetry()

            /**
             * @brief Get the symmetries of this body.
             *        This will also mark all symmetries as modified. Use the const version to avoid this signal. 
             */
            [[nodiscard]] std::vector<symmetry::Symmetry>& get() requires (NONCONST);

            /**
             * @brief Get the symmetries of this body.
             */
            [[nodiscard]] const std::vector<symmetry::Symmetry>& get() const;

            /**
             * @brief Get the symmetry at the specified index.
             *        This will also mark the symmetry as modified. Use the const version to avoid this signal. 
             */
            [[nodiscard]] symmetry::Symmetry& get(unsigned int index) requires (NONCONST);

            /**
             * @brief Get the symmetry at the specified index.
             */
            [[nodiscard]] const symmetry::Symmetry& get(unsigned int index) const;

            [[nodiscard]] const symmetry::Symmetry& back() const;
            [[nodiscard]] symmetry::Symmetry& back() requires (NONCONST);

            [[nodiscard]] const symmetry::Symmetry& front() const;
            [[nodiscard]] symmetry::Symmetry& front() requires (NONCONST);

            /**
             * @brief Get the symmetry storage object.
             *        This will also mark all symmetries as modified. Use the const version to avoid this signal.
             */
            [[nodiscard]] observer_ptr<symmetry::SymmetryStorage> get_obj() requires (NONCONST);

            /**
             * @brief Get the symmetry storage object.
             */
            [[nodiscard]] observer_ptr<const symmetry::SymmetryStorage> get_obj() const;

            /**
             * @brief Get the total number of atoms in the body, including all symmetries.
             */
            std::size_t size_atom_total() const; 

            /**
             * @brief Get the total number of water molecules in the body, including all symmetries.
             */
            std::size_t size_water_total() const;

            /**
             * @brief Set the symmetry storage object.
             */
            void set_obj(std::unique_ptr<symmetry::SymmetryStorage> obj) requires (NONCONST);

            /**
             * @brief Write the Body and all its symmetries to a file.
             */
            void save(const io::File& path) const;

            /**
             * @brief Apply all symmetries to the body, returning a new larger body with all the symmetries applied.
             */
            [[nodiscard]] data::detail::SimpleBody explicit_structure() const;

        private:
            observer_ptr<BODY> body;
    };
}

#if defined(AUSAXS_DANGLING_WARNING)
    #undef AUSAXS_DANGLING_WARNING
    #if defined(__clang__)
        #pragma clang diagnostic pop
    #elif defined(__GNUC__)
        #pragma GCC diagnostic pop
    #endif
#endif