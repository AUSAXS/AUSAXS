#pragma once

#include <data/DataFwd.h>
#include <data/detail/SimpleBody.h>
#include <data/symmetry/Symmetry.h>
#include <data/symmetry/SymmetryStorage.h>
#include <data/symmetry/PredefinedSymmetries.h>
#include <utility/observer_ptr.h>
#include <io/IOFwd.h>

namespace ausaxs::symmetry::detail {
    template<typename BODY, bool NONCONST = !std::is_const_v<BODY>>
    class BodySymmetryFacade {
        public:
            BodySymmetryFacade(observer_ptr<BODY> body) : body(body) {}

            /**
             * @brief Add a symmetry to this body.
             */
            void add(symmetry::Symmetry&& symmetry) requires (NONCONST); //< @copydoc add_symmetry()
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