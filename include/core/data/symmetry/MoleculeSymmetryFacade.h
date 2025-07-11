// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/observer_ptr.h>
#include <data/atoms/AtomFF.h>
#include <data/detail/SimpleBody.h>
#include <data/DataFwd.h>

namespace ausaxs::symmetry::detail {
    class MoleculeSymmetryFacade {
        public:
            MoleculeSymmetryFacade(observer_ptr<const data::Molecule> molecule) : molecule(molecule) {}

            /**
             * @brief Apply all symmetries to the molecule and get all atoms in the resulting explicit structure.
             */
            data::detail::SimpleBody explicit_structure() const;

            /**
             * @brief Check if the molecule has any symmetries defined.
             */
            bool has_symmetries() const;

            /**
             * @brief Save the explicit structure to a file.
             */
            void save(const io::File& path) const;

        private:
            observer_ptr<const data::Molecule> molecule;
    };
}