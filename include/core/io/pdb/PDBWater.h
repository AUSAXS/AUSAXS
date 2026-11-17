// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/pdb/PDBAtom.h>

namespace ausaxs::io::pdb {
    class PDBWater : public PDBAtom {
        public:
            using PDBAtom::PDBAtom;
            PDBWater(PDBAtom&& a) noexcept;
            PDBWater(const PDBAtom& a);

            RecordType get_type() const override;

            std::string get_recName() const override;

            bool is_water() const override;

            /**
             * @brief Create a new default water atom.
             */
            static PDBWater create_new_water();

            /**
             * @brief Create a new water atom.
             * @param coords the coordinates for the new atom.
             */
            static PDBWater create_new_water(const Vector3<double>& coords);

            bool operator==(const PDBWater& rhs) const;
    };

    /**
     * @brief Map the residue names other programs use for water (SOL, WAT, TIP3, and DOD for heavy water) onto HOH.
     */
    void normalize_water_name(std::string& resName);
}