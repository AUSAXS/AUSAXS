// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/IPolyhedralSymmetry.h>

namespace ausaxs::symmetry {
    /**
     * @brief Chiral dihedral point group D_N (order 2N).
     */
    template<int N>
    class DihedralSymmetry : public IPolyhedralSymmetry {
        static_assert(2 <= N, "DihedralSymmetry: order N must be at least 2.");
        public:
            std::unique_ptr<ISymmetry> clone() const override;
            std::string type_name() const override;

        protected:
            const GroupData& group() const override;
    };
}
