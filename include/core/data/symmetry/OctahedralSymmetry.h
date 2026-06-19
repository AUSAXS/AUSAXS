// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/IPolyhedralSymmetry.h>

namespace ausaxs::symmetry {
    class OctahedralSymmetry final : public IPolyhedralSymmetry {
        public:
            std::unique_ptr<ISymmetry> clone() const override;

        private:
            const GroupData& group() const override;
    };
}
