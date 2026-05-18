// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/ISymmetry.h>
#include <math/Matrix.h>
#include <math/Vector3.h>

#include <vector>

namespace ausaxs::symmetry {
    /**
     * @brief A rigid affine map  v -> rotation*v + translation.
     */
    struct AffineTransform {
        Matrix<double> rotation = Matrix<double>::identity(3);
        Vector3<double> translation{0, 0, 0};
    };

    /**
     * @brief Group the C(n,2) unordered pairs among n placements into classes with
     *        geometrically identical inter-placement distances.
     *
     * placements[0] is the original body; placements[1..n-1] are the generated copies.
     * Two pairs are equivalent iff their relative rigid transform is identical (up to
     * inversion, since distance is symmetric). One representative CopyPair is returned per
     * class, with scale set to the class size. Used by the polyhedral, composite, and
     * reference symmetries whose copies do not form a simple cyclic chain.
     */
    std::vector<CopyPair> compute_pair_schedule(const std::vector<AffineTransform>& placements);
}
