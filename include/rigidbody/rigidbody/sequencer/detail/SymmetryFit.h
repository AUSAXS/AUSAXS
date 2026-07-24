// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/symmetry/ISymmetry.h>
#include <math/Vector3.h>

#include <memory>
#include <vector>

namespace ausaxs::rigidbody::sequencer::detail {
    /**
     * @brief The outcome of fitting a symmetry to an existing, assembled structure.
     */
    struct SymmetryFitResult {
        std::unique_ptr<symmetry::ISymmetry> symmetry;  //< the target symmetry with parameters fitted to the input
        double rmsd = 0;                                //< residual RMSD between the input copies and the reconstructed copies
    };

    /**
     * @brief Fit a target symmetry's parameters to an assembled structure.
     *
     * The inverse of expanding one body into many: recover the symmetry parameters that best reproduce
     * @p copies as `copies[0] + template_symmetry`. Supports point, cyclic, polyhedral (dihedral,
     * planar-dihedral, tetrahedral, octahedral, icosahedral) symmetries and composites of these.
     *
     * @param template_symmetry The target symmetry; only its type and repetition count are used, its
     *        parameters are overwritten in the result.
     * @param reference_cm The centre the symmetry is generated about (the reference body's cm).
     * @param copies Atom coordinates of each body: copies[0] is the reference, copies[1..repetitions]
     *        the symmetric copies in order. All must be equal in size, with copies.size() == repetitions()+1.
     */
    SymmetryFitResult fit_symmetry(
        const symmetry::ISymmetry& template_symmetry,
        const Vector3<double>& reference_cm,
        const std::vector<std::vector<Vector3<double>>>& copies
    );

    /**
     * @brief Like @ref fit_symmetry, but for copies given in an unknown order (e.g. PDB chains).
     *
     * Tries the given order and, if the residual exceeds @p accept_rmsd, searches over orderings of the
     * non-reference copies for the best fit. The search is skipped for assemblies too large to enumerate.
     *
     * @param accept_rmsd Residual RMSD (Å) at or below which an ordering is accepted immediately.
     */
    SymmetryFitResult fit_symmetry_best_order(
        const symmetry::ISymmetry& template_symmetry,
        const Vector3<double>& reference_cm,
        const std::vector<std::vector<Vector3<double>>>& copies,
        double accept_rmsd
    );

    /**
     * @brief Expand a reference body into the full set of copies implied by a symmetry.
     *
     * Returns repetitions()+1 point sets: [0] is @p reference, [k] its k-th symmetry image about
     * @p reference_cm. Shared by the fit residual and the API's decomposed structure. This is the
     * coordinate-level counterpart of BodySymmetryFacade::explicit_structure (which works on a Body).
     */
    std::vector<std::vector<Vector3<double>>> reconstruct_copies(
        const symmetry::ISymmetry& symmetry,
        const Vector3<double>& reference_cm,
        const std::vector<Vector3<double>>& reference
    );
}
