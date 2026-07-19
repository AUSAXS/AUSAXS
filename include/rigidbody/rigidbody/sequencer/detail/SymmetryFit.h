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
     * @brief Fit the parameters of a target symmetry to an assembled, symmetric structure.
     *
     * The inverse of expanding one body + a symmetry into many: given the N copies of a single
     * molecule that make up @p copies, this recovers the symmetry parameters (offset + frame/axis)
     * that best reproduce the assembly as `copies[0] + template_symmetry`. Because every copy is the
     * same molecule, atom correspondence is exact (copies[k][i] is the image of copies[0][i]), so
     * each per-copy alignment is a closed-form Kabsch superposition rather than a search.
     *
     * @param template_symmetry The target symmetry type. Only its type (and repetition count) is
     *        used; its parameters are ignored and overwritten in the returned object.
     * @param reference_cm The centre of mass of the reference body copies[0] (atoms only), i.e. the
     *        same value the forward transform is generated about.
     * @param copies The atom coordinates of each body. copies[0] is the reference; copies[k]
     *        (k = 1..repetitions) are the symmetric copies in repetition order. All copies must have
     *        the same number of points, and copies.size() must equal repetitions()+1.
     *
     * Supported target types: PointSymmetry (p2), CyclicSymmetry (c2..c12) and all
     * IPolyhedralSymmetry leaves (dihedral, planar dihedral, tetrahedral, octahedral, icosahedral).
     */
    SymmetryFitResult fit_symmetry(
        const symmetry::ISymmetry& template_symmetry,
        const Vector3<double>& reference_cm,
        const std::vector<std::vector<Vector3<double>>>& copies
    );

    /**
     * @brief Fit a symmetry without assuming the copies are given in repetition order.
     *
     * @ref fit_symmetry requires copies[k] to be the k-th symmetry copy, but chains in a PDB come in
     * arbitrary order. This first tries the given order and, if the residual exceeds @p accept_rmsd,
     * searches over orderings of the non-reference bodies (copies[0] is always kept as the reference,
     * which is valid because the symmetry group acts transitively on the copies) and returns the
     * lowest-RMSD fit found. The search is skipped for assemblies too large to enumerate, in which
     * case the given order is used as-is.
     *
     * @param accept_rmsd Residual RMSD (Å) at or below which an ordering is accepted immediately.
     */
    SymmetryFitResult fit_symmetry_best_order(
        const symmetry::ISymmetry& template_symmetry,
        const Vector3<double>& reference_cm,
        const std::vector<std::vector<Vector3<double>>>& copies,
        double accept_rmsd
    );
}
