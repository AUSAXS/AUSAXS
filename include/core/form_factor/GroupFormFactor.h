// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorConcepts.h>
#include <form_factor/GroupGeometry.h>
#include <math/ConstexprMath.h>

namespace ausaxs::form_factor {
    /**
     * @brief The form factor of a heavy atom X with n equivalent bound hydrogens, built from the form factors of its atoms.
     *        A lone atom is a group without hydrogens, whose form factor is that of the atom.
     *
     * @tparam Atom The form factor of a single atom for the probe.
     */
    template<FormFactorType Atom>
    class GroupFormFactor {
        public:
            constexpr GroupFormFactor(const Atom& heavy, const Atom& hydrogen, const group::Geometry& geometry)
                : heavy(heavy), hydrogen(hydrogen), geometry(geometry)
            {}

            /**
             * @brief Evaluate the amplitude f(q), averaged over the orientations of the hydrogens about the heavy atom.
             *        It enters every pair of distinct groups as f_i(q)*f_j(q)*sinc(q*r_ij).
             */
            constexpr double evaluate(double q) const {
                return heavy.evaluate(q) + geometry.hydrogens*hydrogen.evaluate(q)*math::sinc(q*geometry.d_xh);
            }

            /**
             * @brief Evaluate the self-term s(q), which replaces f(q)^2 for a group paired with itself.
             *        The atoms of a group are fixed relative to each other, so this is the Debye sum over the atom pairs within the group.
             *        It only equals f(q)^2 for a lone atom.
             */
            constexpr double evaluate_self(double q) const {
                const double f_heavy = heavy.evaluate(q);
                const double f_hydrogen = hydrogen.evaluate(q);
                const auto& [n, d_xh, d_hh] = geometry;
                return f_heavy*f_heavy + n*f_hydrogen*f_hydrogen
                    + 2*n*f_heavy*f_hydrogen*math::sinc(q*d_xh)
                    + n*(n-1)*f_hydrogen*f_hydrogen*math::sinc(q*d_hh);
            }

            /**
             * @brief Evaluate the amplitude at q = 0.
             */
            constexpr double I0() const {
                return heavy.I0() + geometry.hydrogens*hydrogen.I0();
            }

        private:
            Atom heavy;
            Atom hydrogen;
            group::Geometry geometry;
    };
}
