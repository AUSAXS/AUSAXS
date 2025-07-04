// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/atoms/Atom.h>
#include <data/atoms/AtomFF.h>
#include <form_factor/FormFactorType.h>

namespace ausaxs::data {
    /**
     * @brief The most basic atomic information extended with form factor information.
     */
    class EMAtom : public detail::AtomForwarder<EMAtom> {
        public:
            EMAtom() = default;
            EMAtom(const Atom& a, double density) : basic(a), density(density) {}
            EMAtom(const Vector3<precision_t>& coords, double weight, double density) : basic(coords, weight), density(density) {}
            Atom& get_atom() {return basic.get_atom();}
            const Atom& get_atom() const {return basic.get_atom();}
            AtomFF get_atom_ff() const {return AtomFF{basic, form_factor::form_factor_t::UNKNOWN};}

            [[nodiscard]] double charge_density() const {return density;}
            [[nodiscard]] double& charge_density() {return density;}

        private:
            Atom basic;
            double density;
    };
    // static_assert(sizeof(AtomFF) == sizeof(Atom) + sizeof(form_factor::form_factor_t), "AtomFF should have no padding.");
    static_assert(std::is_trivial_v<EMAtom>,         "EMAtomFF is not trivial");
    static_assert(std::is_standard_layout_v<EMAtom>, "EMAtomFF is not standard layout");
    static_assert(supports_nothrow_move_v<EMAtom>,   "EMAtomFF should support nothrow move semantics.");
}