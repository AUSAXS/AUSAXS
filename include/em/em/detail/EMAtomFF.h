#pragma once

#include <data/atoms/AtomFF.h>
#include <form_factor/FormFactorType.h>

namespace ausaxs::data {
    /**
     * @brief The most basic atomic information extended with form factor information.
     */
    class EMAtomFF : public detail::AtomForwarder<EMAtomFF> {
        public:
            EMAtomFF() = default;
            EMAtomFF(const AtomFF& a, double density) : basic(a), density(density) {}
            EMAtomFF(const Vector3<precision_t>& coords, precision_t weight, form_factor::form_factor_t t, double density) : basic(coords, weight, t), density(density) {}
            Atom& get_atom() {return basic.get_atom();}
            const Atom& get_atom() const {return basic.get_atom();}
            AtomFF& get_atom_ff() {return basic;}
            const AtomFF& get_atom_ff() const {return basic;}

            [[nodiscard]] form_factor::form_factor_t form_factor_type() const {return basic.form_factor_type();}
            [[nodiscard]] form_factor::form_factor_t& form_factor_type() {return basic.form_factor_type();}
            [[nodiscard]] double get_density() const {return density;}
            [[nodiscard]] double& get_density() {return density;}

        private:
            AtomFF basic;
            double density;
    };
    // static_assert(sizeof(AtomFF) == sizeof(Atom) + sizeof(form_factor::form_factor_t), "AtomFF should have no padding.");
    static_assert(std::is_trivial_v<EMAtomFF>,         "EMAtomFF is not trivial");
    static_assert(std::is_standard_layout_v<EMAtomFF>, "EMAtomFF is not standard layout");
    static_assert(supports_nothrow_move_v<EMAtomFF>,   "EMAtomFF should support nothrow move semantics.");
}