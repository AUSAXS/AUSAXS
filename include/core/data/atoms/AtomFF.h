#pragma once

#include <data/atoms/Atom.h>
#include <form_factor/FormFactorType.h>

namespace ausaxs::data {
    /**
     * @brief The most basic atomic information extended with form factor information.
     */
    class AtomFF : public detail::AtomForwarder<AtomFF> {
        public:
            AtomFF() = default;
            AtomFF(const Atom& a, form_factor::form_factor_t t) : basic(a), type(t) {}
            AtomFF(const Vector3<precision_t>& coords, precision_t weight, form_factor::form_factor_t t) : basic(coords, weight), type(t) {}
            Atom& get_atom() {return basic;}
            const Atom& get_atom() const {return basic;}

            form_factor::form_factor_t get_form_factor_type() const {return type;}
            void set_form_factor_type(form_factor::form_factor_t t) {type = t;}

        private:
            Atom basic;
            form_factor::form_factor_t type;
    };
    // static_assert(sizeof(AtomFF) == sizeof(Atom) + sizeof(form_factor::form_factor_t), "AtomFF should have no padding.");
    static_assert(std::is_trivial_v<AtomFF>,         "AtomFF is not trivial");
    static_assert(std::is_standard_layout_v<AtomFF>, "AtomFF is not standard layout");
    static_assert(supports_nothrow_move_v<AtomFF>,   "AtomFF should support nothrow move semantics.");
}