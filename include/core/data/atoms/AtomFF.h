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
            AtomFF(Atom a, form_factor::form_factor_t t) : basic(std::move(a)), type(t) {}
            AtomFF(Vector3<precision_t> coords, form_factor::form_factor_t t) : basic(coords, constants::charge::nuclear::get_charge(t)), type(t) {}
            AtomFF(Vector3<precision_t> coords, form_factor::form_factor_t t, double weight) : basic(coords, weight), type(t) {}
            Atom& get_atom() {return basic;}
            const Atom& get_atom() const {return basic;}

            [[nodiscard]] form_factor::form_factor_t form_factor_type() const {return type;}
            [[nodiscard]] form_factor::form_factor_t& form_factor_type() {return type;}
            bool operator==(const AtomFF& rhs) const = default;

        private:
            Atom basic;
            form_factor::form_factor_t type;
    };
    // static_assert(sizeof(AtomFF) == sizeof(Atom) + sizeof(form_factor::form_factor_t), "AtomFF should have no padding.");
    static_assert(std::is_trivial_v<AtomFF>,         "AtomFF is not trivial");
    static_assert(std::is_standard_layout_v<AtomFF>, "AtomFF is not standard layout");
    static_assert(supports_nothrow_move_v<AtomFF>,   "AtomFF should support nothrow move semantics.");

	template<typename T> concept AtomVectorFF = std::is_same_v<std::remove_cvref_t<T>, std::vector<data::AtomFF>>;
}