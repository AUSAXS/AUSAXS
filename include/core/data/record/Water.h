#pragma once

#include <data/record/Atom.h>

namespace data::record {
    class Water : public Atom {
        public:
            using Atom::Atom;
            Water(Atom&& a) noexcept;
            Water(const Atom& a);

            RecordType get_type() const override;

            std::string get_recName() const override;

            bool is_water() const override;

            /**
             * @brief Create a new default water atom.
             */
            static Water create_new_water();

            /**
             * @brief Create a new water atom.
             * @param coords the coordinates for the new atom.
             */
            static Water create_new_water(const Vector3<double>& coords);

            bool operator==(const Water& rhs) const;
    };
}