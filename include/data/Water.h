#pragma once

#include <data/Atom.h>

class Water : public Atom {
    public:
        using Atom::Atom; // inherit constructors from Atom
        Water(const Atom&& a) noexcept;
        Water(const Atom& a);
        Water(const Water& w);
        ~Water() override;

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
        static Water create_new_water(Vector3<double> coords);

        Water& operator=(const Water& rhs);

        bool operator==(const Water& rhs) const;
};