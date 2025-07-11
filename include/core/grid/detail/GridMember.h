// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <form_factor/FormFactorType.h>
#include <grid/detail/GridInternalFwd.h>
#include <data/atoms/AtomFF.h>
#include <data/atoms/Water.h>
#include <math/Vector3.h>
#include <constants/Constants.h>
#include <type_traits>

namespace ausaxs::grid {
    /**
     * @brief This class is used internally in Grid for storing all information about a particular member atom. 
     */
    template<grid_member_t T>
    class GridMember {
        public:
            GridMember();

            /**
             * @brief Copy constructor.
             */
            GridMember(const GridMember<T>& gm);

            /**
             * @brief Move constructor.
             */
            GridMember(const GridMember<T>&& gm) noexcept;

            /**
             * @brief Constructor.
             * @param atom The atom itself. 
             * @param loc The grid location of the atom. 
             */
            GridMember(const T& atom, Vector3<int> loc);

            /**
             * @brief Constructor.
             * @param atom The atom itself. 
             * @param loc The grid location of the atom. 
             */
            GridMember(T&& atom, Vector3<int> loc);

            ~GridMember();

            /**
             * @brief Get the bin location of this atom.
             */
            Vector3<int>& get_bin_loc();

            /**
             * @brief Get the bin location of this atom.
             */
            const Vector3<int>& get_bin_loc() const;

            /**
             * @brief Get the absolute location of this atom.
             */
            Vector3<double>& get_absolute_loc(); 

            /**
             * @brief Get the absolute location of this atom.
             */
            const Vector3<double>& get_absolute_loc() const;

            /**
             * @brief Check if the volume of this atom has been expanded.
             */
            bool is_expanded() const;

            /**
             * @brief Mark the volume of this atom as expanded or not.
             */
            void set_expanded(bool b);

            /**
             * @brief Get the atom object itself.
             */
            T& get_atom();

            /**
             * @brief Get the atom object itself.
             */
            const T& get_atom() const;

            /**
             * @brief Get the atom type.
             */
            form_factor::form_factor_t get_atom_type() const;

            GridMember& operator=(const GridMember& gm) = default;
            GridMember& operator=(GridMember&& gm) noexcept = default;
            bool operator==(const GridMember& gm) const = default;
            bool operator==(const T& atom) const {return this->atom == atom;}

        private:
            T atom;
            Vector3<int> loc;               // bin location
            bool expanded_volume = false;   // whether the volume of this atom has been expanded
    };

    template<typename T>
    concept valid_gridmember = 
        std::is_same_v<std::remove_cvref_t<T>, GridMember<data::AtomFF>> 
        || std::is_same_v<std::remove_cvref_t<T>, GridMember<data::Water>
    >;
}