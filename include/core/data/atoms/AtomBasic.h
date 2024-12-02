#pragma once

#include <math/Vector3.h>
#include <constants/ConstantsCoords.h>

#include <type_traits>

namespace ausaxs::data {
    namespace detail {
        template<typename T>
        struct AtomBasicForwarder {
            using precision_t = constants::coords_precision_t;

            template<numeric Q> void translate(const Vector3<Q>& t) {static_cast<T*>(this)->get_atom_basic().coords += t;}

            template<numeric Q> void set_coordinates(Vector3<Q> v) {static_cast<T*>(this)->get_atom_basic().coords = std::move(v);}
            template<numeric Q> void set_position(Vector3<Q> v) {static_cast<T*>(this)->get_atom_basic().coords = std::move(v);}

            [[nodiscard]] const Vector3<precision_t>& get_coordinates() const {return static_cast<T*>(this)->get_atom_basic().coords;}
            [[nodiscard]] Vector3<precision_t>& get_coordinates() {return static_cast<T*>(this)->get_atom_basic().coords;}
            [[nodiscard]] const Vector3<precision_t>& get_position() const {return static_cast<T*>(this)->get_atom_basic().coords;}
            [[nodiscard]] Vector3<precision_t>& get_position() {return static_cast<T*>(this)->get_atom_basic().coords;}

            [[nodiscard]] precision_t x() const {return static_cast<T*>(this)->get_atom_basic().coords.x();}
            [[nodiscard]] precision_t y() const {return static_cast<T*>(this)->get_atom_basic().coords.y();}
            [[nodiscard]] precision_t z() const {return static_cast<T*>(this)->get_atom_basic().coords.z();}
            [[nodiscard]] precision_t& x() {return static_cast<T*>(this)->get_atom_basic().coords.x();}
            [[nodiscard]] precision_t& y() {return static_cast<T*>(this)->get_atom_basic().coords.y();}
            [[nodiscard]] precision_t& z() {return static_cast<T*>(this)->get_atom_basic().coords.z();}
        };
    }

    /**
     * @brief The most basic information of an atom that is needed to calculate a distance histogram.
     */
    struct AtomBasic : detail::AtomBasicForwarder<AtomBasic> {
        using precision_t = constants::coords_precision_t;

        Vector3<precision_t> coords;
        precision_t weight;
    };
    static_assert(sizeof(AtomBasic) == 4*sizeof(constants::coords_precision_t), "AtomBasic size is off");
    static_assert(std::is_trivial_v<AtomBasic>,                                 "AtomBasic is not trivial");
    static_assert(std::is_standard_layout_v<AtomBasic>,                         "AtomBasic is not standard layout");
    static_assert(supports_nothrow_move_v<AtomBasic>,                           "AtomBasic should support nothrow move semantics.");

    /**
     * @brief The most basic information of an atom that is needed to calculate a distance histogram.
     */
    struct WaterBasic : detail::AtomBasicForwarder<WaterBasic> {
        using precision_t = constants::coords_precision_t;

        Vector3<precision_t> coords;
        precision_t weight;
    };
    static_assert(sizeof(WaterBasic) == 4*sizeof(constants::coords_precision_t), "WaterBasic size is off");
    static_assert(std::is_trivial_v<WaterBasic>,                                 "WaterBasic is not trivial");
    static_assert(std::is_standard_layout_v<WaterBasic>,                         "WaterBasic is not standard layout");
    static_assert(supports_nothrow_move_v<WaterBasic>,                           "WaterBasic should support nothrow move semantics.");
}