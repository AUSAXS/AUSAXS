#pragma once

#include <math/Vector3.h>
#include <constants/ConstantsCoords.h>

#include <type_traits>

namespace ausaxs::data {
    /**
     * @brief The most basic information of an atom that is needed to calculate a distance histogram.
     */
    struct AtomBasic {
        using precision_t = constants::coords_precision_t;
        template<numeric T> void translate(const Vector3<T>& t) {coords += t;}

        template<numeric T> void set_coordinates(Vector3<T> v) {coords = std::move(v);}
        template<numeric T> void set_position(Vector3<T> v) {coords = std::move(v);}

        [[nodiscard]] const Vector3<precision_t>& get_coordinates() const {return coords;}
        [[nodiscard]] Vector3<precision_t>&       get_coordinates() {return coords;}
        [[nodiscard]] const Vector3<precision_t>& get_position() const {return coords;}
        [[nodiscard]] Vector3<precision_t>&       get_position() {return coords;}

        [[nodiscard]] precision_t x() const {return coords.x();}
        [[nodiscard]] precision_t y() const {return coords.y();}
        [[nodiscard]] precision_t z() const {return coords.z();}
        [[nodiscard]] precision_t& x() {return coords.x();}
        [[nodiscard]] precision_t& y() {return coords.y();}
        [[nodiscard]] precision_t& z() {return coords.z();}

        Vector3<precision_t> coords;
        precision_t weight;
    };

    template<typename T>
    struct AtomBasicForwarder {
        using precision_t = constants::coords_precision_t;

        template<numeric Q> void translate(const Vector3<Q>& t) {static_cast<T*>(this)->get_atom_basic().translate(t);}

        template<numeric Q> void set_coordinates(Vector3<Q> v) {static_cast<T*>(this)->get_atom_basic().set_coordinates(v);}
        template<numeric Q> void set_position(Vector3<Q> v) {static_cast<T*>(this)->get_atom_basic().set_position(v);}

        [[nodiscard]] const Vector3<precision_t>& get_coordinates() const {return static_cast<T*>(this)->get_atom_basic().get_coordinates();}
        [[nodiscard]] Vector3<precision_t>& get_coordinates() {return static_cast<T*>(this)->get_atom_basic().get_coordinates();}
        [[nodiscard]] const Vector3<precision_t>& get_position() const {return static_cast<T*>(this)->get_atom_basic().get_position();}
        [[nodiscard]] Vector3<precision_t>& get_position() {return static_cast<T*>(this)->get_atom_basic().get_position();}

        [[nodiscard]] precision_t x() const {return static_cast<T*>(this)->get_atom_basic().x();}
        [[nodiscard]] precision_t y() const {return static_cast<T*>(this)->get_atom_basic().y();}
        [[nodiscard]] precision_t z() const {return static_cast<T*>(this)->get_atom_basic().z();}
        [[nodiscard]] precision_t& x() {return static_cast<T*>(this)->get_atom_basic().x();}
        [[nodiscard]] precision_t& y() {return static_cast<T*>(this)->get_atom_basic().y();}
        [[nodiscard]] precision_t& z() {return static_cast<T*>(this)->get_atom_basic().z();}
    };

    static_assert(sizeof(AtomBasic) == 4*sizeof(constants::coords_precision_t), "AtomBasic size is off");
    static_assert(std::is_trivial_v<AtomBasic>,                                 "AtomBasic is not trivial");
    static_assert(std::is_standard_layout_v<AtomBasic>,                         "AtomBasic is not standard layout");
    static_assert(supports_nothrow_move_v<AtomBasic>,                           "AtomBasic should support nothrow move semantics.");
}