// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <cstddef>
#include <iterator>
#include <vector>

namespace ausaxs {
    /**
     * @brief Forward iterator that traverses all atoms across all bodies of a molecule in order.
     *        TBody must be complete at instantiation; a forward declaration suffices to name the type.
     *        Use TBody=Body for a mutable iterator, TBody=const Body for a read-only iterator.
     */
    template<typename TBody>
    class MoleculeAtomIterator {
        using body_t      = std::remove_cv_t<TBody>;
        using ref_body_t  = std::conditional_t<std::is_const_v<TBody>, const body_t, body_t>;
        using body_iter_t = std::conditional_t<std::is_const_v<TBody>, typename std::vector<body_t>::const_iterator, typename std::vector<body_t>::iterator>;
        using atom_vec_t  = std::remove_reference_t<decltype(std::declval<ref_body_t&>().get_atoms())>;
        using atom_iter_t = std::conditional_t<std::is_const_v<atom_vec_t>, typename atom_vec_t::const_iterator, typename atom_vec_t::iterator>;

    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = typename std::remove_cv_t<atom_vec_t>::value_type;
        using pointer           = std::conditional_t<std::is_const_v<atom_vec_t>, const value_type*, value_type*>;
        using reference         = std::conditional_t<std::is_const_v<atom_vec_t>, const value_type&, value_type&>;

        MoleculeAtomIterator() = default;

        MoleculeAtomIterator(body_iter_t first, body_iter_t last) : body_cur_(first), body_end_(last) {
            if (body_cur_ != body_end_) {
                atom_cur_ = body_cur_->get_atoms().begin();
                atom_end_ = body_cur_->get_atoms().end();
                skip_empty();
            }
        }

        reference operator*()  const { return *atom_cur_; }
        pointer   operator->() const { return &*atom_cur_; }

        MoleculeAtomIterator& operator++() {
            ++atom_cur_;
            skip_empty();
            return *this;
        }

        MoleculeAtomIterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }

        bool operator==(const MoleculeAtomIterator& rhs) const {
            if (body_cur_ == body_end_ && rhs.body_cur_ == rhs.body_end_) { return true;  }
            if (body_cur_ == body_end_ || rhs.body_cur_ == rhs.body_end_) { return false; }
            return body_cur_ == rhs.body_cur_ && atom_cur_ == rhs.atom_cur_;
        }

        bool operator!=(const MoleculeAtomIterator& rhs) const { return !(*this == rhs); }

    private:
        body_iter_t body_cur_{}, body_end_{};
        atom_iter_t atom_cur_{}, atom_end_{};

        void skip_empty() {
            while (body_cur_ != body_end_ && atom_cur_ == atom_end_) {
                ++body_cur_;
                if (body_cur_ != body_end_) {
                    atom_cur_ = body_cur_->get_atoms().begin();
                    atom_end_ = body_cur_->get_atoms().end();
                }
            }
        }
    };

    template<typename TBody>
    class MoleculeAtomRange {
        using body_t    = std::remove_cv_t<TBody>;
        using vec_ref_t = std::conditional_t<std::is_const_v<TBody>, const std::vector<body_t>&, std::vector<body_t>&>;
        vec_ref_t bodies_;
    public:
        using iterator = MoleculeAtomIterator<TBody>;

        explicit MoleculeAtomRange(vec_ref_t bodies) : bodies_(bodies) {}

        iterator begin() const { return {bodies_.begin(), bodies_.end()}; }
        iterator end()   const { return {bodies_.end(),   bodies_.end()}; }
    };

    /**
     * @brief Forward iterator that traverses all water molecules across all bodies of a molecule in order,
     *        skipping bodies that carry no hydration layer.
     *        Use TBody=Body for a mutable iterator, TBody=const Body for a read-only iterator.
     */
    template<typename TBody>
    class MoleculeWaterIterator {
        using body_t      = std::remove_cv_t<TBody>;
        using ref_body_t  = std::conditional_t<std::is_const_v<TBody>, const body_t, body_t>;
        using body_iter_t = std::conditional_t<std::is_const_v<TBody>,
                                typename std::vector<body_t>::const_iterator,
                                typename std::vector<body_t>::iterator>;
        // get_waters() returns optional<reference_wrapper<vector<Water>>> or the const variant
        using opt_waters_t  = std::remove_reference_t<decltype(std::declval<ref_body_t&>().get_waters())>;
        // Extract the vector type via reference_wrapper::get()
        using water_vec_t   = std::remove_reference_t<decltype(std::declval<opt_waters_t>()->get())>;
        using water_iter_t  = std::conditional_t<std::is_const_v<water_vec_t>,
                                  typename water_vec_t::const_iterator,
                                  typename water_vec_t::iterator>;

    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = typename std::remove_cv_t<water_vec_t>::value_type;
        using pointer           = std::conditional_t<std::is_const_v<water_vec_t>, const value_type*, value_type*>;
        using reference         = std::conditional_t<std::is_const_v<water_vec_t>, const value_type&, value_type&>;

        MoleculeWaterIterator() = default;

        MoleculeWaterIterator(body_iter_t first, body_iter_t last) : body_cur_(first), body_end_(last) {
            advance_to_valid_body();
        }

        reference operator*()  const { return *water_cur_; }
        pointer   operator->() const { return &*water_cur_; }

        MoleculeWaterIterator& operator++() {
            ++water_cur_;
            skip_empty();
            return *this;
        }

        MoleculeWaterIterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }

        bool operator==(const MoleculeWaterIterator& rhs) const {
            if (body_cur_ == body_end_ && rhs.body_cur_ == rhs.body_end_) { return true;  }
            if (body_cur_ == body_end_ || rhs.body_cur_ == rhs.body_end_) { return false; }
            return body_cur_ == rhs.body_cur_ && water_cur_ == rhs.water_cur_;
        }

        bool operator!=(const MoleculeWaterIterator& rhs) const { return !(*this == rhs); }

    private:
        body_iter_t body_cur_{}, body_end_{};
        water_iter_t water_cur_{}, water_end_{};

        // Advance body_cur_ (starting from current position) until a body with non-empty waters
        // is found, or body_end_ is reached. Sets water_cur_/water_end_ on success.
        void advance_to_valid_body() {
            while (body_cur_ != body_end_) {
                auto w = body_cur_->get_waters();
                if (w.has_value() && !w->get().empty()) {
                    water_cur_ = w->get().begin();
                    water_end_ = w->get().end();
                    return;
                }
                ++body_cur_;
            }
        }

        void skip_empty() {
            while (body_cur_ != body_end_ && water_cur_ == water_end_) {
                ++body_cur_;
                advance_to_valid_body();
            }
        }
    };

    template<typename TBody>
    class MoleculeWaterRange {
        using body_t    = std::remove_cv_t<TBody>;
        using vec_ref_t = std::conditional_t<std::is_const_v<TBody>, const std::vector<body_t>&, std::vector<body_t>&>;
        vec_ref_t bodies_;
    public:
        using iterator = MoleculeWaterIterator<TBody>;

        explicit MoleculeWaterRange(vec_ref_t bodies) : bodies_(bodies) {}

        iterator begin() const { return {bodies_.begin(), bodies_.end()}; }
        iterator end()   const { return {bodies_.end(),   bodies_.end()}; }
    };
}