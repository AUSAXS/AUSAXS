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
        using atom_iter_t = typename atom_vec_t::iterator;

    public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = typename atom_vec_t::value_type;
        using pointer           = typename atom_vec_t::pointer;
        using reference         = typename atom_vec_t::reference;

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
}