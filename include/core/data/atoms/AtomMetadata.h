// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility/Exceptions.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace ausaxs::data {
    /**
     * @brief Backbone classification of an atom.
     */
    enum class backbone_t : std::uint8_t {
        none,       //< not a backbone atom
        n,          //< backbone amide nitrogen
        c_alpha,    //< C-alpha backbone carbon
        c,          //< backbone carbonyl carbon
        o,          //< backbone carbonyl oxygen
    };

    /**
     * @brief Optional per-atom metadata for a Body.
     *
     *        Whenever the atom vector of a Body is reshaped, its metadata must be reshaped identically or the parallel indexing is lost. The operations doing so 
     *        are provided as members here rather than at the call sites, so that adding a field only requires extending _visit below.
     */
    struct AtomMetadata {
        std::optional<std::vector<backbone_t>> backbone;     //< engaged iff store_backbone
        std::optional<std::vector<int>>        residue_seq;  //< residue sequence id; engaged iff store_residue_seq
        std::optional<std::vector<char>>       chain_id;     //< source chain identifier; engaged iff store_chain_id
        std::optional<std::vector<float>>      occupancy;    //< engaged iff store_occupancy

        // which fields are retained when a structure is loaded.
        inline static bool store_backbone    = true;
        inline static bool store_residue_seq = true;
        inline static bool store_chain_id    = true;
        inline static bool store_occupancy   = false;

        // total number of fields, engaged or not. Must be kept in sync with _visit.
        static constexpr std::size_t field_count = 4;

        /**
         * @brief Get the number of currently engaged fields.
         */
        [[nodiscard]] std::size_t active_count() const {
            std::size_t count = 0;
            _visit(*this, [&count] (const auto& field) {count += field.has_value();});
            return count;
        }

        /**
         * @brief Check if no field is engaged. Empty metadata carries no information and does not need to be attached to a Body at all.
         */
        [[nodiscard]] bool empty() const {return active_count() == 0;}

        /**
         * @brief Get the number of atoms this metadata is indexed by, or 0 if no field is engaged.
         */
        [[nodiscard]] std::size_t size() const {
            std::size_t size = 0;
            _visit(*this, [&size] (const auto& field) {
                if (!field) {return;}
                assert((size == 0 || size == field->size()) && "AtomMetadata::size: fields are not parallel-indexed.");
                size = field->size();
            });
            return size;
        }

        /**
         * @brief Shorten every engaged field to its first n elements.
         */
        void resize(std::size_t n) {
            assert(n <= size() && "AtomMetadata::resize: cannot grow metadata, since there is no data to grow it with.");
            _visit(*this, [n] (auto& field) {if (field) {field->resize(n);}});
        }

        /**
         * @brief Get the [begin, end) subrange of every engaged field.
         */
        [[nodiscard]] AtomMetadata subrange(std::size_t begin, std::size_t end) const {
            assert(begin <= end && end <= size() && "AtomMetadata::subrange: range is not contained in the metadata.");
            AtomMetadata result;
            _visit(result, *this, [begin, end] (auto& dst, const auto& src) {
                if (src) {dst.emplace(src->begin()+begin, src->begin()+end);}
            });
            return result;
        }

        /**
         * @brief Append the contents of another metadata instance to this one.
         *
         *        A field consistently absent on both sides is fine and stays absent. A field engaged on only one side means the two were tracked under different
         *        settings, which would desync the result from the atom vector, so that case is an error instead of a silent drop.
         */
        void append(const AtomMetadata& other) {
            _visit(*this, other, [] (auto& dst, const auto& src) {
                if (dst.has_value() != src.has_value()) {
                    throw except::invalid_argument("AtomMetadata::append: cannot append metadata with a different set of engaged fields.");
                }
                if (dst) {dst->insert(dst->end(), src->begin(), src->end());}
            });
        }

        private:
            // Invoke f on every field of self. f must accept any of the field types, i.e. be a generic lambda.
            template<typename S, typename F>
            static void _visit(S& self, F&& f) {
                f(self.backbone);
                f(self.residue_seq);
                f(self.chain_id);
                f(self.occupancy);
            }

            // Invoke f on every field of a paired up with the corresponding field of b. f must accept any of the field types, i.e. be a generic lambda.
            template<typename S1, typename S2, typename F>
            static void _visit(S1& a, S2& b, F&& f) {
                f(a.backbone,    b.backbone);
                f(a.residue_seq, b.residue_seq);
                f(a.chain_id,    b.chain_id);
                f(a.occupancy,   b.occupancy);
            }
    };
}
