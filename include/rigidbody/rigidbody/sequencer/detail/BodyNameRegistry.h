// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>

#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace ausaxs::rigidbody::sequencer::detail {
    /**
     * @brief Tracks the name to encoded index mapping used to resolve body (and symmetry replica) names in a sequencer script.
     *
     * Every registered entity - a base body or one of its symmetry replicas - owns exactly one Entry, keyed by its encoded index. Names are the Entry's
     * contents rather than a store of their own; the name lookup is a pure inverse of the entries, rebuilt from them whenever indices shift, so the two
     * cannot drift apart.
     */
    class BodyNameRegistry {
        public:
            /**
             * @brief The names identifying one entity.
             */
            struct Entry {
                std::string default_name; //< permanent; minted when the entity is registered and never changed afterwards
                std::string alias;        //< user-chosen; empty until the entity is renamed

                /**
                 * @brief The name to address and display the entity by: its alias once renamed, otherwise its permanent default name.
                 */
                const std::string& display_name() const {return alias.empty() ? default_name : alias;}
            };

            /**
             * @brief Register a body's permanent default name ("bN"), together with an optional initial custom alias.
             *
             * N is drawn from a monotonic counter, not from the body index: indices shift when bodies are erased and are reused by bodies created later
             * (split, copy), whereas default names are permanent. A plain load of N bodies therefore still yields exactly "b1".."bN", while bodies created
             * later continue the sequence rather than colliding with a survivor that kept its original name. The counter is never rewound, so a name retired
             * by remove() is never reissued to a different body.
             *
             * @throws std::runtime_error If the alias is already in use by another entity.
             */
            void add_body(int body, const std::string& alias = {});

            /**
             * @brief Register a body under an entry inherited from a body it succeeds, rather than minting a fresh default name.
             *
             * Used when one body is replaced by others that continue its identity, as when `split` hands the original body's entry to its leading fragment.
             * The inherited names must have been released first (by remove()), or they would collide with the body still holding them.
             *
             * @throws std::runtime_error If either inherited name is still in use.
             */
            void add_body(int body, Entry inherited);

            /**
             * @brief Register a symmetry replica's permanent canonical tag ("<base>sYrZ", built from the base body's own default name) for a body's
             * (symmetry, replica) slot. The first replica additionally gets the short "<base>sY" alias.
             *
             * @param body The base body's index.
             * @param isymmetry The symmetry's index within the body (0-based).
             * @param replica The replica number (1-based).
             *
             * @throws std::runtime_error If the base body is not registered, or one of the generated tags is already in use.
             */
            void add_replica(int body, int isymmetry, int replica);

            /**
             * @brief Set an entity's custom alias, identified by any of its currently known names.
             *
             * Only the alias slot changes; the permanent default name is never removed. Renaming an entity back to its own default name clears the alias.
             *
             * @throws std::runtime_error If old_name is unknown, or new_name is already in use by a different entity.
             */
            void rename(std::string_view old_name, std::string_view new_name);

            /**
             * @brief Drop the given base-body indices, along with their replicas, and reindex every surviving entity to its shifted position.
             */
            void remove(std::vector<int> body_indices);

            bool contains(std::string_view name) const;

            /**
             * @brief Resolve a name to the (body, symmetry, replica) selector it refers to.
             *
             * @throws std::runtime_error If the name is unknown.
             */
            BodySymmetrySelector resolve(std::string_view name) const;

            /**
             * @brief Resolve a name to a base body's index, rejecting names that refer to a symmetry replica.
             */
            int resolve_body(std::string_view name) const;

            /**
             * @brief The names of the entity at an encoded index (see to_index).
             *
             * @throws std::out_of_range If nothing is registered at that index.
             */
            const Entry& entry(int index) const;

            /**
             * @brief Every registered entity, keyed and ordered by encoded index, so each body is immediately followed by its own replicas.
             */
            const std::map<int, Entry>& all() const;

            /**
             * @brief Display names of the base bodies, ordered by body index. Symmetry replicas are excluded.
             */
            std::vector<std::string> base_body_names() const;

        private:
            /**
             * @brief Register an entity's names at an already-encoded index.
             *
             * @throws std::runtime_error If either name is already in use. Registering a name twice would silently hide one of the two entities behind the
             *         other, leaving it unaddressable, so it is rejected rather than ignored.
             */
            void add_entity(int index, Entry entry);

            /**
             * @brief Whether anything at all - the body itself or one of its replicas - is registered for a base body index.
             */
            bool has_body(int body) const;

            /**
             * @brief Rebuild the name lookup from the entries, after their indices have shifted.
             */
            void rebuild_lookup();

            std::map<int, Entry> entries;                //< the naming state, keyed and ordered by encoded index
            std::unordered_map<std::string, int> lookup; //< pure inverse of `entries`: every name they hold maps back to its index
            unsigned int bodies_registered = 0;          //< monotonic; source of the "bN" default names. Never rewound, see add_body.
    };
}
