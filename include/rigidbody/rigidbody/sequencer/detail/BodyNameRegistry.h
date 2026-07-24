// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>

#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace ausaxs::rigidbody::sequencer::detail {
    /**
     * @brief Tracks the name to encoded index mapping used to resolve body (and symmetry replica) names in a sequencer script.
     */
    class BodyNameRegistry {
        public:
            /**
             * @brief Register a freshly loaded body's permanent default name ("bN"), together with an optional initial custom alias.
             */
            void add_body(int body, const std::string& alias = {});

            /**
             * @brief Register a symmetry replica's permanent canonical tag ("bXsYrZ", built from the body's own index) for a body's 
             * (symmetry, replica) slot. The first replica additionally gets the short "bXsY" alias.
             *
             * @param body The base body's index.
             * @param isymmetry The symmetry's index within the body (0-based).
             * @param replica The replica number (1-based).
             */
            void add_replica(int body, int isymmetry, int replica);

            /**
             * @brief Create or change an entity's custom alias, identified by any of its currently known names.
             *
             * If old_name is the entity's default name or its current alias, that alias slot is replaced: the default name itself is 
             * never removed. If old_name refers to an untracked (plain) name, it is simply renamed in place.
             */
            void rename(std::string_view old_name, std::string_view new_name);

            /**
             * @brief Drop the given base-body indices and reindex every surviving name (default, alias, and plain) to its shifted position.
             */
            void remove(std::vector<int> body_indices);

            bool contains(std::string_view name) const;
            unsigned int at(std::string_view name) const;

            /**
             * @brief Return the current addressable name for an encoded entity.
             *        This is its alias when renamed, otherwise its permanent default name.
             */
            std::string name(unsigned int index) const;

            /**
             * @brief Resolve a name to the (body, symmetry, replica) selector it refers to. 
             */
            BodySymmetrySelector resolve(std::string_view name) const;

            /**
             * @brief Resolve a name to a base body's index, rejecting names that refer to a symmetry replica.
             */
            int resolve_body(std::string_view name) const;

            std::unordered_map<std::string, unsigned int>::const_iterator begin() const;
            std::unordered_map<std::string, unsigned int>::const_iterator end() const;

            struct Group {
                std::string default_name;
                std::vector<std::string> others;
            };

            /**
             * @brief Group every known name by the exact (body, symmetry, replica) index it refers to.
             */
            std::vector<Group> group_by_index() const;

            /**
             * @brief Display names of the base bodies, ordered by body index. Each entry is the body's custom alias if it has one,
             * otherwise its default name ("bN"). Symmetry replicas are excluded.
             */
            std::vector<std::string> base_body_names() const;

        private:
            /**
             * @brief Register an arbitrary entity's permanent default name for an already-encoded index.
             */
            void add_entity(unsigned int index, const std::string& default_name, const std::string& alias = {});

            std::unordered_map<std::string, unsigned int> names;
            std::unordered_map<unsigned int, std::string> defaults;
            std::unordered_map<unsigned int, std::string> aliases;
    };
}
