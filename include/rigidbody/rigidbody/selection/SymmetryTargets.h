// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

#include <optional>
#include <unordered_map>
#include <vector>

namespace ausaxs::rigidbody::selection {
    /**
     * @brief The set of symmetry slots the optimizer is allowed to drive, indexed for selection.
     *
     * Not every symmetry a body declares can be optimized. A ReferenceSymmetryView forwards to the primary body that owns the shared symmetry and exposes no
     * parameters of its own, so perturbing it is a silent no-op; only the owning slot may be driven. Rather than have every select strategy rediscover that
     * rule, it is applied once here and the strategies draw from the resulting pool.
     *
     * The pool is derived from the molecule, so it goes stale whenever a body or a declared symmetry is added or removed. Rather than rebuild eagerly at every
     * such site - which would also pay for mutations nobody ever reads the result of - those sites call invalidate() and the pool is rebuilt lazily on the next
     * access. Adding a body or a symmetry anywhere without calling invalidate() leaves the optimizer driving a stale set of slots.
     */
    class SymmetryTargets {
        public:
            /**
             * @param molecule The molecule to derive the pool from. Must outlive this object.
             */
            SymmetryTargets(observer_ptr<const data::Molecule> molecule);

            /**
             * @brief One drivable symmetry slot. Named apart from BodySelectStrategy::Target, which is a whole move rather than just the slot it acts on.
             */
            struct Slot {
                unsigned int ibody;     //< the body declaring the symmetry
                unsigned int isymmetry; //< the symmetry's slot within that body
            };

            /**
             * @brief Mark the pool as stale, so the next access rebuilds it.
             *
             * Must be called by anything that adds or removes a body, or adds or removes a symmetry on one.
             */
            void invalidate() {dirty = true;}

            /**
             * @brief Every drivable slot, ordered by body and then by slot.
             */
            const std::vector<Slot>& all() const;

            /**
             * @brief The drivable slots of a single body, or an empty vector if it has none.
             */
            const std::vector<unsigned int>& body_targets(unsigned int ibody) const;

            /**
             * @brief Map a declared symmetry slot onto the drivable slot that actually backs it.
             *
             * Every body participating in a symmetry shared with others declares it, but only the owner holds the parameters; the rest hold non-owning views.
             * A user naming such a slot means the shared symmetry, so views resolve to their owner rather than being rejected. A slot that is drivable in its
             * own right maps to itself.
             *
             * @return The drivable slot, or nullopt if the slot does not exist or nothing can drive it.
             */
            std::optional<Slot> resolve(unsigned int ibody, unsigned int isymmetry) const;

            /**
             * @brief True if the molecule declares no drivable symmetry at all, in which case no symmetry-only step can accomplish anything.
             */
            bool empty() const {return all().empty();}

            std::size_t size() const {return all().size();}

        private:
            /**
             * @brief Rebuild the pool from the molecule's current symmetry set, if it has been invalidated since the last access.
             */
            void refresh() const;

            observer_ptr<const data::Molecule> molecule;

            // the pool is a cache of what the molecule already says, so reading it is a const operation even when it has to be rebuilt first
            mutable std::vector<Slot> targets;
            mutable std::unordered_map<unsigned int, std::vector<unsigned int>> per_body;
            mutable bool dirty = true;
    };
}
