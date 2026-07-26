// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <utility/observer_ptr.h>

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
             * @brief One drivable symmetry slot.
             */
            struct Target {
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
            const std::vector<Target>& all();

            /**
             * @brief The drivable slots of a single body, or an empty vector if it has none.
             */
            const std::vector<unsigned int>& body_targets(unsigned int ibody);

            /**
             * @brief True if the molecule declares no drivable symmetry at all, in which case no symmetry-only step can accomplish anything.
             */
            bool empty() {return all().empty();}

            std::size_t size() {return all().size();}

        private:
            /**
             * @brief Rebuild the pool from the molecule's current symmetry set, if it has been invalidated since the last access.
             */
            void refresh();

            observer_ptr<const data::Molecule> molecule;
            std::vector<Target> targets;
            std::unordered_map<unsigned int, std::vector<unsigned int>> per_body;
            bool dirty = true;
    };
}
