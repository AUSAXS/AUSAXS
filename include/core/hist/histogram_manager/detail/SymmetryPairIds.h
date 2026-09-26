// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <container/Container2D.h>

#include <cassert>
#include <vector>

namespace ausaxs::hist::detail {
	/**
	 * @brief The store ids of the atom-atom results, one for each pair of symmetries of each body pair.
	 *        The symmetry index isym runs over the main body (0) and its symmetries (1..).
	 */
	class SymmetryPairIds {
		public:
			SymmetryPairIds() = default;

			/**
			 * @brief Lay out the pairs of bodies with @a sym_counts symmetry indices each, and give each calculated pair the id returned by @a allocate.
			 */
			template<typename Allocate>
			SymmetryPairIds(const std::vector<int>& sym_counts, Allocate&& allocate) {
				int n = static_cast<int>(sym_counts.size());
				blocks = container::Container2D<container::Container2D<int>>(n, n);
				for (int ibody1 = 0; ibody1 < n; ++ibody1) {
					for (int ibody2 = 0; ibody2 <= ibody1; ++ibody2) {
						auto& block = blocks(ibody1, ibody2);
						block = container::Container2D<int>(sym_counts[ibody1], sym_counts[ibody2]);
						for (int isym1 = 0; isym1 < sym_counts[ibody1]; ++isym1) {
							for (int isym2 = 0; isym2 < sym_counts[ibody2]; ++isym2) {
								block(isym1, isym2) = calculated(ibody1, isym1, ibody2, isym2) ? allocate() : -1;
							}
						}
					}
				}
			}

			int id(int ibody1, int isym1, int ibody2, int isym2) const {
				assert(ibody2 <= ibody1 && "SymmetryPairIds::id: expected a body pair in the lower triangle");
				const auto& block = blocks(ibody1, ibody2);
				assert(0 <= isym1 && isym1 < block.size_x() && 0 <= isym2 && isym2 < block.size_y() && "SymmetryPairIds::id: symmetry index out of range; symmetries may not be added after the first calculation");
				int id = block(isym1, isym2);
				assert(id != -1 && "SymmetryPairIds::id: expected a symmetry pair in the lower triangle");
				return id;
			}

			/**
			 * @brief Call @a f with the id of every calculated pair.
			 */
			template<typename F>
			void for_each_id(F&& f) const {
				for (const auto& block : blocks) {
					for (int id : block) {if (id != -1) {f(id);}}
				}
			}

		private:
			static bool calculated(int ibody1, int isym1, int ibody2, int isym2) {
				return ibody1 != ibody2 || isym2 < isym1 || (isym1 == 0 && isym2 == 0);
			}

			container::Container2D<container::Container2D<int>> blocks; // [ibody1][ibody2][isym1][isym2]; -1 for a pair that is not calculated, and no block above the diagonal
	};
}
