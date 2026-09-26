// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <utility>

namespace ausaxs::utility::indexer {
    /**
     * @brief The index space of a pair (i, j) of equally sized leading dimensions.
     *        Square: every ordered pair has its own element. Triangular: only unordered pairs are stored, so (i, j) and (j, i)
     *        name the same element. A loop over a triangular container must visit each pair once (j >= i), or it double counts.
     */
    enum class Shape {Square, Triangular};

    namespace triangular {
        // the number of unordered pairs (i, j) of an N x N index space, including the diagonal
        constexpr int pair_count(int N) {return N*(N+1)/2;}

        // the position of the unordered pair (i, j), laid out row by row with row i holding j = i, ..., N-1
        constexpr int pair_index(int i, int j, int N) {
            if (j < i) {std::swap(i, j);}
            return i*N - i*(i-1)/2 + (j - i);
        }
    }

    namespace square {
        // the number of ordered pairs (i, j) of an N x M index space
        constexpr int pair_count(int N, int M) {return N*M;}

        // the position of the ordered pair (i, j), laid out row by row with row i holding j = 0, ..., M-1
        constexpr int pair_index(int i, int j, int M) {return i*M + j;}
    }
}
