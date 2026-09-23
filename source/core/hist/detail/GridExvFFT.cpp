// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/detail/GridExvFFT.h>

#include <math/Vector3.h>
#include <math/indexers/Indexer3D.h>
#include <settings/GeneralSettings.h>
#include <utility/Logging.h>
#include <utility/observer_ptr.h>

#include <pocketfft_hdronly.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <initializer_list>
#include <string>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;
using ausaxs::hist::detail::lattice::Correlations;

namespace {
    /**
     * @brief The zero-padded box a correlation is evaluated in.
     */
    struct Box {
        std::array<std::size_t, 3> shape;  // the padded transform shape
        std::array<int32_t, 3> extent;     // the number of lattice sites the points span along each axis
        double spacing;                    // the lattice spacing in Ångström
    };

    // the smallest 5-smooth number >= n. pocketfft is considerably faster on these.
    int next_smooth(int n) {
        for (n = std::max<int>(n, 1);; ++n) {
            int remainder = n;
            for (int factor : {2, 3, 5}) {
                while (remainder % factor == 0) {remainder /= factor;}
            }
            if (remainder == 1) {return n;}
        }
    }

    /**
     * @brief The box spanned by the given lattice sites, which are non-negative by construction.
     *        Padding each axis of the transform box to at least 2*extent-1 removes the circular wrap-around of the correlation.
     */
    Box make_box(std::initializer_list<observer_ptr<const std::vector<Vector3<int>>>> sets, double spacing) {
        Box box{.shape={}, .extent={0, 0, 0}, .spacing=spacing};
        for (const auto* set : sets) {
            for (const auto& p : *set) {
                for (int k = 0; k < 3; ++k) {
                    assert(0 <= p[k] && "lattice::make_box: lattice sites must be non-negative");
                    box.extent[k] = std::max(box.extent[k], p[k] + 1);
                }
            }
        }
        for (int k = 0; k < 3; ++k) {box.shape[k] = next_smooth(2*box.extent[k] - 1);}
        return box;
    }

    /**
     * @brief A 3D indexable view of the real box held in the transform buffer.
     *        Each row is L = padded_row_length reals long, of which only the leading real_dims[2] are live.
     */
    template<typename T>
    struct RealView : utility::indexer::Indexer3D<RealView<T>> {
        RealView(T* data, int N, int M, int L) : data(data), N(N), M(M), L(L) {}
        using utility::indexer::Indexer3D<RealView<T>>::index;

        T* data;
        int N, M, L;
    };

    /**
     * @brief The occupancy box, reused across the correlations of a single point set pair.
     *
     * Both transforms run in place inside a single buffer of one complex value per half-spectrum cell.
     */
    struct Transform {
        Transform(const Box& box) : 
            real_dims{box.shape[0], box.shape[1], box.shape[2]},
            spectrum_dims{box.shape[0], box.shape[1], box.shape[2]/2 + 1},
            padded_row_length(2*spectrum_dims[2]),
            cell_count(static_cast<double>(real_dims[0])*static_cast<double>(real_dims[1])*static_cast<double>(real_dims[2])),
            real_layout{ // stride in x, y, z
                static_cast<std::ptrdiff_t>(real_dims[1]*padded_row_length*sizeof(double)),
                static_cast<std::ptrdiff_t>(padded_row_length*sizeof(double)),
                static_cast<std::ptrdiff_t>(sizeof(double))
              },
            spectrum_layout{ // stride in x, y, z
                static_cast<std::ptrdiff_t>(spectrum_dims[1]*spectrum_dims[2]*sizeof(std::complex<double>)),
                static_cast<std::ptrdiff_t>(spectrum_dims[2]*sizeof(std::complex<double>)),
                static_cast<std::ptrdiff_t>(sizeof(std::complex<double>))
            },
            buffer(spectrum_dims[0]*spectrum_dims[1]*spectrum_dims[2])
        {}

        /**
         * @brief Replace the box contents with the autocorrelation of the indicator function of the given sets.
         *
         * Let f be the occupancy of the box: f(r) = 1 if a point sits on site r, 0 otherwise. Its autocorrelation
         *
         *     A(d) = sum_r f(r) f(r + d)
         *
         * is exactly the number of ordered point pairs separated by the displacement d. Evaluated directly that is the O(N^2) pair loop. 
         * The Wiener-Khinchin theorem instead gives A = F^-1[|F[f]|^2], i.e. one forward transform, a pointwise squared magnitude, and 
         * one inverse transform - O(M log M) in the box volume M.
         *
         * On return, the real box holds A(d) at index d, where a negative displacement component -d is stored at index n-d along its 
         * axis (the transform is periodic). The zero-padding of the box to at least 2*extent-1 guarantees that positive and negative 
         * displacements never land on the same index. See bin() for the readout.
         */
        void autocorrelate(std::initializer_list<observer_ptr<const std::vector<Vector3<int>>>> sets) {
            // build the occupancy f. the whole buffer is cleared, including the padding at the end of each row, since a previous call may 
            // have left spectrum data there. passing several sets gives the occupancy of their union; the sets never share a site, so every 
            // cell stays either 0 or 1.
            std::ranges::fill(buffer, std::complex<double>(0, 0));
            auto real = real_box();
            for (const auto* set : sets) {
                for (const auto& p : *set) {
                    assert(real.index(p[0], p[1], p[2]) == 0 && "lattice::Transform: the sets must not share a site");
                    real.index(p[0], p[1], p[2]) += 1;
                }
            }

            // forward transform, F = F[f]. since f is real, its spectrum is Hermitian, F(-k) = conj(F(k)), so only the non-negative half 
            // along the last axis (n/2+1 values) is computed and stored. this runs in place: the spectrum overwrites the real box in the 
            // same buffer.
            pocketfft::r2c(real_dims, real_layout, spectrum_layout, {0, 1, 2}, pocketfft::FORWARD, data(), buffer.data(), 1., settings::general::threads);

            // the power spectrum |F|^2 is the transform of the autocorrelation. it is real, which is why A comes out symmetric, 
            // A(d) = A(-d): every pair is counted once in each direction.
            for (auto& z : buffer) {z = std::norm(z);}

            // inverse transform, A = F^-1[|F|^2]. the one-call multi-axis c2r allocates a full-size scratch buffer internally, so we spell 
            // out its two stages ourselves:
            //   1. a complex-to-complex inverse along the first two axes, which is genuinely in place.
            //   2. a complex-to-real inverse along the last axis, which expands each half-row of n/2+1 complex values back into n reals. 
            //      it is line-buffered, so it too is safe in place.
            // pocketfft's transforms are unnormalised, so a forward-inverse round trip scales by the number of cells;
            // the final stage divides that back out, which is what turns the result into integer pair counts.
            pocketfft::c2c(spectrum_dims, spectrum_layout, spectrum_layout, {0, 1}, pocketfft::BACKWARD, buffer.data(), buffer.data(), 1., settings::general::threads);
            pocketfft::c2r(real_dims, spectrum_layout, real_layout, 2, pocketfft::BACKWARD, buffer.data(), data(), 1./cell_count, settings::general::threads);
        }

        /**
         * @brief Radially bin the pair counts currently held in the box, adding them to @a out.
         *
         * Every cell of the box holds the number of pairs with one specific integer displacement (dx, dy, dz), so the
         * histogram is filled by visiting each displacement once, computing its length, and adding its count to the
         * matching distance bin. Self-pairs (the zero displacement) are left out, since the callers add them.
         */
        void bin(const Box& box, double inv_bin_width, WeightedDistribution1D& out) const {
            double spacing = box.spacing;
            auto real = real_box();

            // two points spanning extent sites along an axis can differ by at most extent-1 there, so only the
            // displacements in [-(extent-1), extent-1] can hold pairs. everything else is zero-padding.
            // the negative displacements were stored wrapped around to the far end of each axis; wrap() undoes that.
            for (int dx = -(box.extent[0]-1); dx < box.extent[0]; ++dx) {
                int ix = wrap(dx, static_cast<int>(real_dims[0]));
                for (int dy = -(box.extent[1]-1); dy < box.extent[1]; ++dy) {
                    int iy = wrap(dy, static_cast<int>(real_dims[1]));
                    // the extents are bounded far below int overflow by the memory needed for the box itself
                    int dxy2 = dx*dx + dy*dy;
                    for (int dz = -(box.extent[2]-1); dz < box.extent[2]; ++dz) {
                        double pairs = real.index(ix, iy, wrap(dz, static_cast<int>(real_dims[2])));
                        // the counts are integers up to the round-off of the transform, so this is an exact zero test
                        if (pairs < 0.5) {continue;}
                        int d2 = dxy2 + dz*dz;
                        if (d2 == 0) {continue;}
                        double distance = std::sqrt(static_cast<double>(d2))*spacing;
                        // note that count is only 32 bits wide, exactly as in the pair loop this replaces
                        auto count = static_cast<unsigned int>(std::llround(pairs));
                        auto index = static_cast<int32_t>(std::round(distance*inv_bin_width));
                        assert(index < static_cast<int32_t>(out.size()) && "lattice: distance bin out of range");
                        out.add_index(index, WeightedEntry(count, count, count*distance));
                    }
                }
            }
        }

        /**
         * @brief The largest deviation of a pair count from an integer, as a check on the numerical margin.
         */
        double rounding_error() const {
            double worst = 0;
            auto real = real_box();
            // only the leading real_dims[2] reals of each row are live; the rest of the padded row is stale spectrum
            for (int i = 0; i < static_cast<int>(real_dims[0]); ++i) {
                for (int j = 0; j < static_cast<int>(real_dims[1]); ++j) {
                    for (int k = 0; k < static_cast<int>(real_dims[2]); ++k) {
                        worst = std::max(worst, std::abs(real.index(i, j, k) - std::round(real.index(i, j, k))));
                    }
                }
            }
            return worst;
        }

        pocketfft::shape_t real_dims;           // the padded real box
        pocketfft::shape_t spectrum_dims;       // the spectrum left by an r2c along the last axis
        std::size_t padded_row_length;          // reals per row of the aliased real box, 2*spectrum_dims[2]
        double cell_count;                      // the number of real cells, i.e. the inverse transform normalisation
        pocketfft::stride_t real_layout;        // bytes to step one cell along each axis of the real box
        pocketfft::stride_t spectrum_layout;    // bytes to step one cell along each axis of the spectrum
        std::vector<std::complex<double>> buffer;

        // the real box shares the spectrum's allocation; see the struct comment for why that is safe
        double* data() {return reinterpret_cast<double*>(buffer.data());}
        const double* data() const {return reinterpret_cast<const double*>(buffer.data());}
        RealView<double> real_box() {return {data(), static_cast<int>(real_dims[0]), static_cast<int>(real_dims[1]), static_cast<int>(padded_row_length)};}
        RealView<const double> real_box() const {return {data(), static_cast<int>(real_dims[0]), static_cast<int>(real_dims[1]), static_cast<int>(padded_row_length)};}
        static int wrap(int d, int length) {return (d + length) % length;}
    };
}

WeightedDistribution1D hist::detail::lattice::self_correlation(const grid::exv::GridExcludedVolume& exv, double inv_bin_width, int bin_count) {
    WeightedDistribution1D out(bin_count);
    Box box = make_box({&exv.interior_sites}, exv.spacing);
    Transform transform(box);
    transform.autocorrelate({&exv.interior_sites});
    assert(transform.rounding_error() < 0.1 && "lattice::self_correlation: the transform is losing integer precision");
    transform.bin(box, inv_bin_width, out);
    logging::log("lattice::self_correlation: evaluated " + std::to_string(exv.interior_sites.size()) + " excluded volume points by transform.");
    return out;
}

Correlations hist::detail::lattice::correlations(const grid::exv::GridExcludedVolume& exv, double inv_bin_width, int bin_count) {
    Correlations out{
        .first=WeightedDistribution1D(bin_count),
        .second=WeightedDistribution1D(bin_count),
        .cross=WeightedDistribution1D(bin_count)
    };
    Box box = make_box({&exv.interior_sites, &exv.surface_sites}, exv.spacing);
    Transform transform(box);

    transform.autocorrelate({&exv.interior_sites});
    assert(transform.rounding_error() < 0.1 && "lattice::correlations: the transform is losing integer precision");
    transform.bin(box, inv_bin_width, out.first);

    transform.autocorrelate({&exv.surface_sites});
    transform.bin(box, inv_bin_width, out.second);

    // the cross term is whatever the correlation of the combined set holds beyond the two self-correlations.
    // both sides are exact integer counts over the same displacements, so the subtraction is exact.
    transform.autocorrelate({&exv.interior_sites, &exv.surface_sites});
    transform.bin(box, inv_bin_width, out.cross);
    out.cross -= out.first;
    out.cross -= out.second;
    logging::log(
        "lattice::correlations: evaluated " + std::to_string(exv.interior_sites.size()) + " interior and " +
        std::to_string(exv.surface_sites.size()) + " surface excluded volume points by transform."
    );
    return out;
}
