// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <hist/detail/GridExvFFT.h>
#include <math/Vector3.h>
#include <utility/Logging.h>

#include <pocketfft_hdronly.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <string>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;
using ausaxs::hist::detail::lattice::Correlations;

namespace {
    // how far a point may sit from its lattice site, in lattice units, before we refuse to treat the set as a lattice
    constexpr double lattice_tolerance = 1e-6;

    /**
     * @brief The zero-padded box a correlation is evaluated in.
     */
    struct Box {
        std::array<std::size_t, 3> shape;  // the padded transform shape
        std::array<int32_t, 3> extent;     // the number of lattice sites the points span along each axis
        double spacing;                    // the lattice spacing in Ångström
    };

    /**
     * @brief One or two point sets expressed in integer coordinates on their common lattice.
     */
    struct Lattice {
        std::vector<Vector3<int>> first;
        std::vector<Vector3<int>> second;
        Box box;
    };

    // the smallest 5-smooth number >= n. pocketfft handles any length, but is considerably faster on these.
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
     * @brief Express the given point sets in integer coordinates on their common cubic lattice.
     *
     * The origin is the lower corner of the combined bounding box, so both sets end up on a single lattice with
     * non-negative coordinates. Padding each axis of the transform box to at least 2*extent-1 removes the circular
     * wrap-around of the correlation, which is what makes the result exact rather than approximate.
     *
     * @return std::nullopt if any point does not sit on the lattice.
     */
    std::optional<Lattice> project(
        const std::vector<Vector3<double>>& first, const std::vector<Vector3<double>>& second, double spacing)
    {
        if (spacing <= 0 || (first.empty() && second.empty())) {return std::nullopt;}

        std::array<double, 3> origin = {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max()
        };
        for (const auto* set : {&first, &second}) {
            for (const auto& p : *set) {
                for (int k = 0; k < 3; ++k) {origin[k] = std::min(origin[k], p[k]);}
            }
        }

        Lattice res;
        res.box.spacing = spacing;
        res.box.extent = {0, 0, 0};
        const double inv_spacing = 1/spacing;
        auto convert = [&res, &origin, inv_spacing] (const std::vector<Vector3<double>>& set, std::vector<Vector3<int>>& out) {
            out.resize(set.size());
            for (int n = 0; n < static_cast<int>(set.size()); ++n) {
                for (int k = 0; k < 3; ++k) {
                    double site = (set[n][k] - origin[k])*inv_spacing;
                    double rounded = std::round(site);
                    if (lattice_tolerance < std::abs(site - rounded)) {return false;}
                    out[n][k] = static_cast<int32_t>(rounded);
                    res.box.extent[k] = std::max(res.box.extent[k], out[n][k] + 1);
                }
            }
            return true;
        };
        if (!convert(first, res.first) || !convert(second, res.second)) {return std::nullopt;}

        for (int k = 0; k < 3; ++k) {
            res.box.shape[k] = next_smooth(2*res.box.extent[k] - 1);
        }
        return res;
    }

    /**
     * @brief The occupancy box, reused across the correlations of a single point set pair.
     *
     * Both transforms run in place inside a single buffer of one complex value per half-spectrum cell, i.e. a little
     * over 8 bytes per padded real cell. The two obvious spellings each cost a further copy of the box and are
     * deliberately avoided: transforming out-of-place into a separate spectrum array doubles the peak, and pocketfft's
     * multi-axis c2r is a wrapper that allocates a full half-spectrum temporary on every call, which doubles it again.
     *
     * The real box is therefore aliased into the same allocation, laid out FFTW-style: a padded row length of 2*(n/2+1) reals
     * of which the leading n are live. That is what makes the forward r2c safe in place - general_r2c copies each line
     * into scratch before writing the corresponding output line, and with this layout a line's input and output occupy
     * exactly the same row, so no line can clobber another's input on either the scalar or the vectorised path. The
     * inverse is spelled out as the two stages the wrapper would have run: the leading-axis c2c is genuinely in place,
     * and the trailing single-axis c2r is line-buffered just like r2c.
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
         * Passing more than one set gives the autocorrelation of their combined occupancy, which by
         * A_combined = A_first + A_second + A_cross + A_cross^T is how the cross term is recovered without ever
         * holding two spectra at once.
         */
        void autocorrelate(std::initializer_list<const std::vector<Vector3<int>>*> sets) {
            std::ranges::fill(buffer, std::complex<double>(0, 0));
            double* real = data();
            for (const auto* set : sets) {
                for (const auto& p : *set) {
                    real[(static_cast<std::size_t>(p[0])*real_dims[1] + p[1])*padded_row_length + p[2]] += 1;
                }
            }

            pocketfft::r2c(real_dims, real_layout, spectrum_layout, {0, 1, 2}, pocketfft::FORWARD, real, buffer.data(), 1., 1);
            for (auto& z : buffer) {z = std::norm(z);}

            // the multi-axis c2r allocates a full-size scratch buffer internally, so we spell out its two steps ourselves
            pocketfft::c2c(spectrum_dims, spectrum_layout, spectrum_layout, {0, 1}, pocketfft::BACKWARD, buffer.data(), buffer.data(), 1., 1);
            pocketfft::c2r(real_dims, spectrum_layout, real_layout, 2, pocketfft::BACKWARD, buffer.data(), real, 1./cell_count, 1);
        }

        /**
         * @brief Radially bin the pair counts currently held in the box, adding them to @a out.
         */
        void bin(const Box& box, double inv_bin_width, WeightedDistribution1D& out) const {
            double spacing = box.spacing;
            for (int dx = -(box.extent[0]-1); dx < box.extent[0]; ++dx) {
                int ix = wrap(dx, static_cast<int>(real_dims[0]));
                for (int dy = -(box.extent[1]-1); dy < box.extent[1]; ++dy) {
                    int iy = wrap(dy, static_cast<int>(real_dims[1]));
                    const auto* row = data() + (ix*real_dims[1] + iy)*padded_row_length;
                    // the extents are bounded far below int overflow by the memory needed for the box itself
                    int dxy2 = dx*dx + dy*dy;
                    for (int dz = -(box.extent[2]-1); dz < box.extent[2]; ++dz) {
                        double pairs = row[wrap(dz, static_cast<int>(real_dims[2]))];
                        if (pairs < 0.5) {continue;} // exactly zero up to the rounding error of the transform
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
            // only the leading real_dims[2] reals of each row are live; the rest of the padded row is stale spectrum
            for (std::size_t i = 0; i < real_dims[0]; ++i) {
                for (std::size_t j = 0; j < real_dims[1]; ++j) {
                    const double* row = data() + (i*real_dims[1] + j)*padded_row_length;
                    for (std::size_t k = 0; k < real_dims[2]; ++k) {
                        worst = std::max(worst, std::abs(row[k] - std::round(row[k])));
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

        static int wrap(int d, int length) {
            return (d + length) % length;
        }
    };
}

std::optional<WeightedDistribution1D> hist::detail::lattice::self_correlation(
    const std::vector<Vector3<double>>& points, double spacing, double inv_bin_width, int bin_count)
{
    auto projected = project(points, {}, spacing);
    if (!projected.has_value()) {
        logging::log("lattice::self_correlation: the given points are not lattice-supported. Falling back to a pair loop.");
        return std::nullopt;
    }
    WeightedDistribution1D out(bin_count);
    Transform transform(projected->box);
    transform.autocorrelate({&projected->first});
    assert(transform.rounding_error() < 0.1 && "lattice::self_correlation: the transform is losing integer precision");
    transform.bin(projected->box, inv_bin_width, out);
    logging::log("lattice::self_correlation: evaluated " + std::to_string(points.size()) + " excluded volume points by transform.");
    return out;
}

std::optional<Correlations> hist::detail::lattice::correlations(
    const std::vector<Vector3<double>>& first, const std::vector<Vector3<double>>& second,
    double spacing, double inv_bin_width, int bin_count)
{
    auto projected = project(first, second, spacing);
    if (!projected.has_value()) {
        logging::log("lattice::correlations: the given points are not lattice-supported. Falling back to a pair loop.");
        return std::nullopt;
    }
    Correlations out{
        .first=WeightedDistribution1D(bin_count),
        .second=WeightedDistribution1D(bin_count),
        .cross=WeightedDistribution1D(bin_count)
    };
    Transform transform(projected->box);

    transform.autocorrelate({&projected->first});
    assert(transform.rounding_error() < 0.1 && "lattice::correlations: the transform is losing integer precision");
    transform.bin(projected->box, inv_bin_width, out.first);

    transform.autocorrelate({&projected->second});
    transform.bin(projected->box, inv_bin_width, out.second);

    // the cross term is whatever the correlation of the combined set holds beyond the two self-correlations.
    // both sides are exact integer counts over the same displacements, so the subtraction is exact.
    transform.autocorrelate({&projected->first, &projected->second});
    transform.bin(projected->box, inv_bin_width, out.cross);
    out.cross -= out.first;
    out.cross -= out.second;
    logging::log(
        "lattice::correlations: evaluated " + std::to_string(first.size()) + " interior and " +
        std::to_string(second.size()) + " surface excluded volume points by transform."
    );
    return out;
}
