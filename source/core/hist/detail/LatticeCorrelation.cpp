// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

// we drive all parallelism through the global thread pool, so pocketfft must not spawn threads of its own.
// this has to be defined before the header is pulled in.
#define POCKETFFT_NO_MULTITHREADING

#include <hist/detail/LatticeCorrelation.h>
#include <math/Vector3.h>
#include <settings/GridSettings.h>
#include <utility/Console.h>
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
#include <new>
#include <string>

using namespace ausaxs;
using namespace ausaxs::hist;
using namespace ausaxs::hist::detail;
using ausaxs::hist::detail::lattice::Correlations;

namespace {
    using Coordinate = std::array<int32_t, 3>;

    // how far a point may sit from its lattice site, in lattice units, before we refuse to treat the set as a lattice
    constexpr double lattice_tolerance = 1e-6;

    /**
     * @brief The zero-padded box a correlation is evaluated in.
     */
    struct Box {
        std::array<std::size_t, 3> shape;  // the padded transform shape
        std::array<int32_t, 3> extent;     // the number of lattice sites the points span along each axis
        double spacing;                    // the lattice spacing in Ångström

        double cells() const {return static_cast<double>(shape[0])*static_cast<double>(shape[1])*static_cast<double>(shape[2]);}
    };

    /**
     * @brief One or two point sets expressed in integer coordinates on their common lattice.
     */
    struct Lattice {
        std::vector<Coordinate> first;
        std::vector<Coordinate> second;
        Box box;
    };

    // the smallest 5-smooth number >= n. pocketfft handles any length, but is considerably faster on these.
    std::size_t next_smooth(std::size_t n) {
        for (n = std::max<std::size_t>(n, 1);; ++n) {
            std::size_t remainder = n;
            for (std::size_t factor : {2, 3, 5}) {
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

        double origin[3] = {
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
        auto convert = [&res, &origin, inv_spacing] (const std::vector<Vector3<double>>& set, std::vector<Coordinate>& out) {
            out.resize(set.size());
            for (std::size_t n = 0; n < set.size(); ++n) {
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
            res.box.shape[k] = next_smooth(2*static_cast<std::size_t>(res.box.extent[k]) - 1);
        }
        return res;
    }

    /**
     * @brief Whether the transform buffers for the given box stay within settings::grid::exv::max_transform_memory.
     *
     * The buffers scale with the padded bounding box rather than the point count, so this is the one way the lattice
     * path can be worse than the pair loop it replaces. Elongated structures are where to expect it: the box grows
     * with the bounding volume while the pair count grows with the occupied volume.
     */
    bool fits_in_memory(const Box& box) {
        // one real box, plus a half-spectrum holding one complex value per two real cells
        const double bytes = box.cells()*(sizeof(double) + sizeof(std::complex<double>)/2.);
        const double budget = static_cast<double>(settings::grid::exv::max_transform_memory)*1024*1024;
        if (budget < bytes) {
            console::print_warning(
                "lattice::correlations: the excluded volume bounding box would need " + std::to_string(static_cast<int>(bytes/(1024*1024))) +
                "MB of transform buffers, which exceeds the " + std::to_string(settings::grid::exv::max_transform_memory) +
                "MB budget. Falling back to an explicit pair loop, which will be considerably slower."
            );
            return false;
        }
        return true;
    }

    /**
     * @brief The occupancy box and its half-spectrum, reused across the correlations of a single point set pair.
     */
    class Transform {
        public:
            Transform(const Box& box)
                : shape{box.shape[0], box.shape[1], box.shape[2]},
                  half_shape{box.shape[0], box.shape[1], box.shape[2]/2 + 1},
                  real_stride{
                      static_cast<std::ptrdiff_t>(shape[1]*shape[2]*sizeof(double)),
                      static_cast<std::ptrdiff_t>(shape[2]*sizeof(double)),
                      static_cast<std::ptrdiff_t>(sizeof(double))
                  },
                  spectrum_stride{
                      static_cast<std::ptrdiff_t>(half_shape[1]*half_shape[2]*sizeof(std::complex<double>)),
                      static_cast<std::ptrdiff_t>(half_shape[2]*sizeof(std::complex<double>)),
                      static_cast<std::ptrdiff_t>(sizeof(std::complex<double>))
                  },
                  real(shape[0]*shape[1]*shape[2]),
                  spectrum(half_shape[0]*half_shape[1]*half_shape[2])
            {}

            /**
             * @brief Replace the box contents with the autocorrelation of the indicator function of the given sets.
             *
             * Passing more than one set gives the autocorrelation of their combined occupancy, which by
             * A_combined = A_first + A_second + A_cross + A_cross^T is how the cross term is recovered without ever
             * holding two spectra at once.
             */
            void autocorrelate(std::initializer_list<const std::vector<Coordinate>*> sets) {
                std::fill(real.begin(), real.end(), 0.);
                for (const auto* set : sets) {
                    for (const auto& p : *set) {
                        real[(static_cast<std::size_t>(p[0])*shape[1] + p[1])*shape[2] + p[2]] += 1;
                    }
                }

                pocketfft::r2c(
                    shape, real_stride, spectrum_stride, {0, 1, 2}, pocketfft::FORWARD, real.data(), spectrum.data(), 1., 1
                );
                for (auto& z : spectrum) {z = std::norm(z);}

                // the multi-axis c2r allocates a full-size scratch buffer internally, so we spell out its two steps ourselves
                pocketfft::c2c(
                    half_shape, spectrum_stride, spectrum_stride, {0, 1}, pocketfft::BACKWARD, spectrum.data(), spectrum.data(), 1., 1
                );
                pocketfft::c2r(
                    shape, spectrum_stride, real_stride, 2, pocketfft::BACKWARD, spectrum.data(), real.data(),
                    1./static_cast<double>(real.size()), 1
                );
            }

            /**
             * @brief Radially bin the pair counts currently held in the box, adding them to @a out.
             *
             * Self-pairs are left out, since the callers add that term themselves. Every displacement has an exactly
             * known integer squared length, so the distances carry no accumulated error; the bin index is formed in
             * double precision, unlike the float32 of the pair loop, so the two can disagree on a bin edge.
             */
            void bin(const Box& box, double inv_bin_width, WeightedDistribution1D& out) const {
                for (int dx = -(box.extent[0]-1); dx < box.extent[0]; ++dx) {
                    const std::size_t ix = wrap(dx, shape[0]);
                    for (int dy = -(box.extent[1]-1); dy < box.extent[1]; ++dy) {
                        const std::size_t iy = wrap(dy, shape[1]);
                        const double* row = real.data() + (ix*shape[1] + iy)*shape[2];
                        const std::int64_t dxy2 = std::int64_t(dx)*dx + std::int64_t(dy)*dy;
                        for (int dz = -(box.extent[2]-1); dz < box.extent[2]; ++dz) {
                            const double pairs = row[wrap(dz, shape[2])];
                            if (pairs < 0.5) {continue;} // exactly zero up to the rounding error of the transform
                            const std::int64_t d2 = dxy2 + std::int64_t(dz)*dz;
                            if (d2 == 0) {continue;}
                            const double distance = std::sqrt(static_cast<double>(d2))*box.spacing;
                            // note that count is only 32 bits wide, exactly as in the pair loop this replaces
                            const unsigned int count = static_cast<unsigned int>(std::llround(pairs));
                            const auto index = static_cast<int32_t>(std::round(distance*inv_bin_width));
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
                for (double pairs : real) {worst = std::max(worst, std::abs(pairs - std::round(pairs)));}
                return worst;
            }

        private:
            pocketfft::shape_t shape;            // the padded real box
            pocketfft::shape_t half_shape;       // the spectrum left by an r2c along the last axis
            pocketfft::stride_t real_stride;
            pocketfft::stride_t spectrum_stride;
            std::vector<double> real;
            std::vector<std::complex<double>> spectrum;

            static std::size_t wrap(int d, std::size_t length) {
                return static_cast<std::size_t>((d + static_cast<int>(length)) % static_cast<int>(length));
            }
    };
}

std::optional<WeightedDistribution1D> hist::detail::lattice::self_correlation(
    const std::vector<Vector3<double>>& points, double spacing, double inv_bin_width, unsigned int bin_count)
{
    auto projected = project(points, {}, spacing);
    if (!projected.has_value()) {
        logging::log("lattice::self_correlation: the given points are not lattice-supported. Falling back to a pair loop.");
        return std::nullopt;
    }
    if (!fits_in_memory(projected->box)) {return std::nullopt;}

    WeightedDistribution1D out(bin_count);
    try {
        Transform transform(projected->box);
        transform.autocorrelate({&projected->first});
        assert(transform.rounding_error() < 0.1 && "lattice::self_correlation: the transform is losing integer precision");
        transform.bin(projected->box, inv_bin_width, out);
    } catch (const std::bad_alloc&) {
        console::print_warning(
            "lattice::self_correlation: could not allocate the transform buffers. Falling back to an explicit pair loop."
        );
        return std::nullopt;
    }
    logging::log("lattice::self_correlation: evaluated " + std::to_string(points.size()) + " excluded volume points by transform.");
    return out;
}

std::optional<Correlations> hist::detail::lattice::correlations(
    const std::vector<Vector3<double>>& first, const std::vector<Vector3<double>>& second,
    double spacing, double inv_bin_width, unsigned int bin_count)
{
    auto projected = project(first, second, spacing);
    if (!projected.has_value()) {
        logging::log("lattice::correlations: the given points are not lattice-supported. Falling back to a pair loop.");
        return std::nullopt;
    }
    if (!fits_in_memory(projected->box)) {return std::nullopt;}

    Correlations out{
        WeightedDistribution1D(bin_count),
        WeightedDistribution1D(bin_count),
        WeightedDistribution1D(bin_count)
    };
    try {
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
    } catch (const std::bad_alloc&) {
        console::print_warning(
            "lattice::correlations: could not allocate the transform buffers. Falling back to an explicit pair loop."
        );
        return std::nullopt;
    }
    logging::log(
        "lattice::correlations: evaluated " + std::to_string(first.size()) + " interior and " +
        std::to_string(second.size()) + " surface excluded volume points by transform."
    );
    return out;
}
