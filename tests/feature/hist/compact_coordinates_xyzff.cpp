#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/detail/CompactCoordinatesFF.h>
#include <hist/detail/data/CompactCoordinatesXYZFF.h>
#include <constants/Constants.h>
#include <math/Vector3.h>
#include <form_factor/FormFactorType.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

using namespace ausaxs;
using namespace hist::detail;
using namespace hist::detail::xyzff;

TEST_CASE("CompactCoordinatesFF<vbw>: component storage") {
    SECTION("positions and form factor types round-trip") {
        std::vector<data::AtomFF> atoms = {
            data::AtomFF({1, 2, 3}, form_factor::form_factor_t::C),
            data::AtomFF({4, 5, 6}, form_factor::form_factor_t::O)
        };
        CompactCoordinatesFF<false> data(atoms);
        REQUIRE(data.size() == 2);
        CHECK(data.x(0) == 1);
        CHECK(data.y(0) == 2);
        CHECK(data.z(0) == 3);

        // the atom/block views the kernels consume must agree with the accessors
        auto a = data.atom(1);
        CHECK(a.x == 4);
        CHECK(a.y == 5);
        CHECK(a.z == 6);
        CHECK(a.ff == data.get_ff_type(1));

        auto b = data.block(0);
        CHECK(b.x[1] == 4);
        CHECK(b.ff[1] == data.get_ff_type(1));
        // the two species must map to different indices
        CHECK(data.get_ff_type(0) != data.get_ff_type(1));
    }
}

// The kernels take one broadcast atom plus one pointer per component of the opposing block,
// so the test inputs are built the same way. Note the form factor indices live in their own
// int32 array rather than being smuggled through a float lane.
template<std::size_t N>
struct Others {
    alignas(64) std::array<float, N> x{}, y{}, z{};
    alignas(64) std::array<int32_t, N> ff{};
    Block block() const {return Block{x.data(), y.data(), z.data(), ff.data()};}
};

template<std::size_t N>
Others<N> make_others(const std::array<std::pair<Vector3<double>, int32_t>, N>& in) {
    Others<N> o;
    for (std::size_t k = 0; k < N; ++k) {
        o.x[k] = static_cast<float>(in[k].first.x());
        o.y[k] = static_cast<float>(in[k].first.y());
        o.z[k] = static_cast<float>(in[k].first.z());
        o.ff[k] = in[k].second;
    }
    return o;
}

// SIMD backends may reorder output elements; sort by distance and compare as sets
template<std::size_t N>
void check_unordered(
    const std::array<float, N>& actual_dist,
    const std::array<int32_t, N>& actual_ff,
    std::vector<std::pair<double, int32_t>> expected,
    double tol)
{
    std::sort(expected.begin(), expected.end());
    std::array<std::size_t, N> idx;
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](auto a, auto b) { return actual_dist[a] < actual_dist[b]; });
    for (std::size_t k = 0; k < N; ++k) {
        CHECK_THAT(static_cast<double>(actual_dist[idx[k]]), Catch::Matchers::WithinAbs(expected[k].first, tol));
        CHECK(actual_ff[idx[k]] == expected[k].second);
    }
}

template<std::size_t N>
void check_unordered_rounded(
    const std::array<int32_t, N>& actual_dist,
    const std::array<int32_t, N>& actual_ff,
    std::vector<std::pair<int32_t, int32_t>> expected)
{
    std::sort(expected.begin(), expected.end());
    std::array<std::size_t, N> idx;
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](auto a, auto b) { return actual_dist[a] < actual_dist[b]; });
    for (std::size_t k = 0; k < N; ++k) {
        CHECK(actual_dist[idx[k]] == expected[k].first);
        CHECK(actual_ff[idx[k]] == expected[k].second);
    }
}

namespace {
    const Atom self{1, 1, 1, 2};

    // shared test geometry: atom n sits at (n+1, n+1, n+1) apart from the first
    const std::array<std::pair<Vector3<double>, int32_t>, 16> geometry = {{
        {{2, 1, 1}, 4},   {{2, 2, 2}, 8},     {{3, 3, 3}, 16},    {{4, 4, 4}, 32},
        {{5, 5, 5}, 64},  {{6, 6, 6}, 128},   {{7, 7, 7}, 15},    {{8, 8, 8}, 5},
        {{9, 9, 9}, 7},   {{10, 10, 10}, 11}, {{11, 11, 11}, 13}, {{12, 12, 12}, 17},
        {{13, 13, 13}, 19}, {{14, 14, 14}, 23}, {{15, 15, 15}, 29}, {{16, 16, 16}, 31}
    }};

    const std::array<double, 16> distances = {
        1.0, std::sqrt(3.0), std::sqrt(12.0), std::sqrt(27.0),
        std::sqrt(48.0), std::sqrt(75.0), std::sqrt(108.0), std::sqrt(147.0),
        std::sqrt(192.0), std::sqrt(243.0), std::sqrt(300.0), std::sqrt(363.0),
        std::sqrt(432.0), std::sqrt(507.0), std::sqrt(588.0), std::sqrt(675.0)
    };

    template<std::size_t N>
    Others<N> first_n() {
        std::array<std::pair<Vector3<double>, int32_t>, N> in;
        std::copy(geometry.begin(), geometry.begin()+N, in.begin());
        return make_others<N>(in);
    }

    template<std::size_t N>
    std::vector<std::pair<double, int32_t>> expected_n() {
        std::vector<std::pair<double, int32_t>> out;
        for (std::size_t k = 0; k < N; ++k) {out.emplace_back(distances[k], ff_bin_index(self.ff, geometry[k].second));}
        return out;
    }

    template<std::size_t N>
    std::vector<std::pair<int32_t, int32_t>> expected_n_rounded() {
        const double width = constants::axes::d_axis.width();
        std::vector<std::pair<int32_t, int32_t>> out;
        for (std::size_t k = 0; k < N; ++k) {
            out.emplace_back(static_cast<int32_t>(std::round(distances[k]/width)), ff_bin_index(self.ff, geometry[k].second));
        }
        return out;
    }
}

template<bool vbw>
void single_tests() {
    SECTION("single distance") {
        auto o = first_n<1>();
        auto result = evaluate<vbw>(self, o.block());
        CHECK(result.distance == 1);
        CHECK(result.ff_bin == ff_bin_index(2, 4));

        auto o2 = make_others<1>({{{{2, 2, 2}, 8}}});
        result = evaluate<vbw>(self, o2.block());
        CHECK_THAT(result.distance, Catch::Matchers::WithinAbs(std::sqrt(3), 1e-6));
        CHECK(result.ff_bin == ff_bin_index(2, 8));
    }
}

template<bool vbw>
void single_tests_rounded() {
    SECTION("single distance") {
        const double width = constants::axes::d_axis.width();
        auto o = first_n<1>();
        auto result = evaluate_rounded<vbw>(self, o.block());
        CHECK(result.distance == std::round(1./width));
        CHECK(result.ff_bin == ff_bin_index(2, 4));

        auto o2 = make_others<1>({{{{2, 2, 2}, 8}}});
        result = evaluate_rounded<vbw>(self, o2.block());
        CHECK(result.distance == std::round(std::sqrt(3)/width));
        CHECK(result.ff_bin == ff_bin_index(2, 8));
    }
}

template<std::size_t N, typename F>
void block_tests(F&& evaluate_block, double tol) {
    auto o = first_n<N>();
    auto result = evaluate_block(self, o.block());
    check_unordered<N>(result.distances, result.ff_bins, expected_n<N>(), tol);
}

template<std::size_t N, typename F>
void block_tests_rounded(F&& evaluate_block) {
    auto o = first_n<N>();
    auto result = evaluate_block(self, o.block());
    check_unordered_rounded<N>(result.distances, result.ff_bins, expected_n_rounded<N>());
}

template<bool vbw>
void run_tests() {
    SECTION("scalar") {
        single_tests<vbw>();
        single_tests_rounded<vbw>();
        block_tests<4>([](Atom s, Block b) {QuadEvaluatedResult r; evaluate_N_scalar<vbw, 4>(s, b, r); return r;}, 1e-6);
        block_tests<8>([](Atom s, Block b) {OctoEvaluatedResult r; evaluate_N_scalar<vbw, 8>(s, b, r); return r;}, 1e-5);
        block_tests<16>([](Atom s, Block b) {HexaEvaluatedResult r; evaluate_N_scalar<vbw, 16>(s, b, r); return r;}, 1e-3);
        block_tests_rounded<4>([](Atom s, Block b) {QuadEvaluatedResultRounded r; evaluate_rounded_N_scalar<vbw, 4>(s, b, r); return r;});
        block_tests_rounded<8>([](Atom s, Block b) {OctoEvaluatedResultRounded r; evaluate_rounded_N_scalar<vbw, 8>(s, b, r); return r;});
        block_tests_rounded<16>([](Atom s, Block b) {HexaEvaluatedResultRounded r; evaluate_rounded_N_scalar<vbw, 16>(s, b, r); return r;});
    }

    #if defined AUSAXS_USE_SSE2
        SECTION("sse") {
            block_tests<4>([](Atom s, Block b) {
                QuadEvaluatedResult r;
                evaluate_4_sse_into<vbw>(s, b, r.distances.data(), r.distance_bins.data(), r.ff_bins.data());
                return r;
            }, 1e-6);
            block_tests_rounded<4>([](Atom s, Block b) {
                QuadEvaluatedResultRounded r;
                evaluate_rounded_4_sse_into<vbw>(s, b, r.distances.data(), r.ff_bins.data());
                return r;
            });
        }
    #endif

    #if defined AUSAXS_USE_AVX2
        SECTION("avx") {
            block_tests<8>([](Atom s, Block b) {
                OctoEvaluatedResult r;
                evaluate_8_avx_into<vbw>(s, b, r.distances.data(), r.distance_bins.data(), r.ff_bins.data());
                return r;
            }, 1e-5);
            block_tests_rounded<8>([](Atom s, Block b) {
                OctoEvaluatedResultRounded r;
                evaluate_rounded_8_avx_into<vbw>(s, b, r.distances.data(), r.ff_bins.data());
                return r;
            });
        }
    #endif

    #if defined AUSAXS_USE_AVX512
        SECTION("avx512") {
            block_tests<16>([](Atom s, Block b) {
                HexaEvaluatedResult r;
                evaluate_16_avx512_into<vbw>(s, b, r.distances.data(), r.distance_bins.data(), r.ff_bins.data());
                return r;
            }, 1e-3);
            block_tests_rounded<16>([](Atom s, Block b) {
                HexaEvaluatedResultRounded r;
                evaluate_rounded_16_avx512_into<vbw>(s, b, r.distances.data(), r.ff_bins.data());
                return r;
            });
        }
    #endif

    SECTION("dispatch") {
        block_tests<4>([](Atom s, Block b) {return evaluate_4<vbw>(s, b);}, 1e-6);
        block_tests<8>([](Atom s, Block b) {return evaluate_8<vbw>(s, b);}, 1e-5);
        block_tests<16>([](Atom s, Block b) {return evaluate_16<vbw>(s, b);}, 1e-3);
        block_tests_rounded<4>([](Atom s, Block b) {return evaluate_rounded_4<vbw>(s, b);});
        block_tests_rounded<8>([](Atom s, Block b) {return evaluate_rounded_8<vbw>(s, b);});
        block_tests_rounded<16>([](Atom s, Block b) {return evaluate_rounded_16<vbw>(s, b);});
    }
}

TEST_CASE("CompactCoordinatesXYZFF<vbw>::evaluate") {
    SECTION("variable bin width") {
        run_tests<true>();
    }
    SECTION("fixed bin width") {
        run_tests<false>();
    }
}
