#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <hist/detail/CompactCoordinates.h>
#include <hist/detail/data/CompactCoordinatesXYZW.h>
#include <constants/Constants.h>
#include <math/Vector3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

using namespace ausaxs;
using namespace hist::detail;
using namespace hist::detail::xyzw;

TEST_CASE("CompactCoordinates<vbw>: component storage") {
    SECTION("positions and weights round-trip") {
        std::vector<data::AtomFF> atoms = {
            data::AtomFF({1, 2, 3}, form_factor::form_factor_t::C),
            data::AtomFF({4, 5, 6}, form_factor::form_factor_t::O)
        };
        CompactCoordinates<false> data(atoms);
        REQUIRE(data.size() == 2);
        CHECK(data.x(0) == 1);
        CHECK(data.y(0) == 2);
        CHECK(data.z(0) == 3);
        CHECK(data.x(1) == 4);
        CHECK(data.y(1) == 5);
        CHECK(data.z(1) == 6);

        // the atom/block views the kernels consume must agree with the accessors
        auto a = data.atom(1);
        CHECK(a.x == 4);
        CHECK(a.y == 5);
        CHECK(a.z == 6);
        CHECK(a.w == data.get_weight(1));

        auto b = data.block(0);
        CHECK(b.x[1] == 4);
        CHECK(b.y[1] == 5);
        CHECK(b.z[1] == 6);
    }

    SECTION("shuffle_order permutes atoms without splitting them") {
        std::vector<data::AtomFF> atoms;
        for (int i = 0; i < 64; ++i) {
            atoms.emplace_back(Vector3<double>{double(i), 2.0*i, 3.0*i}, form_factor::form_factor_t::C);
        }
        CompactCoordinates<false> data(atoms);
        data.shuffle_order(12345);
        REQUIRE(data.size() == 64);
        // every atom must still satisfy y == 2x and z == 3x, i.e. the components were
        // permuted by one common permutation rather than independently
        std::vector<float> seen;
        for (unsigned int i = 0; i < data.size(); ++i) {
            CHECK_THAT(data.y(i), Catch::Matchers::WithinAbs(2*data.x(i), 1e-6));
            CHECK_THAT(data.z(i), Catch::Matchers::WithinAbs(3*data.x(i), 1e-6));
            seen.push_back(data.x(i));
        }
        // and the set of atoms is unchanged
        std::sort(seen.begin(), seen.end());
        for (int i = 0; i < 64; ++i) {CHECK_THAT(seen[i], Catch::Matchers::WithinAbs(i, 1e-6));}
    }
}

// The kernels take one broadcast atom plus one pointer per component of the opposing block,
// so the test inputs are built the same way.
template<std::size_t N>
struct Others {
    alignas(64) std::array<float, N> x{}, y{}, z{}, w{};
    Block block() const {return Block{x.data(), y.data(), z.data(), w.data()};}
};

template<std::size_t N>
Others<N> make_others(const std::array<std::pair<Vector3<double>, float>, N>& in) {
    Others<N> o;
    for (std::size_t k = 0; k < N; ++k) {
        o.x[k] = static_cast<float>(in[k].first.x());
        o.y[k] = static_cast<float>(in[k].first.y());
        o.z[k] = static_cast<float>(in[k].first.z());
        o.w[k] = in[k].second;
    }
    return o;
}

// SIMD backends may reorder output elements; sort by distance and compare as sets
template<std::size_t N>
void check_unordered(
    const std::array<float, N>& actual_dist,
    const std::array<float, N>& actual_wt,
    std::vector<std::pair<double, float>> expected,
    double tol)
{
    std::sort(expected.begin(), expected.end());
    std::array<std::size_t, N> idx;
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](auto a, auto b) { return actual_dist[a] < actual_dist[b]; });
    for (std::size_t k = 0; k < N; ++k) {
        CHECK_THAT(static_cast<double>(actual_dist[idx[k]]), Catch::Matchers::WithinAbs(expected[k].first, tol));
        CHECK(actual_wt[idx[k]] == expected[k].second);
    }
}

template<std::size_t N>
void check_unordered_rounded(
    const std::array<int32_t, N>& actual_dist,
    const std::array<float, N>& actual_wt,
    std::vector<std::pair<int32_t, float>> expected)
{
    std::sort(expected.begin(), expected.end());
    std::array<std::size_t, N> idx;
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](auto a, auto b) { return actual_dist[a] < actual_dist[b]; });
    for (std::size_t k = 0; k < N; ++k) {
        CHECK(actual_dist[idx[k]] == expected[k].first);
        CHECK(actual_wt[idx[k]] == expected[k].second);
    }
}

namespace {
    const Atom self{1, 1, 1, 2};

    // shared test geometry: atom n sits at (n+1, n+1, n+1) apart from the first
    const std::array<std::pair<Vector3<double>, float>, 16> geometry = {{
        {{2, 1, 1}, 4},   {{2, 2, 2}, 8},     {{3, 3, 3}, 16},    {{4, 4, 4}, 32},
        {{5, 5, 5}, 64},  {{6, 6, 6}, 128},   {{7, 7, 7}, 15},    {{8, 8, 8}, 5},
        {{9, 9, 9}, 7},   {{10, 10, 10}, 11}, {{11, 11, 11}, 13}, {{12, 12, 12}, 17},
        {{13, 13, 13}, 19}, {{14, 14, 14}, 23}, {{15, 15, 15}, 29}, {{16, 16, 16}, 31}
    }};

    // the squared distances from self to each of the above, and the product weights
    const std::vector<std::pair<double, float>> expected_16 = {
        {1.0, 8.0f}, {std::sqrt(3.0), 16.0f}, {std::sqrt(12.0), 32.0f}, {std::sqrt(27.0), 64.0f},
        {std::sqrt(48.0), 128.0f}, {std::sqrt(75.0), 256.0f}, {std::sqrt(108.0), 30.0f}, {std::sqrt(147.0), 10.0f},
        {std::sqrt(192.0), 14.0f}, {std::sqrt(243.0), 22.0f}, {std::sqrt(300.0), 26.0f}, {std::sqrt(363.0), 34.0f},
        {std::sqrt(432.0), 38.0f}, {std::sqrt(507.0), 46.0f}, {std::sqrt(588.0), 58.0f}, {std::sqrt(675.0), 62.0f}
    };

    template<std::size_t N>
    Others<N> first_n() {
        std::array<std::pair<Vector3<double>, float>, N> in;
        std::copy(geometry.begin(), geometry.begin()+N, in.begin());
        return make_others<N>(in);
    }

    template<std::size_t N>
    std::vector<std::pair<double, float>> expected_n() {
        return {expected_16.begin(), expected_16.begin()+N};
    }

    template<std::size_t N>
    std::vector<std::pair<int32_t, float>> expected_n_rounded() {
        const double width = constants::axes::d_axis.width();
        std::vector<std::pair<int32_t, float>> out;
        for (std::size_t k = 0; k < N; ++k) {
            out.emplace_back(static_cast<int32_t>(std::round(expected_16[k].first/width)), expected_16[k].second);
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
        CHECK(result.weight == 8);

        auto o2 = make_others<1>({{{{2, 2, 2}, 8}}});
        result = evaluate<vbw>(self, o2.block());
        CHECK_THAT(result.distance, Catch::Matchers::WithinAbs(std::sqrt(3), 1e-6));
        CHECK(result.weight == 16);
    }
}

template<bool vbw>
void single_tests_rounded() {
    SECTION("single distance") {
        const double width = constants::axes::d_axis.width();
        auto o = first_n<1>();
        auto result = evaluate_rounded<vbw>(self, o.block());
        CHECK(result.distance == std::round(1./width));
        CHECK(result.weight == 8);

        auto o2 = make_others<1>({{{{2, 2, 2}, 8}}});
        result = evaluate_rounded<vbw>(self, o2.block());
        CHECK(result.distance == std::round(std::sqrt(3)/width));
        CHECK(result.weight == 16);
    }
}

template<std::size_t N, typename F>
void block_tests(F&& evaluate_block, double tol) {
    auto o = first_n<N>();
    auto result = evaluate_block(self, o.block());
    check_unordered<N>(result.distances, result.weights, expected_n<N>(), tol);
}

template<std::size_t N, typename F>
void block_tests_rounded(F&& evaluate_block) {
    auto o = first_n<N>();
    auto result = evaluate_block(self, o.block());
    check_unordered_rounded<N>(result.distances, result.weights, expected_n_rounded<N>());
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
                evaluate_4_sse_into<vbw>(s, b, r.distances.data(), r.distance_bins.data(), r.weights.data());
                return r;
            }, 1e-6);
            block_tests_rounded<4>([](Atom s, Block b) {
                QuadEvaluatedResultRounded r;
                evaluate_rounded_4_sse_into<vbw>(s, b, r.distances.data(), r.weights.data());
                return r;
            });
        }
    #endif

    #if defined AUSAXS_USE_AVX2
        SECTION("avx") {
            block_tests<8>([](Atom s, Block b) {
                OctoEvaluatedResult r;
                evaluate_8_avx_into<vbw>(s, b, r.distances.data(), r.distance_bins.data(), r.weights.data());
                return r;
            }, 1e-5);
            block_tests_rounded<8>([](Atom s, Block b) {
                OctoEvaluatedResultRounded r;
                evaluate_rounded_8_avx_into<vbw>(s, b, r.distances.data(), r.weights.data());
                return r;
            });
        }
    #endif

    #if defined AUSAXS_USE_AVX512
        SECTION("avx512") {
            block_tests<16>([](Atom s, Block b) {
                HexaEvaluatedResult r;
                evaluate_16_avx512_into<vbw>(s, b, r.distances.data(), r.distance_bins.data(), r.weights.data());
                return r;
            }, 1e-3);
            block_tests_rounded<16>([](Atom s, Block b) {
                HexaEvaluatedResultRounded r;
                evaluate_rounded_16_avx512_into<vbw>(s, b, r.distances.data(), r.weights.data());
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

TEST_CASE("CompactCoordinatesXYZW<vbw>::evaluate") {
    SECTION("variable bin width") {
        run_tests<true>();
    }
    SECTION("fixed bin width") {
        run_tests<false>();
    }
}
