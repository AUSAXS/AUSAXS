#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/symmetry/SymmetryManagerMT.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/symmetry/Symmetry.h>
#include <hist/intensity_calculator/ICompositeDistanceHistogramExv.h>
#include <hist/distribution/Distribution1D.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"

#include <random>

using namespace ausaxs;
using namespace ausaxs::data;

struct RES {
    RES(double d, int v) : index(std::round(d*constants::axes::d_inv_width)), val(v) {}
    int index;
    int val;
};

void check_hist(const std::vector<double>& h, std::vector<RES> checks) {
    std::sort(checks.begin(), checks.end(), [](const RES& a, const RES& b) {return a.index < b.index;});
    std::vector<RES> tmp;
    for (int i = 0; i < static_cast<int>(checks.size()); ++i) {
        if (i == 0 || checks[i].index != checks[i-1].index) {
            tmp.push_back(checks[i]);
        } else {
            tmp.back().val += checks[i].val;
        }
    }
    checks = tmp;
    REQUIRE(checks.back().index < static_cast<int>(h.size()));
    int j = 0;
    for (int i = 0; i < static_cast<int>(h.size()); ++i) {
        if (i == checks[j].index) {
            if (h[i] != checks[j].val) {
                INFO("i = " << i << ", dist = " << i*constants::axes::d_axis.width());
                INFO("h[i] = " << h[i] << ", checks[j].val = " << checks[j].val);
                CHECK(false);
            }
            ++j;
        } else {
            if (h[i] != 0) {
                INFO("i = " << i << ", dist = " << i*constants::axes::d_axis.width());
                INFO("h[i] = " << h[i] << ", checks[j].val = " << 0);
                CHECK(false);
            }
        }
    }
    SUCCEED();
}

TEST_CASE("SymmetryManager: translations") {
    settings::molecule::implicit_hydrogens = false;

    SECTION("one body with one atom") {
        AtomFF a({0, 0, 0}, form_factor::form_factor_t::C);
        Molecule m({Body{std::vector{a}}});
        set_unity_charge(m);

        SECTION("no copies") {
            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {RES(0, 1)});
        }

        SECTION("one copy") {
            m.get_body(0).symmetry().add({Vector3<double>(1, 0, 0)});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 2), 
                RES(1, 2)
            });
        }

        SECTION("two copies") {
            m.get_body(0).symmetry().add({Vector3<double>(-1, 0, 0)});
            m.get_body(0).symmetry().add({Vector3<double>( 1, 0, 0)});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 3), 
                RES(1, 4), 
                RES(2, 2)
            });
        }

        SECTION("four copies") {
            m.get_body(0).symmetry().add({Vector3<double>(-1, 0, 0)});
            m.get_body(0).symmetry().add({Vector3<double>( 1, 0, 0)});
            m.get_body(0).symmetry().add({Vector3<double>( 0,-1, 0)});
            m.get_body(0).symmetry().add({Vector3<double>( 0, 1, 0)});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 5), 
                RES(1, 8), 
                RES(std::sqrt(2), 8), 
                RES(2, 4)
            });
        }
    }

    SECTION("two atoms") {
        AtomFF a1({1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);        
        Molecule m({Body{std::vector{a1, a2}}});
        set_unity_charge(m);

        SECTION("no copies") {
            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 2), 
                RES(1, 2)
            });
        }

        SECTION("one copy") {
            m.get_body(0).symmetry().add({Vector3<double>(0, 1, 0)});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 4), 
                RES(1, 8), 
                RES(std::sqrt(2), 4)
            });
        }

        SECTION("two copies") {
            m.get_body(0).symmetry().add({Vector3<double>(0, 1, 0)});
            m.get_body(0).symmetry().add({Vector3<double>(0, 2, 0)});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 6), 
                RES(1, 14), 
                RES(std::sqrt(2), 8), 
                RES(2, 4),
                RES(std::sqrt(5), 4)
            });
        }
    }

    SECTION("two bodies with one atom") {
        AtomFF a1({1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);        
        Molecule m({Body{std::vector{a1}}, Body{std::vector{a2}}});
        set_unity_charge(m);

        SECTION("no copies") {
            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 2), 
                RES(1, 2)
            });
        }

        SECTION("one copy of body1") {
            m.get_body(0).symmetry().add({Vector3<double>(0, 1, 0)});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 3), 
                RES(1, 4), 
                RES(std::sqrt(2), 2)
            });
        }

        SECTION("one copy of each #1") {
            m.get_body(0).symmetry().add({Vector3<double>(0, 1, 0)});
            m.get_body(1).symmetry().add({Vector3<double>(0, 1, 0)});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 4), 
                RES(1, 8), 
                RES(std::sqrt(2), 4)
            });
        }

        SECTION("one copy of each #2") {
            m.get_body(0).symmetry().add({Vector3<double>( 1, 0, 0)});
            m.get_body(1).symmetry().add({Vector3<double>(-1, 0, 0)});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 4), 
                RES(1, 6),
                RES(2, 4),
                RES(3, 2)
            });
        }
    }

    SECTION("one body with waters") {
        AtomFF a1({0, 0, 0}, form_factor::form_factor_t::C);
        Water  w1({1, 0, 0});
        Body b1(std::vector<AtomFF>{a1}, std::vector{w1});
        Molecule m({b1});
        set_unity_charge(m);

        SECTION("no copies") {
            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<true>(m);
            auto htot = h->get_total_counts();
            auto haa = h->get_aa_counts();
            auto haw = h->get_aw_counts();
            auto hww = h->get_ww_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            REQUIRE(bin1 < static_cast<int>(htot.size()));
            CHECK(htot[0] == 2);
            CHECK(haa[0] == 1);
            CHECK(haw[0] == 0);
            CHECK(hww[0] == 1);
            for (int i = 1; i < bin1; ++i) {
                CHECK(htot[i] == 0);
            }
            CHECK(htot[bin1] == 2);
            CHECK(haa[bin1] == 0);
            CHECK(haw[bin1] == 2);
            CHECK(hww[bin1] == 0);
            for (int i = bin1+1; i < static_cast<int>(htot.size()); ++i) {
                CHECK(htot[i] == 0);
            }
        }

        SECTION("one copy") {
            m.get_body(0).symmetry().add({Vector3<double>(0, 1, 0)});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<true>(m);
            auto htot = h->get_total_counts();
            auto haa = h->get_aa_counts();
            auto haw = h->get_aw_counts();
            auto hww = h->get_ww_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(std::sqrt(2)*constants::axes::d_inv_width);
            REQUIRE(bin1 < static_cast<int>(htot.size()));
            CHECK(htot[0] == 4);
            CHECK(haa[0] == 2);
            CHECK(haw[0] == 0);
            CHECK(hww[0] == 2);
            for (int i = 1; i < bin1; ++i) {
                CHECK(htot[i] == 0);
            }
            CHECK(htot[bin1] == 8);
            CHECK(haa[bin1] == 2);
            CHECK(haw[bin1] == 4);
            CHECK(hww[bin1] == 2);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(htot[i] == 0);
            }
            CHECK(htot[bin2] == 4);
            CHECK(haa[bin2] == 0);
            CHECK(haw[bin2] == 4);
            CHECK(hww[bin2] == 0);
            for (int i = bin2+1; i < static_cast<int>(htot.size()); ++i) {
                CHECK(htot[i] == 0);
            }
        }
    }

    SECTION("real data") {
        settings::hist::histogram_manager = settings::hist::HistogramManagerChoice::HistogramManagerMT;
        settings::general::verbose = false;

        data::Molecule m("tests/files/2epe.pdb");
        m.generate_new_hydration();
        m.get_body(0).get_waters() = m.get_waters();
        m.clear_hydration();
        set_unity_charge(m);

        symmetry::Symmetry s{{10, 0, 0}}; 
        SECTION("single copy") {
            m.get_body(0).symmetry().add(symmetry::Symmetry(s));

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<true>(m);

            // manually perform the transformation for comparison
            data::Molecule m_copy({m.get_body(0)});
            data::Body& b_copy = m_copy.get_body(0);
            std::vector<AtomFF> a_copy = b_copy.get_atoms();
            for (auto& a : b_copy.get_atoms()) {
                a.coordinates() += s.translate;
                a_copy.push_back(a);
            }
            std::vector<Water> w_copy = b_copy.get_waters();
            for (auto& w : b_copy.get_waters()) {
                w.coordinates() += s.translate;
                w_copy.push_back(w);
            }
            b_copy.get_atoms() = std::move(a_copy);
            m_copy.get_waters() = std::move(w_copy);
            // ###

            REQUIRE(m_copy.size_atom() == 2*m.size_atom());
            REQUIRE(m_copy.size_water() == 2*m.get_body(0).size_water());

            auto h2 = m_copy.get_histogram();
            CHECK(compare_hist_approx(h->get_aa_counts(), h2->get_aa_counts()));
            CHECK(compare_hist_approx(h->get_ww_counts(), h2->get_ww_counts()));
            CHECK(compare_hist_approx(h->get_aw_counts(), h2->get_aw_counts()));
            CHECK(compare_hist_approx(h->get_total_counts(), h2->get_total_counts()));
        }
    }
}

TEST_CASE("SymmetryManager: repeating symmetries") {
    settings::molecule::implicit_hydrogens = false;

    SECTION("one body with one atom") {
        AtomFF a({0, 0, 0}, form_factor::form_factor_t::C);
        Molecule m({Body{std::vector{a}}});
        set_unity_charge(m);

        SECTION("two repeats") {
            m.get_body(0).symmetry().add({{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, 2});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 3), 
                RES(1, 4),
                RES(2, 2)
            });
        }

        SECTION("three repeats") {
            m.get_body(0).symmetry().add({{1, 0, 0}, {0, 0, 0}, {0, 0, 0}, 3});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 4), 
                RES(1, 6),
                RES(2, 4),
                RES(3, 2)
            });
        }
    }

    SECTION("two bodies") {
        AtomFF a1({1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({0, 0, 0}, form_factor::form_factor_t::C);        
        Molecule m({Body{std::vector{a1}}, Body{std::vector{a2}}});
        set_unity_charge(m);

        SECTION("two repeats") {
            m.get_body(0).symmetry().add({{0, 1, 0}, {0, 0, 0}, {0, 0, 0}, 2});
            m.get_body(1).symmetry().add({{0, 1, 0}, {0, 0, 0}, {0, 0, 0}, 2});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 6), 
                RES(1, 14),
                RES(std::sqrt(2), 8),
                RES(2, 4),
                RES(std::sqrt(5), 4)
            });
        }

        SECTION("different repeats") {
            m.get_body(0).symmetry().add({{0, 1, 0}, {0, 0, 0}, {0, 0, 0}, 1});
            m.get_body(1).symmetry().add({{0, 1, 0}, {0, 0, 0}, {0, 0, 0}, 2});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 5), 
                RES(1, 10),
                RES(std::sqrt(2), 6),
                RES(2, 2),
                RES(std::sqrt(5), 2)
            });
        }
    }
}

TEST_CASE("SymmetryManager: rotations") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;

    SECTION("one body with one atom") {
        AtomFF a({1, 0, 0}, form_factor::form_factor_t::C);
        Molecule m({Body{std::vector{a}}});
        set_unity_charge(m);

        SECTION("one copy") {
            m.get_body(0).symmetry().add({{0, 0, 0}, {0, std::numbers::pi/2, 0}});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<true>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 2), 
                RES(std::sqrt(2), 2)
            });
        }

        SECTION("three copies") {
            m.get_body(0).symmetry().add({{0, 0, 0}, {0, 0, 0}, {0, std::numbers::pi/2, 0}, 3});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<true>(m)->get_total_counts();
            check_hist(h, {
                RES(0, 4), 
                RES(std::sqrt(2), 8),
                RES(2, 4)
            });
        }
    }
}

TEST_CASE("SymmetryManager: multi-atom systems") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;

    SECTION("cross") {
        AtomFF a1({1, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a2({2, 0, 0}, form_factor::form_factor_t::C);
        AtomFF a3({3, 0, 0}, form_factor::form_factor_t::C);
        Molecule m({Body{std::vector{a1, a2, a3}}});
        set_unity_charge(m);

        // external rotate pi/2 around the y-axis and replicate thrice
        // this gives the structure
        //
        //         x
        //         x
        //         x
        //   x x x   x x x
        //         x
        //         x
        //         x
        //
        m.get_body(0).symmetry().add({{0, 0, 0}, {0, 0, 0}, {0, 0, std::numbers::pi/2}, 3});

        symmetry::SymmetryManagerMT sm;
        auto h = sm.calculate<true>(m)->get_total_counts();
        check_hist(h, {
            {0, 12},
            {1, 16},
            {2, 12},
            {3, 8},
            {4, 12},
            {5, 8},
            {6, 4},
            {std::sqrt(2), 8},
            {std::sqrt(5), 16},
            {std::sqrt(8), 8},
            {std::sqrt(10), 16},
            {std::sqrt(13), 16},
            {std::sqrt(18), 8}
        });
    }

    SECTION("helix") {
        AtomFF a1({1, 0, 0}, form_factor::form_factor_t::C);
        Molecule m({Body{std::vector{a1}}});
        set_unity_charge(m);
        m.get_body(0).symmetry().add({{0, 0, 1}, {0, 0, 0}, {0, 0, std::numbers::pi/2}, 4});

        symmetry::SymmetryManagerMT sm;
        auto h = sm.calculate<true>(m)->get_total_counts();

        std::vector<RES> checks = {
            {0, 5},
            {std::sqrt(3),  8}, // 1.73
            {std::sqrt(8),  6}, // 2.83
            {std::sqrt(11), 4}, // 3.32
            {4, 2},
        };
        check_hist(h, checks);
    }
}

TEST_CASE("SymmetryManager: random tests") {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> d(-10, 10);
    static std::uniform_real_distribution<> r(-2*std::numbers::pi, 2*std::numbers::pi);
    static std::uniform_int_distribution<> n(1, 10);

    int n_atoms = GENERATE(1, 3, 5, 7, 9);
    SECTION("single symmetry") {
        for (int i = 0; i < 10; ++i) {

            std::vector<AtomFF> atoms;
            for (int j = 0; j < n_atoms; ++j) {
                atoms.push_back(
                    AtomFF(
                        {d(gen), d(gen), d(gen)},
                        form_factor::form_factor_t::C
                    )
                );
            }

            Vector3<double> translate{d(gen), d(gen), d(gen)};
            Vector3<double> axis{d(gen), d(gen), d(gen)};
            Vector3<double> angles{r(gen), r(gen), r(gen)};
            int repeats = n(gen);

            Molecule m({Body{atoms}});
            set_unity_charge(m);
            m.get_body(0).symmetry().add({translate, axis, angles, repeats});

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<true>(m)->get_total_counts();

            auto b = m.get_body(0).symmetry().explicit_structure();
            auto m2 = Molecule({Body{std::move(b.atoms), std::move(b.waters)}});
            set_unity_charge(m2);
            auto h2 = m2.get_histogram()->get_total_counts();

            CHECK(compare_hist_approx(h, h2));
        }
    }

    SECTION("multiple symmetries") {
        int n_atoms = 5;
        for (int i = 0; i < 10; ++i) {

            std::vector<AtomFF> atoms;
            for (int j = 0; j < n_atoms; ++j) {
                atoms.push_back(
                    AtomFF(
                        {d(gen), d(gen), d(gen)},
                        form_factor::form_factor_t::C
                    )
                );
            }

            Molecule m({Body{atoms}});
            set_unity_charge(m);
            for (int j = 0; j < n(gen); ++j) {
                Vector3<double> translate{d(gen), d(gen), d(gen)};
                Vector3<double> axis{d(gen), d(gen), d(gen)};
                Vector3<double> angles{r(gen), r(gen), r(gen)};
                int repeats = n(gen);
                m.get_body(0).symmetry().add({translate, axis, angles, repeats});
            }

            symmetry::SymmetryManagerMT sm;
            auto h = sm.calculate<true>(m)->get_total_counts();

            auto b = m.get_body(0).symmetry().explicit_structure();
            auto m2 = Molecule({Body{std::move(b.atoms), std::move(b.waters)}});
            set_unity_charge(m2);
            auto h2 = m2.get_histogram()->get_total_counts();

            compare_hist_approx(h, h2);
        }
    }
}