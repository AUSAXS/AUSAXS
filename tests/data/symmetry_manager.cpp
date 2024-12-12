#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <data/SymmetryManager.h>
#include <data/Body.h>
#include <data/Molecule.h>
#include <data/Symmetry.h>
#include <data/record/Atom.h>
#include <data/record/Water.h>
#include <settings/All.h>

#include "hist/hist_test_helper.h"

using namespace ausaxs;
using namespace ausaxs::data;

TEST_CASE("SymmetryManager: translations") {
    settings::molecule::implicit_hydrogens = false;

    SECTION("one body with one atom") {
        record::Atom a(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1);
        a.set_effective_charge(1);
        Molecule m({a});

        SECTION("no copies") {
            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();
            REQUIRE(1 <= h.size());
            CHECK(h[0] == 1);
            for (int i = 1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }

        SECTION("one copy") {
            m.get_body(0).add_symmetry({Vector3<double>(1, 0, 0)});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin = std::round(1*constants::axes::d_inv_width);
            REQUIRE(bin < static_cast<int>(h.size()));
            CHECK(h[0] == 2);
            for (int i = 1; i < bin; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[std::round(1*constants::axes::d_inv_width)] == 2);
            for (int i = bin+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }

        SECTION("two copies") {
            m.get_body(0).add_symmetry({Vector3<double>(-1, 0, 0)});
            m.get_body(0).add_symmetry({Vector3<double>( 1, 0, 0)});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(2*constants::axes::d_inv_width);
            REQUIRE(bin2 < static_cast<int>(h.size()));
            CHECK(h[0] == 3);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 4);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 2);
            for (int i = bin2+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }

        SECTION("four copies") {
            m.get_body(0).add_symmetry({Vector3<double>(-1, 0, 0)});
            m.get_body(0).add_symmetry({Vector3<double>( 1, 0, 0)});
            m.get_body(0).add_symmetry({Vector3<double>( 0,-1, 0)});
            m.get_body(0).add_symmetry({Vector3<double>( 0, 1, 0)});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(std::sqrt(2)*constants::axes::d_inv_width);
            int bin3 = std::round(2*constants::axes::d_inv_width);
            REQUIRE(bin3 < static_cast<int>(h.size()));
            CHECK(h[0] == 5);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 8);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 8);
            for (int i = bin2+1; i < bin3; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin3] == 4);
            for (int i = bin3+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }
    }

    SECTION("two atoms") {
        record::Atom a1(Vector3<double>(1, 0, 0), 1, constants::atom_t::C, "C", 1);
        record::Atom a2(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1);        
        a1.set_effective_charge(1);
        a2.set_effective_charge(1);
        Molecule m({a1, a2});

        SECTION("no copies") {
            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            REQUIRE(bin1 < static_cast<int>(h.size()));
            CHECK(h[0] == 2);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 2);
            for (int i = bin1+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }

        SECTION("one copy") {
            m.get_body(0).add_symmetry({Vector3<double>(0, 1, 0)});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(std::sqrt(2)*constants::axes::d_inv_width);
            REQUIRE(bin1 < static_cast<int>(h.size()));
            CHECK(h[0] == 4);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 8);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 4);
            for (int i = bin2+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }

        SECTION("two copies") {
            m.get_body(0).add_symmetry({Vector3<double>(0, 1, 0)});
            m.get_body(0).add_symmetry({Vector3<double>(0, 2, 0)});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(std::sqrt(2)*constants::axes::d_inv_width);
            int bin3 = std::round(2*constants::axes::d_inv_width);
            int bin4 = std::round(std::sqrt(5)*constants::axes::d_inv_width);
            REQUIRE(bin4 < static_cast<int>(h.size()));
            CHECK(h[0] == 6);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 14);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 8);
            for (int i = bin2+1; i < bin3; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin3] == 4);
            for (int i = bin3+1; i < bin4; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin4] == 4);
            for (int i = bin4+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }
    }

    SECTION("two bodies with one atom") {
        record::Atom a1(Vector3<double>(1, 0, 0), 1, constants::atom_t::C, "C", 1);
        record::Atom a2(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1);        
        a1.set_effective_charge(1);
        a2.set_effective_charge(1);
        Body b1({a1});
        Body b2({a2});
        Molecule m({b1, b2});

        SECTION("no copies") {
            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            REQUIRE(bin1 < static_cast<int>(h.size()));
            CHECK(h[0] == 2);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 2);
            for (int i = bin1+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }

        SECTION("one copy of body1") {
            m.get_body(0).add_symmetry({Vector3<double>(0, 1, 0)});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(std::sqrt(2)*constants::axes::d_inv_width);
            REQUIRE(bin1 < static_cast<int>(h.size()));
            CHECK(h[0] == 3);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 4);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 2);
            for (int i = bin2+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }

        SECTION("one copy of each #1") {
            m.get_body(0).add_symmetry({Vector3<double>(0, 1, 0)});
            m.get_body(1).add_symmetry({Vector3<double>(0, 1, 0)});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(std::sqrt(2)*constants::axes::d_inv_width);

            REQUIRE(bin2 < static_cast<int>(h.size()));
            CHECK(h[0] == 4);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 8);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 4);
            for (int i = bin2+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }

        SECTION("one copy of each #2") {
            m.get_body(0).add_symmetry({Vector3<double>( 1, 0, 0)});
            m.get_body(1).add_symmetry({Vector3<double>(-1, 0, 0)});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(2*constants::axes::d_inv_width);
            int bin3 = std::round(3*constants::axes::d_inv_width);

            REQUIRE(bin3 < static_cast<int>(h.size()));
            CHECK(h[0] == 4);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 6);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 4);
            for (int i = bin2+1; i < bin3; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin3] == 2);
            for (int i = bin3+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }        
    }

    SECTION("one body with waters") {
        record::Atom  a1(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1);
        record::Water w1(Vector3<double>(1, 0, 0), 1, constants::atom_t::O, "O", 1);
        a1.set_effective_charge(1);
        w1.set_effective_charge(1);
        Body b1({a1}, {w1});
        Molecule m({b1});

        SECTION("no copies") {
            hist::detail::SymmetryManager sm;
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
            m.get_body(0).add_symmetry({Vector3<double>(0, 1, 0)});

            hist::detail::SymmetryManager sm;
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
        m.get_body(0).get_waters() = std::move(m.get_waters());
        m.clear_hydration();

        data::detail::Symmetry s{{10, 0, 0}}; 
        SECTION("single copy") {
            m.get_body(0).add_symmetry(s);

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<true>(m);

            // manually perform the transformation for comparison
            data::Molecule m_copy({m.get_body(0)});
            data::Body& b_copy = m_copy.get_body(0);
            std::vector<record::Atom> a_copy = b_copy.get_atoms();
            for (auto& a : b_copy.get_atoms()) {
                a.get_coordinates() += s.translate;
                a_copy.push_back(a);
            }
            std::vector<record::Water> w_copy = b_copy.get_waters();
            for (auto& w : b_copy.get_waters()) {
                w.get_coordinates() += s.translate;
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
        record::Atom a(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1);
        a.set_effective_charge(1);
        Molecule m({a});

        SECTION("two repeats") {
            m.get_body(0).add_symmetry({Vector3<double>(1, 0, 0), Vector3<double>(0, 0, 0), 0, Vector3<double>(0, 0, 0), 2});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(2*constants::axes::d_inv_width);
            REQUIRE(bin2 < static_cast<int>(h.size()));
            CHECK(h[0] == 3);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 4);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 2);
            for (int i = bin2+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }

        SECTION("three repeats") {
            m.get_body(0).add_symmetry({Vector3<double>(1, 0, 0), Vector3<double>(0, 0, 0), 0, Vector3<double>(0, 0, 0), 3});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(2*constants::axes::d_inv_width);
            int bin3 = std::round(3*constants::axes::d_inv_width);
            REQUIRE(bin2 < static_cast<int>(h.size()));
            CHECK(h[0] == 4);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 6);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 4);
            for (int i = bin2+1; i < bin3; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin3] == 2);
            for (int i = bin3+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }
    }

    SECTION("two bodies") {
        record::Atom a1(Vector3<double>(1, 0, 0), 1, constants::atom_t::C, "C", 1);
        record::Atom a2(Vector3<double>(0, 0, 0), 1, constants::atom_t::C, "C", 1);        
        a1.set_effective_charge(1);
        a2.set_effective_charge(1);
        Body b1({a1});
        Body b2({a2});
        Molecule m({b1, b2});

        SECTION("two repeats") {
            m.get_body(0).add_symmetry({Vector3<double>(0, 1, 0), Vector3<double>(0, 0, 0), 0, Vector3<double>(0, 0, 0), 2});
            m.get_body(1).add_symmetry({Vector3<double>(0, 1, 0), Vector3<double>(0, 0, 0), 0, Vector3<double>(0, 0, 0), 2});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(std::sqrt(2)*constants::axes::d_inv_width);
            int bin3 = std::round(2*constants::axes::d_inv_width);
            int bin4 = std::round(std::sqrt(5)*constants::axes::d_inv_width);
            REQUIRE(bin3 < static_cast<int>(h.size()));
            CHECK(h[0] == 6);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 14);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 8);
            for (int i = bin2+1; i < bin3; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin3] == 4);
            for (int i = bin3+1; i < bin4; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin4] == 4);
            for (int i = bin4+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }

        SECTION("different repeats") {
            m.get_body(0).add_symmetry({Vector3<double>(0, 1, 0), Vector3<double>(0, 0, 0), 0, Vector3<double>(0, 0, 0), 1});
            m.get_body(1).add_symmetry({Vector3<double>(0, 1, 0), Vector3<double>(0, 0, 0), 0, Vector3<double>(0, 0, 0), 2});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<false>(m)->get_total_counts();

            int bin1 = std::round(1*constants::axes::d_inv_width);
            int bin2 = std::round(std::sqrt(2)*constants::axes::d_inv_width);
            int bin3 = std::round(2*constants::axes::d_inv_width);
            int bin4 = std::round(std::sqrt(5)*constants::axes::d_inv_width);
            REQUIRE(bin3 < static_cast<int>(h.size()));
            CHECK(h[0] == 5);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 10);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 6);
            for (int i = bin2+1; i < bin3; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin3] == 2);
            for (int i = bin3+1; i < bin4; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin4] == 2);
            for (int i = bin4+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }
    }
}

TEST_CASE("SymmetryManager: rotations") {
    settings::molecule::implicit_hydrogens = false;
    settings::molecule::center = false;

    SECTION("one body with one atom") {
        record::Atom a(Vector3<double>(1, 0, 0), 1, constants::atom_t::C, "C", 1);
        a.set_effective_charge(1);
        Molecule m({a});

        SECTION("one copy") {
            m.get_body(0).add_symmetry({{0., 1., 0.}, std::numbers::pi/2});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<true>(m)->get_total_counts();

            int bin1 = std::round(std::sqrt(2)*constants::axes::d_inv_width);
            REQUIRE(bin1 < static_cast<int>(h.size()));
            CHECK(h[0] == 2);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 2);
            for (int i = bin1+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }

        SECTION("three copies") {
            m.get_body(0).add_symmetry({{0., 0., 0.}, {0., 1., 0.}, std::numbers::pi/2, {0., 0., 0.}, 3});

            hist::detail::SymmetryManager sm;
            auto h = sm.calculate<true>(m)->get_total_counts();

            int bin1 = std::round(std::sqrt(2)*constants::axes::d_inv_width);
            int bin2 = std::round(2*constants::axes::d_inv_width);
            REQUIRE(bin1 < static_cast<int>(h.size()));
            CHECK(h[0] == 4);
            for (int i = 1; i < bin1; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin1] == 8);
            for (int i = bin1+1; i < bin2; ++i) {
                CHECK(h[i] == 0);
            }
            CHECK(h[bin2] == 4);
            for (int i = bin2+1; i < static_cast<int>(h.size()); ++i) {
                CHECK(h[i] == 0);
            }
        }
    }
}