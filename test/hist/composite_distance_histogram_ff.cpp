#include <catch2/catch_test_macros.hpp>

#include <hist/CompositeDistanceHistogramFF.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/HistogramManagerMT.h>
#include <hist/HistogramManagerMTFF.h>
#include <hist/Histogram.h>
#include <hist/detail/FormFactor.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Atom.h>
#include <data/Water.h>
#include <io/ExistingFile.h>
#include <utility/Utility.h>
#include <settings/All.h>

bool compare_hist(Vector<double> p1, Vector<double> p2) {
    unsigned int pmin = std::min(p1.size(), p2.size());
    for (unsigned int i = 0; i < pmin; i++) {
        if (!utility::approx(p1[i], p2[i], 1e-6, 0.02)) {
            std::cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << std::endl;
            return false;
        }
    }

    if (p1.size() < p2.size()) {
        for (unsigned int i = p1.size(); i < p2.size(); i++) {
            if (!utility::approx(p2[i], 0)) {
                std::cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << std::endl;
                return false;
            }
        }
    }

    else {
        for (unsigned int i = p2.size(); i < p1.size(); i++) {
            if (!utility::approx(p1[i], 0)) {
                std::cout << "Failed on index " << i << ". Values: " << p1[i] << ", " << p2[i] << std::endl;
                return false;
            }
        }
    }

    return true;
}

void set_unity_charge(Protein& protein) {
    // set the weights to 1 so we can analytically determine the result
    // waters
    for (auto& atom : protein.get_waters()) {
        atom.set_effective_charge(1);
    }
    // atoms
    for (auto& body : protein.get_bodies()) {
        for (auto& atom : body.get_atoms()) {
            atom.set_effective_charge(1);
        }
    }
}

TEST_CASE("CompositeDistanceHistogramFF::debye_transform") {
    settings::protein::use_effective_charge = false;
    std::vector<Atom> b1 = {Atom(Vector3<double>(-1, -1, -1), 1, "C", "C", 1),  Atom(Vector3<double>(-1, 1, -1), 1, "C", "C", 1)};
    std::vector<Atom> b2 = {Atom(Vector3<double>( 1, -1, -1), 1, "C", "C", 1),  Atom(Vector3<double>( 1, 1, -1), 1, "C", "C", 1)};
    std::vector<Atom> b3 = {Atom(Vector3<double>(-1, -1,  1), 1, "C", "C", 1),  Atom(Vector3<double>(-1, 1,  1), 1, "C", "C", 1)};
    std::vector<Atom> b4 = {Atom(Vector3<double>( 1, -1,  1), 1, "C", "C", 1),  Atom(Vector3<double>( 1, 1,  1), 1, "C", "C", 1)};
    std::vector<Body> a = {Body(b1), Body(b2), Body(b3), Body(b4)};
    Protein protein(a);

    set_unity_charge(protein);
    double Z = protein.get_excluded_volume()*constants::charge::density::water/8;
    protein.set_excluded_volume_scaling(1./Z);

    // check distances
    // calculation: 8 identical points. 
    //      each point has:
    //          1 line  of length 0
    //          3 lines of length 2
    //          3 lines of length sqrt(2*2^2) = sqrt(8) = 2.82
    //          1 line  of length sqrt(3*2^2) = sqrt(12) = 3.46
    std::vector<double> p_exp = {8, 0, 2*8*3, 8, 0, 0, 0, 0, 0, 0};
    std::vector<double> Iq_exp;
    {
        auto ff_carbon = hist::detail::FormFactorStorage::get_form_factor(hist::detail::form_factor_t::NEUTRAL_CARBON);
        auto ff_exv = hist::detail::FormFactorStorage::get_form_factor(hist::detail::form_factor_t::EXCLUDED_VOLUME);
        std::vector<double> q_axis = Axis(settings::axes::qmin, settings::axes::qmax, settings::axes::bins).as_vector();
        Iq_exp.resize(q_axis.size(), 0);

        for (unsigned int q = 0; q < q_axis.size(); ++q) {
            double dsum = 
                8 + 
                24*std::sin(q_axis[q]*2)/(q_axis[q]*2) + 
                24*std::sin(q_axis[q]*std::sqrt(8))/(q_axis[q]*std::sqrt(8)) + 
                8*std::sin(q_axis[q]*std::sqrt(12))/(q_axis[q]*std::sqrt(12));
            if (q==940) {
                std::cout << "DSUM[" << q << "] = " << dsum << std::endl;
                std::cout << "\t" << std::sin(q_axis[q]*2)/(q_axis[q]*2) << std::endl;
                std::cout << "\t" << std::sin(q_axis[q]*std::sqrt(8))/(q_axis[q]*std::sqrt(8)) << std::endl;
                std::cout << "\t" << std::sin(q_axis[q]*std::sqrt(12))/(q_axis[q]*std::sqrt(12)) << std::endl;
            }
            Iq_exp[q] += dsum*std::pow(ff_carbon.evaluate(q_axis[q]), 2);
            Iq_exp[q] -= 2*dsum*ff_carbon.evaluate(q_axis[q])*ff_exv.evaluate(q_axis[q]);
            Iq_exp[q] += dsum*std::pow(ff_exv.evaluate(q_axis[q]), 2);
        }
    }

    auto Iq = hist::HistogramManagerMTFF(&protein).calculate_all()->debye_transform();
    REQUIRE(compare_hist(Iq_exp, Iq.get_counts()));
}

// TEST_CASE("CompositeDistanceHistogramFF::debye_transform") {
//     Protein protein("test/files/2epe.pdb");
//     protein.generate_new_hydration();

//     // change all atoms to Cl since these have the default form factor of 1
//     for (auto& body : protein.get_bodies()) {
//         for (auto& atom : body.get_atoms()) {
//             atom.set_element("Cl");
//         }
//     }

//     auto hm_mt = hist::HistogramManagerMT(&protein).calculate_all();
//     auto hm_mt_ff = hist::HistogramManagerMTFF(&protein).calculate_all();
//     REQUIRE(compare_hist(hm_mt->p, hm_mt_ff->p));

//     auto hm_mt_debye = hm_mt->debye_transform();
//     auto hm_mt_ff_debye = hm_mt_ff->debye_transform();
//     REQUIRE(compare_hist(hm_mt_debye.p, hm_mt_ff_debye.p));
// }