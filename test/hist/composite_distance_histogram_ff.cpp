#include <catch2/catch_test_macros.hpp>

#include <hist/CompositeDistanceHistogramFF.h>
#include <hist/CompositeDistanceHistogram.h>
#include <hist/HistogramManagerMT.h>
#include <hist/HistogramManagerMTFF.h>
#include <hist/Histogram.h>
#include <data/Protein.h>
#include <data/Body.h>
#include <data/Atom.h>
#include <io/ExistingFile.h>
#include <utility/Utility.h>

bool compare_hist(Vector<double> p1, Vector<double> p2) {
    unsigned int pmin = std::min(p1.size(), p2.size());
    for (unsigned int i = 0; i < pmin; i++) {
        if (!utility::approx(p1[i], p2[i])) {
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