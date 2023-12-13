#include <crystal/miller/ReducedMillers.h>
#include <crystal/miller/Miller.h>
#include <crystal/Fval.h>
#include <settings/CrystalSettings.h>

#include <iostream>

using namespace crystal;

ReducedMillers::ReducedMillers() : h(0), k(0), l(0) {}

ReducedMillers::ReducedMillers(unsigned int h, unsigned int k, unsigned int l) : h(h), k(k), l(l) {}

std::vector<Miller> ReducedMillers::generate() const {
    std::vector<Miller> bases = generate_independent_bases();
    std::cout << "Generated " << bases.size() << " independent bases" << std::endl;
    std::vector<Miller> millers;

    // now generate all millers indices
    // we can do this by multiplying the base pairs with integers
    // int abs_h = std::abs(h), abs_k = std::abs(k), abs_l = std::abs(l);
    for (const auto& base : bases) {
        int multiplier = 1;
        while (multiplier++ < 10000) { // hard limit to prevent infinite loop
            Miller miller(base.h*multiplier, base.k*multiplier, base.l*multiplier);
            double q = crystal::Fval::Q(miller).norm();
            if (q > settings::crystal::max_q) {break;}
            millers.emplace_back(miller);
        }
    }

    return millers;
}

std::vector<Miller> ReducedMillers::generate_independent_bases(double limit) const {
    if (limit < 0) {limit = settings::crystal::reduced::basis_q;}
    int limit2 = static_cast<int>(std::round(limit*limit));
    int ilimit = static_cast<int>(std::round(limit));

    std::vector<Miller> millers;
    for (int h = -ilimit; h <= ilimit; h++) {
        for (int k = -ilimit; k <= ilimit; k++) {
            for (int l = -ilimit; l <= ilimit; l++) {
                if (h*h + k*k + l*l <= limit2) {
                    millers.emplace_back(h, k, l);
                }
            }
        }
    }

    // filter out Friedel symmetry equivalent pairs
    std::vector<Miller> friedel_independent;
    for (unsigned int i = 0; i < millers.size(); i++) {
        bool is_independent = true;
        for (unsigned int j = 0; j < i; j++) {
            if (millers[i].friedel_equivalent(millers[j])) {
                is_independent = false;
                break;
            }
        }
        if (is_independent) {friedel_independent.push_back(millers[i]);}
    }

    // sort the millers indices by their length
    std::sort(friedel_independent.begin(), friedel_independent.end(), [](const Miller& a, const Miller& b) {
        return a.length2() < b.length2();
    });

    // since the indices are sorted, we only check the previous indices for linear independence
    // thus we will keep the shortest bases and discard the longer ones
    std::vector<Miller> linearly_independent;
    for (unsigned int i = 1; i < friedel_independent.size(); i++) {
        bool is_independent = true;
        for (unsigned int j = 0; j < i; j++) {
            for (unsigned int k = 1; k <= limit; k++) {
                if (friedel_independent[i] == friedel_independent[j]*k) {
                    is_independent = false;
                    break;
                }
                if (!is_independent) {break;}
            }
        }
        if (is_independent) {linearly_independent.push_back(friedel_independent[i]);}
    }

    return linearly_independent;
}