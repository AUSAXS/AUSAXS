// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <crystal/miller/AllMillers.h>
#include <crystal/miller/Miller.h>
#include <settings/CrystalSettings.h>
#include <crystal/Fval.h>

using namespace ausaxs::crystal;

AllMillers::AllMillers(unsigned int h, unsigned int k, unsigned int l) : h(h), k(k), l(l) {}

std::vector<Miller> AllMillers::generate() const {
    std::vector<Miller> millers;
    for (int h = 0; h <= this->h; h++) {
        for (int k = -this->k; k <= this->k; k++) {
            for (int l = -this->l; l <= this->l; l++) {
                Miller m(h, k, l);
                double q = crystal::Fval::Q(m).norm();
                if (q < settings::crystal::max_q) {
                    millers.push_back(Miller(h, k, l));
                }
            }
        }
    }
    return millers;
}