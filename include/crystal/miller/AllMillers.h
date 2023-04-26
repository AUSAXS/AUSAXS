#pragma once

#include <crystal/miller/MillerGenerationStrategy.h>
#include <settings/CrystalSettings.h>

namespace crystal {
    /**
     * @brief Generates all miller indices within the range specified by
     *          settings::crystal::h
     *          settings::crystal::k
     *          settings::crystal::l
     * 
     *        The maximum allowed length of the indices is specified by
     *          settings::crystal::max_q
     *        which guarantees that the indices spans a spherical volume.
     */
    class AllMillers : public MillerGenerationStrategy {
        public: 
            AllMillers(unsigned int h, unsigned int k, unsigned int l) : h(h), k(k), l(l) {}

            std::vector<Miller> generate() const override {
                std::vector<Miller> millers;
                for (int h = 0; h <= this->h; h++) {
                    for (int k = -this->k; k <= this->k; k++) {
                        for (int l = -this->l; l <= this->l; l++) {
                            Miller m(h, k, l);
                            if (m.length() < settings::crystal::max_q) {
                                millers.push_back(Miller(h, k, l));
                            }
                        }
                    }
                }
                return millers;
            }
        private: 
            int h, k, l;
    };
}