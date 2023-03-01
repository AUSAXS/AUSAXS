#pragma once

#include <crystal/miller/MillerGenerationStrategy.h>

namespace crystal {
    class AllMillers : public MillerGenerationStrategy {
        public: 
            AllMillers(unsigned int h, unsigned int k, unsigned int l) : h(h), k(k), l(l) {}

            std::vector<Miller> generate() const override {
                std::vector<Miller> millers;
                millers.reserve((h + 1)*(2*k + 1)*(2*l + 1));
                for (int h = 0; h <= this->h; h++) {
                    for (int k = -this->k; k <= this->k; k++) {
                        for (int l = -this->l; l <= this->l; l++) {
                            millers.push_back(Miller(h, k, l));
                        }
                    }
                }
                return millers;
            }
        private: 
            int h, k, l;
    };
}