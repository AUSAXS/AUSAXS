#pragma once

#include <crystal/miller/MillerGenerationStrategy.h>

namespace crystal {
    class AllMillers : public MillerGenerationStrategy {
        public: 
            AllMillers(unsigned int max_h, unsigned int max_k, unsigned int max_l) : max_h(max_h), max_k(max_k), max_l(max_l) {}

            std::vector<Miller<>> generate() override {
                std::vector<Miller<>> millers;
                millers.reserve((max_h + 1)*(2*max_k + 1)*(2*max_l + 1));
                for (int h = 0; h <= max_h; h++) {
                    for (int k = -max_k; k <= max_k; k++) {
                        for (int l = -max_l; l <= max_l; l++) {
                            millers.push_back(Miller<>(h, k, l));
                        }
                    }
                }
                return millers;
            }
        private: 
            int max_h, max_k, max_l;
    };
}