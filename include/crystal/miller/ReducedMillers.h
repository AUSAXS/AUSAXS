#pragma once

#include <crystal/miller/MillerGenerationStrategy.h>

namespace crystal {
    class ReducedMillers : public MillerGenerationStrategy {
        public:
            ReducedMillers();
            ReducedMillers(unsigned int h, unsigned int k, unsigned int l);

            std::vector<Miller> generate() const override;

        protected:
            std::vector<Miller> generate_independent_bases(double limit = -1) const;

        private:
            int h, k, l;
    };
}