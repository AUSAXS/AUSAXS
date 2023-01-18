#pragma once

#include <mini/Minimizer.h>

#include <concepts>

namespace mini {
    struct column_vector;
 
    template<mini::type algo>
    class dlibMinimizer : public Minimizer {
        public:
            dlibMinimizer();

            dlibMinimizer(std::function<double(std::vector<double>)> function, std::vector<Parameter> param = {});

            dlibMinimizer(std::function<double(double)> function, Parameter param = Parameter());

            ~dlibMinimizer() override;

            void max_evals(unsigned int max_evals);

        private: 
            unsigned int _max_evals = 100;

            /**
             * @brief Perform the minimization.
             */
            Result minimize_override() override;
    };
}