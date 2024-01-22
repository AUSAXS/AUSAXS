#pragma once

#include <mini/Minimizer.h>
#include <mini/detail/Parameter.h>

#include <concepts>

#if defined(DLIB_AVAILABLE)
    namespace mini {
        struct column_vector;
    
        template<mini::type algo>
        class dlibMinimizer : public Minimizer {
            public:
                dlibMinimizer();

                dlibMinimizer(std::function<double(std::vector<double>)> function, std::vector<Parameter> param = {});

                dlibMinimizer(std::function<double(double)> function, const Parameter& param = Parameter());

                ~dlibMinimizer() override;

            private: 
                /**
                 * @brief Perform the minimization.
                 */
                Result minimize_override() override;
        };
    }
#endif