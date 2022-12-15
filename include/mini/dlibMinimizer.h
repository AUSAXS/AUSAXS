#pragma once

#include <mini/Minimizer.h>
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

namespace mini {
    class dlibMinimizer : public Minimizer {
        typedef dlib::matrix<double,0,1> column_vector;
        public:
            dlibMinimizer();

            dlibMinimizer(std::function<double(std::vector<double>)> function, std::vector<Parameter> param = {});

            dlibMinimizer(std::function<double(double)> function, Parameter param = Parameter());

            ~dlibMinimizer() override;

        private: 
            std::function<double(column_vector)> dlib_fwrapper;

            /**
             * @brief Perform the minimization.
             */
            Result minimize_override() override;
            
            /**
             * @brief Prepare the minimizer for fitting. 
             *        This fills it with the added parameters and sets its function.
             */
            void prepare_minimizer();
    };
}