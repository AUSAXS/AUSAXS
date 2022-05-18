#pragma once

#include <minimizer/Minimizer.h>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

namespace mini {
    /**
     * @brief A wrapper class for the minimization algorithms provided by ROOT. 
     */
    class ROOTMinimizer : public Minimizer {
        public:
            /**
             * @brief Default constructor.
             */
            ROOTMinimizer() {}

            /**
             * @brief Constructor.
             * 
             * Construct a new minimizer using one of the algorithms provided by ROOT. 
             * 
             * @param package The package containing the algorithm.
             * @param algorithm The name of the algorithm. 
             */
            ROOTMinimizer(std::string package, std::string algorithm);

            /**
             * @brief Constructor.
             * 
             * Construct a new minimizer using one of the algorithms provided by ROOT. 
             * 
             * @param package The package containing the algorithm.
             * @param algorithm The name of the algorithm. 
             * @param function The function to minimize.
             * @param param A single parameter. 
             */
            ROOTMinimizer(std::string package, std::string algorithm, double(&function)(const double*), Parameter param = Parameter());

            /**
             * @brief Constructor.
             * 
             * Construct a new minimizer using one of the algorithms provided by ROOT. 
             * 
             * @param package The package containing the algorithm.
             * @param algorithm The name of the algorithm. 
             * @param function The function to minimize.
             * @param param A single parameter. 
             */
            ROOTMinimizer(std::string package, std::string algorithm, std::function<double(const double*)> function, Parameter param = Parameter());

            /**
             * @brief Constructor.
             * 
             * Construct a new minimizer using one of the algorithms provided by ROOT. 
             * 
             * @param package The package containing the algorithm.
             * @param algorithm The name of the algorithm. 
             * @param function The function to minimize. 
             * @param params A list of parameters.
             */
            ROOTMinimizer(std::string package, std::string algorithm, double(&function)(const double*), std::vector<Parameter> params = {});

            /**
             * @brief Constructor.
             * 
             * Construct a new minimizer using one of the algorithms provided by ROOT. 
             * 
             * @param package The package containing the algorithm.
             * @param algorithm The name of the algorithm. 
             * @param function The function to minimize. 
             * @param params A list of parameters.
             */
            ROOTMinimizer(std::string package, std::string algorithm, std::function<double(const double*)> function, std::vector<Parameter> params = {});

            /**
             * @brief Perform the minimization.
             */
            Result minimize() override;

            /**
             * @brief Generate a landscape of the function.
             */
            Dataset landscape(unsigned int evals = 100) const override;

        private: 
            ROOT::Math::Minimizer* mini;
            ROOT::Math::Functor functor;
            
            /**
             * @brief Create a ROOT minimizer.
             */
            void create_minimizer(std::string package, std::string algorithm);

            /**
             * @brief Prepare the minimizer for fitting. 
             *        This fills it with the added parameters and sets its function.
             */
            void prepare_minimizer();
    };
}