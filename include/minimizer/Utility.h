#pragma once

#include <utility/Limit.h>

#include <optional>

namespace mini {
    /**
     * @brief A representation of a parameter.
     */
    struct Parameter {
        /**
         * @brief Create a Parameter with a guess value and bounds.
         * 
         * @param name The name of the parameter.
         * @param guess The guess value.
         * @param bounds The bounds. 
         */
        Parameter(std::string name, double guess = 0, Limit bounds = {0, 0});

        /**
         * @brief Create a Parameter without a guess value.
         * 
         * @param name The name of the parameter.
         * @param bounds The bounds.
         */
        Parameter(std::string name, Limit bounds);

        /**
         * @brief Check if this parameter is bounded.
         */
        [[nodiscard]]
        bool has_bounds() const noexcept;

        /**
         * @brief Check if this parameter has a guess value.
         */
        [[nodiscard]]
        bool has_guess() const noexcept;

        std::string name;            // The name of this parameter.
        std::optional<double> guess; // The guess value.
        std::optional<Limit> bounds; // The bounds of this parameter. 
    };

    struct FittedParameter {
        /**
         * @brief Create a FittedParameter with asymmetric errors.
         * 
         * @param name Name of this parameter. 
         * @param val Optimal value of this parameter.
         * @param error Asymmetrical errors of this parameter.
         */
        FittedParameter(std::string name, double val, Limit error);

        /**
         * @brief Create a FittedParameter with symmetric errors.
         * 
         * @param name Name of this parameter. 
         * @param val Optimal value of this parameter.
         * @param error Symmetrical errors of this parameter.
         */
        FittedParameter(std::string name, double val, double error);

        /**
         * @brief Create a FittedParameter from a Parameter with asymmetric errors.
         * 
         * @param param The fitted parameter.
         * @param val Optimal value of this parameter.
         * @param error Asymmetrical errors of this parameter.
         */
        FittedParameter(const Parameter& param, double val, Limit error);

        /**
         * @brief Create a FittedParameter from a Parameter with symmetric errors.
         * 
         * @param param The fitted parameter.
         * @param val Optimal value of this parameter.
         * @param error Symmetrical errors of this parameter.
         */
        FittedParameter(const Parameter& param, double val, double error);

        std::string name; // The name of this parameter.
        double val;       // The optimal value of this parameter.
        Limit error;      // The error on this parameter. 
    };

    struct Result {
        /**
         * @brief Default constructor.
         */
        Result() {}

        /**
         * @brief Construct a Result. 
         * 
         * @param params The fitted parameters.
         * @param fval The function value.
         */
        Result(const FittedParameter& param, double fval);

        /**
         * @brief Construct a Result. 
         * 
         * @param params The fitted parameters.
         * @param fval The function value.
         */
        Result(const std::vector<FittedParameter>& params, double fval);

        FittedParameter get_parameter(std::string name) const;

        std::vector<FittedParameter> parameters;
        double fval;
    };

    struct Evaluation {
        Evaluation(std::vector<double> vals, double fval);

        std::vector<double> vals;
        double fval;
    };
}