#pragma once

#include <utility/Limit.h>

#include <optional>

namespace mini {
    /**
     * @brief A representation of a parameter.
     */
    struct Parameter {
        Parameter(std::string name, double guess = 0, Limit bounds = {0, 0}) : name(name), guess(guess), bounds(bounds) {}

        bool has_bounds() const noexcept {return bounds.has_value();}
        bool has_guess() const noexcept {return guess.has_value();}

        std::string name;
        std::optional<double> guess;
        std::optional<Limit> bounds;
    };

    struct FittedParameter {
        /**
         * @brief Create a FittedParameter with asymmetric errors.
         * 
         * @param name Name of this parameter. 
         * @param val Optimal value of this parameter.
         * @param error Asymmetrical errors of this parameter.
         */
        FittedParameter(std::string name, double val, Limit error) : name(name), val(val), error(error) {}

        /**
         * @brief Create a FittedParameter with symmetric errors.
         * 
         * @param name Name of this parameter. 
         * @param val Optimal value of this parameter.
         * @param error Symmetrical errors of this parameter.
         */
        FittedParameter(std::string name, double val, double error) : name(name), val(val), error({-error, +error}) {}

        /**
         * @brief Create a FittedParameter from a Parameter with asymmetric errors.
         * 
         * @param param The fitted parameter.
         * @param val Optimal value of this parameter.
         * @param error Asymmetrical errors of this parameter.
         */
        FittedParameter(const Parameter& param, double val, Limit error) : name(param.name), val(val), error(error) {}

        /**
         * @brief Create a FittedParameter from a Parameter with symmetric errors.
         * 
         * @param param The fitted parameter.
         * @param val Optimal value of this parameter.
         * @param error Symmetrical errors of this parameter.
         */
        FittedParameter(const Parameter& param, double val, double error) : name(param.name), val(val), error(-error, +error) {}

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

        std::vector<FittedParameter> parameters;
        double fval;
    };

    struct Evaluation {
        Evaluation(std::vector<double> vals, double fval) : vals(vals), fval(fval) {}

        std::vector<double> vals;
        double fval;
    };
}