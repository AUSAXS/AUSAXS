#pragma once

#include <utility/Limit.h>

#include <optional>

namespace mini {
    /**
     * @brief A representation of a parameter.
     */
    struct Parameter {
        Parameter() {}

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

        /**
         * @brief Check if this parameter is named.
         */
        [[nodiscard]]
        bool has_name() const noexcept;

        /**
         * @brief Check if this parameter has been initialized.
         */
        [[nodiscard]]
        bool empty() const noexcept;

        /**
         * @brief Get a string representation of this parameter.
         */
        std::string to_string() const noexcept;

        /**
		 * @brief Stream output operator. 
		 * 
		 * Allows this object to easily be output to a given stream. 
		 */
		friend std::ostream& operator<<(std::ostream& os, const Parameter& param) {os << param.to_string(); return os;}

        std::string name;            // The name of this parameter.
        std::optional<double> guess; // The guess value.
        std::optional<Limit> bounds; // The bounds of this parameter. 
    };

    struct FittedParameter {
        /**
         * @brief Default constructor.
         */
        FittedParameter() {}

        /**
         * @brief Create a FittedParameter with asymmetric errors.
         * 
         * @param name Name of this parameter. 
         * @param value Optimal value of this parameter.
         * @param error Asymmetrical errors of this parameter.
         */
        FittedParameter(std::string name, double value, Limit error);

        /**
         * @brief Create a FittedParameter with symmetric errors.
         * 
         * @param name Name of this parameter. 
         * @param value Optimal value of this parameter.
         * @param error Symmetrical errors of this parameter.
         */
        FittedParameter(std::string name, double value, double error);

        /**
         * @brief Create a FittedParameter from a Parameter with asymmetric errors.
         * 
         * @param param The fitted parameter.
         * @param value Optimal value of this parameter.
         * @param error Asymmetrical errors of this parameter.
         */
        FittedParameter(const Parameter& param, double value, Limit error);

        /**
         * @brief Create a FittedParameter from a Parameter with symmetric errors.
         * 
         * @param param The fitted parameter.
         * @param value Optimal value of this parameter.
         * @param error Symmetrical errors of this parameter.
         */
        FittedParameter(const Parameter& param, double value, double error);

        /**
         * @brief Get a string representation of this parameter.
         */
        std::string to_string() const noexcept;

        /**
		 * @brief Stream output operator. 
		 * 
		 * Allows this object to easily be output to a given stream. 
		 */
		friend std::ostream& operator<<(std::ostream& os, const FittedParameter& param) {os << param.to_string(); return os;}

        std::string name; // The name of this parameter.
        double value;     // The optimal value of this parameter.
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
         * @param fevals The number of function evaluations.
         */
        Result(const FittedParameter& param, double fval, unsigned int fevals);

        /**
         * @brief Construct a Result. 
         * 
         * @param params The fitted parameters.
         * @param fval The function value.
         * @param fevals The number of function evaluations.
         */
        Result(const std::vector<FittedParameter>& params, double fval, unsigned int fevals);

        /**
         * @brief Get a parameter based on its name from this result.
         */
        FittedParameter get_parameter(std::string name) const;

        /**
         * @brief Get a parameter based on its index from this result.
         */
        FittedParameter get_parameter(unsigned int index) const;

        /**
         * @brief Add a parameter to this result.
         */
        void add_parameter(const FittedParameter& param);

        /**
         * @brief Get the number of parameters in this result.
         */
        size_t size() const noexcept;

        /**
         * @brief Get the number of parameters in this result.
         */
        size_t dim() const noexcept;

        std::vector<FittedParameter> parameters; // The fitted parameters
        double fval;                             // The minimum function value
        unsigned int fevals;                     // The number of function evaluations
        int status = 0;                          // The minimization status. Anything but 0 indicates an error with the fit. 
    };

    struct Evaluation {
        Evaluation(std::vector<double> vals, double fval);

        std::vector<double> vals;
        double fval;
    };
}