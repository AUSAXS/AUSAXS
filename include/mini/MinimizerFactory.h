#pragma once

#include <Symbols.h>
#include <mini/dlibMinimizer.h>
#include <mini/Golden.h>
#include <mini/MinimumExplorer.h>
#include <mini/Scan.h>
#include <mini/LimitedScan.h>

#include <memory>
#include <functional>
#include <concepts>

namespace mini {
    namespace detail {
        inline std::shared_ptr<Minimizer> create_minimizer(type t) {
            switch (t) {
            case type::DLIB_GLOBAL:
                return std::make_shared<dlibMinimizer<type::DLIB_GLOBAL>>();
            case type::BFGS:
                return std::make_shared<dlibMinimizer<type::BFGS>>();
            case type::GOLDEN:
                return std::make_shared<Golden>();
            case type::MINIMUM_EXPLORER:
                return std::make_shared<MinimumExplorer>();
            case type::SCAN:
                return std::make_shared<Scan>();
            case type::LIMITED_SCAN:
                return std::make_shared<LimitedScan>();
            default:
                throw except::invalid_argument("all::create_minimizer: Unknown minimizer type.");
            }
        }
    }

    /**
     * @brief Create a new minimizer.
     *
     * @param t The algorithm to use.
     * @param func The function to minimize.
     */
    [[maybe_unused]] inline std::shared_ptr<Minimizer> create_minimizer(type t, const std::function<double(std::vector<double>)>& func) {
        auto minimizer = detail::create_minimizer(t);
        minimizer->set_function(func);
        return minimizer;
    }

    /**
     * @brief Create a new minimizer.
     *
     * @param t The algorithm to use.
     * @param func The function to minimize.
     * @param param The first parameter.
     */
    [[maybe_unused]] inline std::shared_ptr<Minimizer> create_minimizer(type t, const std::function<double(std::vector<double>)>& func, const Parameter& param) {
        auto minimizer = create_minimizer(t, func);
        minimizer->add_parameter(param);
        return minimizer;
    }

    /**
     * @brief Create a new minimizer.
     *
     * @param t The algorithm to use.
     * @param func The function to minimize.
     * @param param The parameter list.
     */
    [[maybe_unused]] inline std::shared_ptr<Minimizer> create_minimizer(type t, const std::function<double(std::vector<double>)>& func, const std::vector<Parameter>& param) {
        auto minimizer = create_minimizer(t, func);
        std::for_each(param.begin(), param.end(), [&](const Parameter& p) { minimizer->add_parameter(p); });
        return minimizer;
    }

    /**
     * @brief Create a new minimizer.
     *
     * @param t The algorithm to use.
     * @param func The function to minimize.
     * @param param The first parameter.
     * @param evals The number of evaluations to perform. Not supported by all minimizers.
     */
    [[maybe_unused]] inline std::shared_ptr<Minimizer> create_minimizer(type t, const std::function<double(std::vector<double>)>& func, const Parameter& param, unsigned int evals) {
        auto minimizer = create_minimizer(t, func, param);
        minimizer->set_max_evals(evals);
        return minimizer;
    }
}