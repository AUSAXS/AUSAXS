// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <mini/dlibMinimizer.h>
#include <mini/Golden.h>
#include <mini/MinimumExplorer.h>
#include <mini/Scan.h>
#include <mini/LimitedScan.h>

#include <memory>
#include <functional>

namespace ausaxs::mini {
    namespace detail {
        inline std::shared_ptr<Minimizer> create_minimizer(algorithm t) {
            switch (t) {
                #if defined(DLIB_AVAILABLE)
                    case algorithm::DLIB_GLOBAL:
                        return std::make_shared<dlibMinimizer<algorithm::DLIB_GLOBAL>>();
                    case algorithm::BFGS:
                        return std::make_shared<dlibMinimizer<algorithm::BFGS>>();
                #endif
                case algorithm::GOLDEN:
                    return std::make_shared<Golden>();
                case algorithm::MINIMUM_EXPLORER:
                    return std::make_shared<MinimumExplorer>();
                case algorithm::SCAN:
                    return std::make_shared<Scan>();
                case algorithm::LIMITED_SCAN:
                    return std::make_shared<LimitedScan>();
                default:
                    throw std::invalid_argument("all::create_minimizer: Unknown minimizer type.");
            }
        }

        [[maybe_unused]] inline std::shared_ptr<Minimizer> create_minimizer(algorithm t, std::function<double(std::vector<double>)>&& func) {
            auto minimizer = create_minimizer(t);
            minimizer->set_function(std::move(func));
            return minimizer;
        }
    }

    /**
     * @brief Create a new minimizer.
     *
     * @param t The algorithm to use.
     * @param func The function to minimize.
     * @param param The first parameter.
     */
    [[maybe_unused]] inline std::shared_ptr<Minimizer> create_minimizer(algorithm t, std::function<double(std::vector<double>)> func, const Parameter& param) {
        auto minimizer = detail::create_minimizer(t, std::move(func));
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
    [[maybe_unused]] inline std::shared_ptr<Minimizer> create_minimizer(algorithm t, std::function<double(std::vector<double>)> func, const std::vector<Parameter>& param) {
        auto minimizer = detail::create_minimizer(t, std::move(func));
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
    [[maybe_unused]] inline std::shared_ptr<Minimizer> create_minimizer(algorithm t, std::function<double(std::vector<double>)> func, const Parameter& param, unsigned int evals) {
        auto minimizer = detail::create_minimizer(t, std::move(func));
        minimizer->add_parameter(param);
        minimizer->set_max_evals(evals);
        return minimizer;
    }
}