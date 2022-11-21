#pragma once

#include <preprocessor.h>
#include <mini/dlibMinimizer.h>
#include <mini/Golden.h>
#include <mini/MinimumExplorer.h>
#include <mini/Scan.h>
#include <mini/LimitedScan.h>
#include <mini/Utility.h>

#include <memory>
#include <functional>
#include <concepts>

namespace mini {
    enum class type {
        BFGS,
        GOLDEN,
        MINIMUM_EXPLORER,
        SCAN,
        LIMITED_SCAN
    };

    namespace detail {
        static std::shared_ptr<Minimizer> create_minimizer(type t) {
            switch (t) {
                case type::BFGS:
                    return std::make_shared<dlibMinimizer>();
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

    [[maybe_unused]] static std::shared_ptr<Minimizer> create_minimizer(type t, const std::function<double(std::vector<double>)>& func) {
        auto minimizer = detail::create_minimizer(t);
        minimizer->set_function(func);
        return minimizer;
    }

    [[maybe_unused]] static std::shared_ptr<Minimizer> create_minimizer(type t, const std::function<double(std::vector<double>)>& func, const Parameter& param) {
        auto minimizer = create_minimizer(t, func);
        minimizer->add_parameter(param);
        return minimizer;
    }

    [[maybe_unused]] static std::shared_ptr<Minimizer> create_minimizer(type t, const std::function<double(std::vector<double>)>& func, const Parameter& param, unsigned int evals) {
        switch (t) {
            case type::MINIMUM_EXPLORER:
                return std::make_shared<MinimumExplorer>(func, param, evals);
            case type::SCAN:
                return std::make_shared<Scan>(func, param, evals);
            case type::LIMITED_SCAN:
                return std::make_shared<LimitedScan>(func, param, evals);
            default:
                debug_print("all::create_minimizer: evals is only used for MinimumExplorer, Scan, and LimitedScan.");
                return create_minimizer(t, func, param);
        }
    }
}