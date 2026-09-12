// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <fitter/FitResult.h>
#include <io/IOFwd.h>
#include <utility/observer_ptr.h>

#include <functional>
#include <string>

namespace ausaxs::fitter {
    class FitReporter {
        public:
            /**
             * @brief Print the fit result to the console.
             */
            static void report(observer_ptr<FitResult> fit);

            /**
             * @brief Print multiple fit results to the console, with optional titles.
             */
            static void report(const std::vector<FitResult>& fits, const std::vector<std::string>& titles = {});

            /**
             * @brief Save the fit result to a file, storing the command line arguments.
             */
            static void save(observer_ptr<FitResult> fit, const io::File& path, int argc, char const* const* argv);

            /**
             * @brief Save the fit result to a file, with an optional header.
             */
            static void save(observer_ptr<FitResult> fit, const io::File& path, const std::string& header = "");

            /**
             * @brief Save multiple fit results to a file, with optional titles.
             */
            static void save(const std::vector<FitResult>& fits, const io::File& path, const std::vector<std::string>& titles = {});

        private:
            [[nodiscard]] static std::function<std::string(std::string)> get_title_reporter(const std::vector<std::string>& titles); 
    };
}