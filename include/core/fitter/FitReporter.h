#pragma once

#include <io/IOFwd.h>
#include <fitter/FitResult.h>
#include <utility/observer_ptr.h>

#include <string>
#include <functional>

namespace fitter {
    class FitReporter {
        public:
            static void report(const observer_ptr<FitResult> fit);
            static void report(const std::vector<FitResult>& fits, const std::vector<std::string>& titles = {});

            static void save(const observer_ptr<FitResult> fit, const io::File& path, const std::string& header = "");
            static void save(const std::vector<FitResult>& fits, const io::File& path, const std::vector<std::string>& titles = {});

        private:
            [[nodiscard]] static std::function<std::string(std::string)> get_title_reporter(const std::vector<std::string>& titles); 
    };
}