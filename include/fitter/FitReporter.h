#pragma once

#include <fitter/Fit.h>
#include <io/ExistingFile.h>

#include <string>

namespace fitter {
    class FitReporter {
        public:
            template<FitType T>
            static void report(const T& fit);

            template<FitType T>
            static void report(const std::shared_ptr<T> fit);

            template<FitType T>
            static void report(const std::vector<T>& fits, std::vector<std::string> titles = {});

            template<FitType T>
            static void save(const T& fit, const io::File& path);

            template<FitType T>
            static void save(const std::shared_ptr<T> fit, const io::File& path);

            template<FitType T>
            static void save(const std::vector<T>& fits, const io::File& path, std::vector<std::string> titles = {});

    private:
            [[nodiscard]] static std::function<std::string(std::string)> get_title_reporter(std::vector<std::string> titles); 
    };
}