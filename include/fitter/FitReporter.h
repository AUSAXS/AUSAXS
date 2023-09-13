#pragma once

#include <fitter/Fit.h>

#include <string>
#include <functional>

namespace io {class File;}
namespace fitter {
    class FitReporter {
        public:
            template<FitType T>
            static void report(const T& fit);

            template<FitType T>
            static void report(const std::shared_ptr<T> fit);

            template<FitType T>
            static void report(const std::vector<T>& fits, const std::vector<std::string>& titles = {});

            template<FitType T>
            static void save(const T& fit, const io::File& path);

            template<FitType T>
            static void save(const std::shared_ptr<T> fit, const io::File& path);

            template<FitType T>
            static void save(const std::shared_ptr<T> fit, const io::File& path, const std::string& header);

            template<FitType T>
            static void save(const std::vector<T>& fits, const io::File& path, const std::vector<std::string>& titles = {});

    private:
            [[nodiscard]] static std::function<std::string(std::string)> get_title_reporter(const std::vector<std::string>& titles); 
    };
}