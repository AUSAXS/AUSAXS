#pragma once

#include <fitter/Fit.h>
#include <concepts>

class FitReporter {
    public:
        template<FitType T>
        static void report(const T& fit);

        template<FitType T>
        static void report(const std::shared_ptr<T> fit);

        template<FitType T>
        static void report(const std::vector<T>& fits, std::vector<std::string> titles = {});

        template<FitType T>
        static void save(std::string path, const T& fit);

        template<FitType T>
        static void save(std::string path, const std::shared_ptr<T> fit);

        template<FitType T>
        static void save(std::string path, const std::vector<T>& fits, std::vector<std::string> titles = {});

  private:
        static std::function<std::string(std::string)> get_title_reporter(std::vector<std::string> titles); 
};