#pragma once

#include <fitter/Fit.h>

class FitReporter {
    public:
        static void report(const Fit& fit);
        static void report(const std::vector<Fit>& fits, std::vector<std::string> titles = {});

        static void save(const Fit& fit, std::string path);
        static void save(const std::vector<Fit>& fits, std::string path);

    private:
        static std::function<std::string(std::string)> get_title_reporter(std::vector<std::string> titles); 
};