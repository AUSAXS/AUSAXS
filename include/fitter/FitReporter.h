#pragma once

#include <fitter/Fit.h>

class FitReporter {
    public:
        static void report(const Fit& fit);
        static void report(const std::shared_ptr<Fit> fit);
        static void report(const std::vector<Fit>& fits, std::vector<std::string> titles = {});

        static void save(std::string path, const Fit& fit);
        static void save(std::string path, const std::shared_ptr<Fit> fit);
        static void save(std::string path, const std::vector<Fit>& fits, std::vector<std::string> titles = {});

    private:
        static std::function<std::string(std::string)> get_title_reporter(std::vector<std::string> titles); 
};