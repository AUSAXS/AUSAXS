#pragma once

#include <fitter/Fit.h>
#include <fitter/FitReporter.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

#include <concepts>
#include <fstream>
#include <iostream>

template<typename T>
concept FitType = 
    std::is_base_of_v<Fit, T> && 
    requires (T t) {
        {t.to_string()} -> std::convertible_to<std::string>;
    };

class FitReporter {
    public:
        template<typename FitType>
        static void report(const FitType& fit) {
            std::cout << fit.to_string() << std::endl;
        }

        template<typename FitType>
        static void report(const std::shared_ptr<FitType> fit) {
            report(*fit);
        }

        static void report(const std::vector<Fit>& fits, std::vector<std::string> titles = {});

        template<typename FitType>
        static void save(const FitType& fit, std::string path) {
            utility::create_directory(path);

            std::ofstream out(path);
            if (!out.is_open()) {throw except::io_error("FitReporter::save: Could not open file path \"" + path + "\".");}
            out << fit.to_string() << std::endl;
            out.close();
        }

        template<typename FitType>
        static void save(const std::shared_ptr<FitType> fit, std::string path) {
            save(*fit, path);
        }

        static void save(const std::vector<Fit>& fits, std::string path, std::vector<std::string> titles = {});

    private:
        static std::function<std::string(std::string)> get_title_reporter(std::vector<std::string> titles); 
};