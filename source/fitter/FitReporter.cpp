#include <fitter/FitReporter.h>
#include <utility/Exceptions.h>
#include <utility/Utility.h>

#include <string>
#include <fstream>
#include <iostream>

template<FitType T>
void FitReporter::report(const T& fit) {
    std::cout << fit.to_string() << std::endl;
}

template<FitType T>
void FitReporter::report(const std::shared_ptr<T> fit) {report(*fit);}

template<FitType T>
void FitReporter::report(const std::vector<T>& fits, std::vector<std::string> titles) {
    if (!titles.empty() && titles.size() != fits.size()) {throw except::size_error("FitReporter::report: Size of fits and titles must be equal.");}

    auto title_reporter = get_title_reporter(titles);
    for (unsigned int i = 0; i < fits.size(); i++) {
        std::cout << title_reporter(titles[i]);
        std::cout << fits[i].to_string() << std::endl;
    }
}

template<FitType T>
void FitReporter::save(const T& fit, std::string path) {
    utility::create_directory(path);

    std::ofstream out(path);
    if (!out.is_open()) {throw except::io_error("FitReporter::save: Could not open file path \"" + path + "\".");}
    out << fit.to_string() << std::endl;
    out.close();
}

template<FitType T>
void FitReporter::save(const std::shared_ptr<T> fit, std::string path) {
    save(*fit, path);
}

template<FitType T>
void FitReporter::save(const std::vector<T>& fits, std::string path, std::vector<std::string> titles) {
    if (!titles.empty() && titles.size() != fits.size()) {throw except::size_error("FitReporter::report: Size of fits and titles must be equal.");}
    utility::create_directory(path);

    std::ofstream out(path);
    if (!out.is_open()) {throw except::io_error("FitReporter::save: Could not open file path \"" + path + "\".");}

    auto title_reporter = get_title_reporter(titles);
    for (unsigned int i = 0; i < fits.size(); i++) {
        out << title_reporter(titles[i]);
        out << fits[i].to_string() << std::endl;
    }
    out.close();
}

std::function<std::string(std::string)> FitReporter::get_title_reporter(std::vector<std::string> titles) {
    std::function<std::string(std::string)> title_reporter;
    if (titles.empty()) {
        title_reporter = [] (std::string) {return "";};
    } else {
        title_reporter = [] (std::string title) {
            string output;
            output += "\n+----------------------------------------------------------+";

            int spaces = (60 - title.size())/2 - 1;
            string spacing(spaces, ' ');
            output += "\n|" + spacing + title + (title.size() % 2 == 0 ? spacing : spacing + " ") + "|\n";
            return output;
        };
    }

    return title_reporter;
}

template void FitReporter::report(const Fit& fit);
template void FitReporter::report(const std::shared_ptr<Fit> fit);
template void FitReporter::report(const std::vector<Fit>& fits, std::vector<std::string> titles);
template void FitReporter::save(const Fit& fit, std::string path);
template void FitReporter::save(const std::shared_ptr<Fit> fit, std::string path);
template void FitReporter::save(const std::vector<Fit>& fits, std::string path, std::vector<std::string> titles);

template void FitReporter::report(const EMFit& fit);
template void FitReporter::report(const std::shared_ptr<EMFit> fit);
template void FitReporter::report(const std::vector<EMFit>& fits, std::vector<std::string> titles);
template void FitReporter::save(const EMFit& fit, std::string path);
template void FitReporter::save(const std::shared_ptr<EMFit> fit, std::string path);
template void FitReporter::save(const std::vector<EMFit>& fits, std::string path, std::vector<std::string> titles);