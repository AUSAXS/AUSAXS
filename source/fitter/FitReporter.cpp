#include <fitter/FitReporter.h>
#include <utility/Exceptions.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>

#include <fstream>
#include <iostream>

using namespace fitter;

template<FitType T>
void FitReporter::report(const T& fit) {
    std::cout << fit.to_string() << std::endl;
}

template<FitType T>
void FitReporter::report(const std::shared_ptr<T> fit) {report(*fit);}

template<FitType T>
void FitReporter::report(const std::vector<T>& fits, const std::vector<std::string>& titles) {
    if (!titles.empty() && titles.size() != fits.size()) {throw except::size_error("FitReporter::report: Size of fits and titles must be equal.");}

    auto title_reporter = get_title_reporter(titles);
    for (unsigned int i = 0; i < fits.size(); i++) {
        std::cout << title_reporter(titles[i]);
        std::cout << fits[i].to_string() << std::endl;
    }
}

template<FitType T>
void FitReporter::save(const T& fit, const io::File& path) {
    path.directory().create();

    std::ofstream out(path);
    if (!out.is_open()) {throw except::io_error("FitReporter::save: Could not open file path \"" + path + "\".");}
    out << fit.to_string() << std::endl;
    out.close();
}

template<FitType T>
void FitReporter::save(const std::shared_ptr<T> fit, const io::File& path, const std::string& header) {
    path.directory().create();

    std::ofstream out(path);
    if (!out.is_open()) {throw except::io_error("FitReporter::save: Could not open file path \"" + path + "\".");}
    out << header << std::endl;
    out << fit->to_string() << std::endl;
    out.close();
}

template<FitType T>
void FitReporter::save(const std::shared_ptr<T> fit, const io::File& path) {
    save(*fit, path);
}

template<FitType T>
void FitReporter::save(const std::vector<T>& fits, const io::File& path, const std::vector<std::string>& titles) {
    if (!titles.empty() && titles.size() != fits.size()) {throw except::size_error("FitReporter::report: Size of fits and titles must be equal.");}
    path.directory().create();

    std::ofstream out(path);
    if (!out.is_open()) {throw except::io_error("FitReporter::save: Could not open file path \"" + path + "\".");}

    auto title_reporter = get_title_reporter(titles);
    for (unsigned int i = 0; i < fits.size(); i++) {
        out << title_reporter(titles[i]);
        out << fits[i].to_string() << std::endl;
    }
    out.close();
}

std::function<std::string(std::string)> FitReporter::get_title_reporter(const std::vector<std::string>& titles) {
    std::function<std::string(const std::string&)> title_reporter;
    if (titles.empty()) {
        title_reporter = [] (const std::string&) {return "";};
    } else {
        title_reporter = [] (const std::string& title) {
            std::string output;
            output += "\n+----------------------------------------------------------+";

            int spaces = (60 - title.size())/2 - 1;
            std::string spacing(spaces, ' ');
            output += "\n|" + spacing + title + (title.size() % 2 == 0 ? spacing : spacing + " ") + "|\n";
            return output;
        };
    }

    return title_reporter;
}

template void FitReporter::report(const Fit&);
template void FitReporter::report(const std::shared_ptr<Fit>);
template void FitReporter::report(const std::vector<Fit>&, const std::vector<std::string>&);
template void FitReporter::save(const Fit&, const io::File&);
template void FitReporter::save(const std::shared_ptr<Fit>, const io::File&);
template void FitReporter::save(const std::shared_ptr<Fit>, const io::File&, const std::string&);
template void FitReporter::save(const std::vector<Fit>&, const io::File&, const std::vector<std::string>&);

template void FitReporter::report(const EMFit&);
template void FitReporter::report(const std::shared_ptr<EMFit>);
template void FitReporter::report(const std::vector<EMFit>&, const std::vector<std::string>&);
template void FitReporter::save(const EMFit&, const io::File&);
template void FitReporter::save(const std::shared_ptr<EMFit>, const io::File&);
template void FitReporter::save(const std::shared_ptr<EMFit>, const io::File&, const std::string&);
template void FitReporter::save(const std::vector<EMFit>&, const io::File&, const std::vector<std::string>&);