#include <fitter/FitReporter.h>
#include <Exceptions.h>
#include <string>

void FitReporter::report(const Fit& fit) {
    std::cout << fit.to_string() << std::endl;
}

void FitReporter::report(const std::vector<Fit>& fits, std::vector<std::string> titles) {
    if (!titles.empty() && titles.size() != fits.size()) {throw except::size_error("Error in FitReporter::report: Size of fits and titles must be equal.");}

    auto title_reporter = get_title_reporter(titles);
    for (unsigned int i = 0; i < fits.size(); i++) {
        std::cout << title_reporter(titles[i]);
        std::cout << fits[i].to_string();
    }
}

void FitReporter::save(const Fit& fit, std::string path) {}
void FitReporter::save(const std::vector<Fit>& fits, std::string path) {}

std::function<std::string(std::string)> FitReporter::get_title_reporter(std::vector<std::string> titles) {
    std::function<void(std::string)> title_reporter;
    if (titles.empty()) {
        title_reporter = [] (std::string) {};
    } else {
        title_reporter = [] (std::string title) {
            string output;
            output += "\n+----------------------------------------------------------+";

            int spacing = std::floor((60 - title.size())/2) - 1;
            string space(spacing, ' ');
            output += "\n|" + spacing + title + spacing + "|";
            output += "\n+----------------------------------------------------------+";
        };
    }
}