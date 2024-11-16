/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <fitter/FitReporter.h>
#include <utility/Exceptions.h>
#include <mini/detail/FittedParameter.h>
#include <mini/detail/Evaluation.h>

#include <fstream>
#include <iostream>

using namespace ausaxs;
using namespace ausaxs::fitter;

void FitReporter::report(const observer_ptr<FitResult> fit) {
    std::cout << fit->to_string() << std::endl;
}

void FitReporter::report(const std::vector<FitResult>& fits, const std::vector<std::string>& titles) {
    if (!titles.empty() && titles.size() != fits.size()) {throw except::size_error("FitReporter::report: Size of fits and titles must be equal.");}

    auto title_reporter = get_title_reporter(titles);
    for (unsigned int i = 0; i < fits.size(); i++) {
        std::cout << title_reporter(titles[i]);
        std::cout << fits[i].to_string() << std::endl;
    }
}

void FitReporter::save(const observer_ptr<FitResult> fit, const io::File& path, int argc, char const* argv[]) {
    std::string cmd_line;
    for (int i = 0; i < argc; ++i) {cmd_line.append(argv[i]).append(" ");}
    save(fit, path, cmd_line);
}

void FitReporter::save(const observer_ptr<FitResult> fit, const io::File& path, const std::string& header) {
    path.directory().create();

    std::ofstream out(path);
    if (!out.is_open()) {throw except::io_error("FitReporter::save: Could not open file path \"" + path.str() + "\".");}
    if (!header.empty()) {out << header << std::endl;}
    out << fit->to_string() << std::endl;
    out.close();
}

void FitReporter::save(const std::vector<FitResult>& fits, const io::File& path, const std::vector<std::string>& titles) {
    if (!titles.empty() && titles.size() != fits.size()) {throw except::size_error("FitReporter::report: Size of fits and titles must be equal.");}
    path.directory().create();

    std::ofstream out(path);
    if (!out.is_open()) {throw except::io_error("FitReporter::save: Could not open file path \"" + path.str() + "\".");}

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