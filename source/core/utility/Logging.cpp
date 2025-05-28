/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <utility/Logging.h>

#include <io/File.h>
#include <constants/Version.h>
#include <settings/GeneralSettings.h>

#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>

using namespace ausaxs;

struct {
    std::ofstream log_file;
    void start(const std::string& name) {
        std::cout << "Logging to " << settings::general::cache + name + "_log.txt" << std::endl;
        auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::string current_time = std::ctime(&t);
        log_file = std::ofstream(settings::general::cache + name + "_log.txt");
        log_file << "AUSAXS " << constants::version << " log file\nstarted at " << current_time << "\n";
    }

    void log(std::string_view msg) {
        if (!log_file.is_open()) {return;}
        log_file << msg << "\n";
    }

    void log_console(std::string_view msg) {
        if (!log_file.is_open()) {return;}
        log_file << "CONSOLE: " << msg << std::endl;
    }

    void log_critical(std::string_view msg) {
        log(msg);
        log_file.flush();
    }
} logobj;

void logging::start(std::string_view name) {
    logobj.start(std::string(name));
}

void logging::log_critical(std::string_view msg) {
    logobj.log_console(msg);
    logobj.log("CRITICAL: " + std::string(msg));
}

void logging::log(std::string_view msg) {
    logobj.log(msg);
}

void logging::log_console(std::string_view msg) {
    logobj.log_console(msg);
}