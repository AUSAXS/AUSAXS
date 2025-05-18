#include <utility/Logging.h>

#include <io/File.h>
#include <constants/Version.h>
#include <settings/GeneralSettings.h>

#include <iostream>
#include <fstream>
#include <chrono>
#include <ctime>
#include <memory>

using namespace ausaxs;

struct ILogger {
    virtual void start(const std::string& name) = 0;
    virtual void log(std::string_view msg) = 0;
    virtual void log_console(std::string_view msg) = 0;
    virtual void log_critical(std::string_view msg) = 0;
};

template<bool to_console>
struct Logger : ILogger {
    std::ofstream log_file;
    void start(const std::string& name) override {
        if constexpr (to_console) {
            return;
        }

        std::cout << "Logging to " << settings::general::cache + name + "_log.txt" << std::endl;
        auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::string current_time = std::ctime(&t);
        log_file = std::ofstream(settings::general::cache + name + "_log.txt");
        log_file << "AUSAXS " << constants::version << " log file\nstarted at " << current_time << "\n";
    }

    void log(std::string_view msg) override {
        if constexpr (to_console) {
            std::cout << "LOG: " << msg << std::endl;
            return;
        }

        if (!log_file.is_open()) {return;}
        log_file << msg << "\n";
    }

    void log_console(std::string_view msg) override {
        if (!log_file.is_open()) {return;}
        log_file << "CONSOLE: " << msg << std::endl;
    }

    void log_critical(std::string_view msg) override {
        log(msg);
        log_file.flush();
    }
}; 
std::unique_ptr<ILogger> logobj = std::make_unique<Logger<false>>();

void logging::start(std::string_view name) {
    logobj->start(std::string(name));
}

void logging::log_critical(std::string_view msg) {
    logobj->log_console(msg);
    logobj->log("CRITICAL: " + std::string(msg));
}

void logging::log(std::string_view msg) {
    logobj->log(msg);
}

void logging::log_console(std::string_view msg) {
    logobj->log_console(msg);
}

void logging::log_to_console(bool enable) {
    if (enable) {
        logobj = std::make_unique<Logger<true>>();
    } else {
        logobj = std::make_unique<Logger<false>>();
    }
}