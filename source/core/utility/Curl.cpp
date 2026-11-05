// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#ifdef _MSC_VER
    #pragma warning(disable:4996) // disable fopen deprecation warning on MSVC
#endif

#include <utility/Curl.h>

#include <constants/Version.h>
#include <io/File.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Exceptions.h>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string_view>

#include <curl/curl.h>

using namespace ausaxs;

namespace {
    struct Progress {
        int last_step = -1;
    };

    // reported in steps of 10%, as whole lines, so the progress reads the same in a terminal and in a log
    int report_progress(void* data, curl_off_t total, curl_off_t now, curl_off_t /* void */, curl_off_t /* void */) {
        if (total <= 0) {return 0;} // size not known yet
        auto* progress = static_cast<Progress*>(data);
        int step = static_cast<int>(10*now/total);
        if (step == progress->last_step) {return 0;}
        progress->last_step = step;
        console::print_text(std::to_string(now/1'000'000) + " / " + std::to_string(total/1'000'000) + " MB (" + std::to_string(10*step) + "%)");
        return 0;
    }
}

bool curl::download(const std::string& url, const io::File& path, bool show_progress) {
    if (settings::general::offline) {
        console::print_warning("curl::download: Offline mode is enabled. Skipping download of \"" + url + "\".");
        return false;
    }

    static bool curl_inited = [](){
        CURLcode cres = curl_global_init(CURL_GLOBAL_DEFAULT);
        if (cres != CURLE_OK) {throw except::runtime_error(std::string("curl::download: curl_global_init failed: ") + curl_easy_strerror(cres));}
        std::atexit([](){ curl_global_cleanup(); });
        return true;
    }();
    (void)curl_inited; // suppress unused variable warning

    CURL* raw_curl = curl_easy_init();
    if (raw_curl == nullptr) {
        console::print_warning("curl::download: Failed to create CURL handle.");
        return false;
    }
    std::unique_ptr<CURL, decltype(&curl_easy_cleanup)> curl_ptr(raw_curl, &curl_easy_cleanup);

    // written under a temporary name and only moved into place once the transfer succeeded, so that a failed or
    // truncated transfer can never leave a file behind at the destination (or overwrite a good one)
    io::File partial = path.append(".part");
    {
        FILE* raw_fp = fopen(partial.path().c_str(), "wb");
        if (raw_fp == nullptr) {throw ausaxs::except::runtime_error("curl::download: Failed to open destination file: \"" + partial.path() + "\"");}
        std::unique_ptr<FILE, int(*)(FILE*)> fp(raw_fp, &fclose);

        std::string user_agent = "ausaxs/" + std::string(constants::version);
        Progress progress;
        auto set = [&] (CURLoption option, auto value, std::string_view what) {
            if (curl_easy_setopt(raw_curl, option, value) != CURLE_OK) {
                fp.reset();
                partial.remove();
                throw ausaxs::except::runtime_error("curl::download: Failed to set " + std::string(what) + " for \"" + url + "\".");
            }
        };
        set(CURLOPT_URL, url.c_str(), "URL");
        set(CURLOPT_WRITEDATA, fp.get(), "write data");
        set(CURLOPT_USERAGENT, user_agent.c_str(), "user agent");
        set(CURLOPT_FOLLOWLOCATION, 1L, "redirect following"); // e.g. GitHub release assets redirect to a CDN
        set(CURLOPT_MAXREDIRS, 10L, "redirect limit");
        set(CURLOPT_FAILONERROR, 1L, "error handling");        // otherwise an HTTP error page (e.g. a 404) is saved as if it were the file
        set(CURLOPT_CONNECTTIMEOUT, 30L, "connection timeout");
        set(CURLOPT_LOW_SPEED_LIMIT, 1024L, "stall limit");    // give up on a transfer that has stalled, but never on one that is merely slow
        set(CURLOPT_LOW_SPEED_TIME, 60L, "stall time");
        if (show_progress) {
            set(CURLOPT_XFERINFOFUNCTION, &report_progress, "progress function");
            set(CURLOPT_XFERINFODATA, static_cast<void*>(&progress), "progress data");
            set(CURLOPT_NOPROGRESS, 0L, "progress reporting");
        }

        if (show_progress) {console::indent();}
        CURLcode res = curl_easy_perform(raw_curl);
        if (show_progress) {console::unindent();}
        bool closed = fclose(fp.release()) == 0; // a failed close means the data may not be on disk
        if (res != CURLE_OK || !closed) {
            partial.remove();
            console::print_warning(std::string("curl::download: Failed to download \"") + url + "\": " + (res != CURLE_OK ? curl_easy_strerror(res) : "could not write the file"));
            return false;
        }
    }

    try {
        partial.rename(path.stem());
    } catch (const std::filesystem::filesystem_error& e) {
        console::print_warning("curl::download: Failed to move the download into place at \"" + path.path() + "\": " + e.code().message());
        partial.remove();
        return false;
    }
    if (settings::general::verbose) {console::print_success("Successfully downloaded " + url + " to " + path.str());}
    return true;
}