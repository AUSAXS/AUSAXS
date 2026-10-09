// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#ifdef _MSC_VER
    #pragma warning(disable:4996) // disable fopen deprecation warning on MSVC
#endif

#include <constants/Version.h>
#include <io/File.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Curl.h>
#include <utility/Exceptions.h>

#include <curl/curl.h>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>

using namespace ausaxs;

namespace {
    struct Progress {
        int last_percent = -1;
    };

    int report_progress(void* data, curl_off_t total, curl_off_t now, curl_off_t, curl_off_t) {
        if (total <= 0) {return 0;} // size not known yet
        auto* progress = static_cast<Progress*>(data);
        int percent = static_cast<int>(100*now/total);
        if (percent == progress->last_percent) {return 0;}
        progress->last_percent = percent;

        char line[64];
        std::snprintf(line, sizeof(line), "\r    %.1f / %.1f MB (%d%%)", static_cast<double>(now)/1e6, static_cast<double>(total)/1e6, percent);
        std::cout << line << std::flush;
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

    // written under a temporary name, so that only a complete download ever appears at the destination
    std::string partial = path.path() + ".part";
    FILE* raw_fp = fopen(partial.c_str(), "wb");
    if (raw_fp == nullptr) {throw ausaxs::except::runtime_error("curl::download: Failed to open destination file: \"" + partial + "\"");}
    std::unique_ptr<FILE, int(*)(FILE*)> fp(raw_fp, &fclose);

    std::string user_agent = "ausaxs/" + std::string(constants::version);
    Progress progress;
    auto set = [&] (CURLoption option, auto value, std::string_view what) {
        if (curl_easy_setopt(raw_curl, option, value) != CURLE_OK) {
            throw ausaxs::except::runtime_error("curl::download: Failed to set " + std::string(what) + " for \"" + url + "\".");
        }
    };
    set(CURLOPT_URL, url.c_str(), "URL");
    set(CURLOPT_WRITEDATA, fp.get(), "write data");
    set(CURLOPT_USERAGENT, user_agent.c_str(), "user agent");
    set(CURLOPT_FOLLOWLOCATION, 1L, "redirect following"); // e.g. GitHub release assets redirect to a CDN
    set(CURLOPT_MAXREDIRS, 10L, "redirect limit");
    set(CURLOPT_FAILONERROR, 1L, "error handling");        // otherwise an HTTP error page is saved as if it were the file
    set(CURLOPT_CONNECTTIMEOUT, 30L, "connection timeout");
    set(CURLOPT_LOW_SPEED_LIMIT, 1024L, "stall limit");    // give up on a transfer that has stalled, but never on one that is merely slow
    set(CURLOPT_LOW_SPEED_TIME, 60L, "stall time");
    if (show_progress) {
        set(CURLOPT_XFERINFOFUNCTION, &report_progress, "progress function");
        set(CURLOPT_XFERINFODATA, static_cast<void*>(&progress), "progress data");
        set(CURLOPT_NOPROGRESS, 0L, "progress reporting");
    }

    CURLcode res = curl_easy_perform(raw_curl);
    if (show_progress && progress.last_percent >= 0) {std::cout << std::endl;}
    bool written = (fflush(fp.get()) == 0);
    fp.reset();

    std::error_code ec;
    if (res == CURLE_OK && written) {
        std::filesystem::rename(partial, path.path(), ec);
        if (!ec) {
            if (settings::general::verbose) {console::print_success("Successfully downloaded " + url + " to " + path.str());}
            return true;
        }
        console::print_warning("curl::download: Failed to move the download into place at \"" + path.path() + "\": " + ec.message());
    } else if (res != CURLE_OK) {
        console::print_warning(std::string("curl::download: Failed to download \"") + url + "\": " + curl_easy_strerror(res));
    } else {
        console::print_warning("curl::download: Failed to write \"" + partial + "\".");
    }
    std::filesystem::remove(partial, ec);
    return false;
}
