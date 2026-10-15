// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#ifdef _MSC_VER
    #pragma warning(disable:4996) // disable fopen deprecation warning on MSVC
#endif

#include <io/File.h>
#include <settings/GeneralSettings.h>
#include <utility/Console.h>
#include <utility/Curl.h>
#include <utility/Exceptions.h>

#include <curl/curl.h>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <memory>

using namespace ausaxs;

bool curl::download(const std::string& url, const io::File& path) {
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
    // download to a temporary file and only move it into place once the transfer succeeded, so that a failed or
    // truncated transfer can never leave a file behind at the destination (or overwrite a good one)
    io::File tmp(path.path() + ".part");

    {
        FILE* raw_fp = fopen(tmp.path().c_str(), "wb");
        if (raw_fp == nullptr) {throw ausaxs::except::runtime_error("curl::download: Failed to open destination file: \"" + tmp.path() + "\"");}
        std::unique_ptr<FILE, int(*)(FILE*)> fp(raw_fp, &fclose);

        CURLcode res = curl_easy_setopt(raw_curl, CURLOPT_URL, url.c_str());
        if (res != CURLE_OK) {fp.reset(); tmp.remove(); throw ausaxs::except::runtime_error("curl::download: Failed to set URL: \"" + url + "\".");}

        res = curl_easy_setopt(raw_curl, CURLOPT_WRITEDATA, fp.get());
        if (res != CURLE_OK) {fp.reset(); tmp.remove(); throw ausaxs::except::runtime_error("curl::download: Failed to set write data: \"" + path.path() + "\".");}

        // treat HTTP errors (e.g. a 404 page) as failures instead of saving the error page as the file
        res = curl_easy_setopt(raw_curl, CURLOPT_FAILONERROR, 1L);
        if (res != CURLE_OK) {fp.reset(); tmp.remove(); throw ausaxs::except::runtime_error("curl::download: Failed to set fail-on-error.");}

        res = curl_easy_perform(raw_curl);
        bool closed = fclose(fp.release()) == 0; // a failed close means the data may not be on disk
        if (res != CURLE_OK || !closed) {
            tmp.remove();
            console::print_warning(std::string("curl::download: Failed to download \"") + url + "\": " + (res != CURLE_OK ? curl_easy_strerror(res) : "could not write the file"));
            return false;
        }
    }

    std::filesystem::rename(tmp.path(), path.path());
    if (settings::general::verbose) {console::print_success("Successfully downloaded " + url + " to " + path.str());}
    return true;
}
