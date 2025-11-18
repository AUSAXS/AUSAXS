// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#ifdef _MSC_VER
    #pragma warning(disable:4996) // disable fopen deprecation warning on MSVC
#endif

#include <utility/Curl.h>
#include <utility/Console.h>
#include <io/File.h>
#include <settings/GeneralSettings.h>

#include <curl/curl.h>

#include <stdexcept>
#include <cstdlib>
#include <memory>
#include <cstdio>

using namespace ausaxs;

bool curl::download(const std::string& url, const io::File& path) {
    if (settings::general::offline) {
        console::print_warning("curl::download: Offline mode is enabled. Skipping download of \"" + url + "\".");
        return false;
    }

    static bool curl_inited = [](){
        CURLcode cres = curl_global_init(CURL_GLOBAL_DEFAULT);
        if (cres != CURLE_OK) {throw std::runtime_error(std::string("curl::download: curl_global_init failed: ") + curl_easy_strerror(cres));}
        std::atexit([](){ curl_global_cleanup(); });
        return true;
    }();
    (void)curl_inited; // suppress unused variable warning

    CURL* raw_curl = curl_easy_init();
    if (!raw_curl) {
        console::print_warning("curl::download: Failed to create CURL handle.");
        return false;
    }
    std::unique_ptr<CURL, decltype(&curl_easy_cleanup)> curl_ptr(raw_curl, &curl_easy_cleanup);
    FILE* raw_fp = fopen(path.path().c_str(), "wb");
    if (!raw_fp) {throw std::runtime_error("curl::download: Failed to open destination file: \"" + path.path() + "\"");}
    std::unique_ptr<FILE, int(*)(FILE*)> fp(raw_fp, &fclose);

    CURLcode res = curl_easy_setopt(raw_curl, CURLOPT_URL, url.c_str());
    if (res != CURLE_OK) {throw std::runtime_error("curl::download: Failed to set URL: \"" + url + "\".");}

    res = curl_easy_setopt(raw_curl, CURLOPT_WRITEDATA, fp.get());
    if (res != CURLE_OK) {throw std::runtime_error("curl::download: Failed to set write data: \"" + path.path() + "\".");}

    res = curl_easy_perform(raw_curl);
    if (res == CURLE_OK) {
        if (settings::general::verbose) {console::print_success("Successfully downloaded " + url + " to " + path.str());}
        return true;
    }

    console::print_warning(std::string("curl::download: Failed to download \"") + url + "\": " + curl_easy_strerror(res));
    return false;
}