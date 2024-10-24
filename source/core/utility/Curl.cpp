/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#ifdef _MSC_VER
    #pragma warning(disable:4996) // disable fopen deprecation warning on MSVC
#endif

#include <utility/Curl.h>
#include <utility/Console.h>
#include <io/File.h>
#include <settings/GeneralSettings.h>

#include <curl/curl.h>

#include <stdexcept>

using namespace ausaxs;

// based on the example from https://curl.se/libcurl/c/url2file.html
void curl::download(const std::string& url, const io::File& path) {
    CURL *curl;
    CURLcode res = CURLE_FAILED_INIT;
    curl = curl_easy_init();
    if (curl) {
        FILE* fp = fopen(path.path().c_str(), "wb");
        res = curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        if (res != CURLE_OK) {throw std::runtime_error("curl::download: Failed to set URL: " + url);}

        res = curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
        if (res != CURLE_OK) {throw std::runtime_error("curl::download: Failed to set write data: " + path.path());}

        res = curl_easy_perform(curl);
        curl_easy_cleanup(curl);
        curl_global_cleanup();
        fclose(fp);
    }

    if (res == CURLE_OK) {
        if (settings::general::verbose) {console::print_success("\tSuccessfully downloaded " + url + " to " + path);}
        return;
    }
    console::print_warning("\tcurl::download: Failed to download " + url);
    throw std::runtime_error("Failed to download " + url);
}