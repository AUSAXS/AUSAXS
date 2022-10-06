#include <vector>
#include <string>
#include <fstream>

#include <curl/curl.h>

namespace curl {
    /**
     * @brief Download the given URL as a file. 
     * 
     * Based on the example from https://curl.se/libcurl/c/url2file.html
     */ 
    static void download(std::string url, std::string path) {
        CURL *curl;
        CURLcode res = CURLE_FAILED_INIT;
        FILE *fp;
        curl = curl_easy_init();
        if (curl) {
            fp = fopen(path.c_str(), "wb");
            curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
            res = curl_easy_perform(curl);
            curl_easy_cleanup(curl);
            curl_global_cleanup();
            fclose(fp);
        }

        if (res == CURLE_OK) {
            utility::print_success("\tSuccessfully downloaded " + url + " to " + path);
        } else {
            utility::print_warning("\tFailed to download " + url);
            exit(1);
        }
    }
}