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
        CURLcode res;
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
    }
}

// static size_t write_data(void *ptr, size_t size, size_t nmemb, void *stream) {
//     size_t written = fwrite(ptr, size, nmemb, (FILE *)stream);
//     return written;
// }

// int main(int, char const**) {
//     std::string file = "ASH";

//     CURL *curl_handle;
//     FILE *pagefile;
//     std::string url = "https://files.rcsb.org/ligands/view/" + file + ".cif";
//     std::string out = file + ".cif";

//     curl_global_init(CURL_GLOBAL_ALL);
    
//     /* init the curl session */
//     curl_handle = curl_easy_init();
    
//     /* set URL to get here */
//     curl_easy_setopt(curl_handle, CURLOPT_URL, url.c_str());
    
//     /* Switch on full protocol/debug output while testing */
//     curl_easy_setopt(curl_handle, CURLOPT_VERBOSE, 1L);
    
//     /* disable progress meter, set to 0L to enable it */
//     curl_easy_setopt(curl_handle, CURLOPT_NOPROGRESS, 0L);
    
//     /* send all data to this function  */
//     curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);
    
//     /* open the file */
//     pagefile = fopen(out.c_str(), "wb");
//     if (pagefile) {
//         /* write the page body to this file handle */
//         curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, pagefile);
    
//         /* get it! */
//         curl_easy_perform(curl_handle);
    
//         /* close the header file */
//         fclose(pagefile);
//     }

//     curl_easy_cleanup(curl_handle);
//     curl_global_cleanup();
//     return 0;
// }