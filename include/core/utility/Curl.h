#include <curl/curl.h>
#include <io/IOFwd.h>

#include <string>

namespace curl {
    /**
     * @brief Download the given URL as a file. 
     * 
     * Based on the example from https://curl.se/libcurl/c/url2file.html
     */ 
    void download(const std::string& url, const io::File& path);
}