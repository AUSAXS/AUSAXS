#include <curl/curl.h>

#include <string>

namespace io {class File;}
namespace curl {
    /**
     * @brief Download the given URL as a file. 
     * 
     * Based on the example from https://curl.se/libcurl/c/url2file.html
     */ 
    void download(const std::string& url, const io::File& path);
}