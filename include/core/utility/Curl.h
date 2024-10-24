#include <io/IOFwd.h>

#include <string>

namespace ausaxs::curl {
    /**
     * @brief Download the given URL as a file. 
     */ 
    void download(const std::string& url, const io::File& path);
}