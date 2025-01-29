#include <io/IOFwd.h>

#include <string>

namespace ausaxs::curl {
    /**
     * @brief Download the given URL as a file. 
     * @return True if the download was successful, false otherwise.
     */ 
    bool download(const std::string& url, const io::File& path);
}