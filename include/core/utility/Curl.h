#include <io/IOFwd.h>

#include <string>

namespace curl {
    /**
     * @brief Download the given URL as a file. 
     */ 
    void download(const std::string& url, const io::File& path);
}