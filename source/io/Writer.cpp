#include <io/Writer.h>
#include <utility/Exceptions.h>

void Writer::write(std::string) {
    throw except::unexpected("Error in Writer::write: This code should be unreachable.");
}