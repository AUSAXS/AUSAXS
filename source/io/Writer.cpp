#include <io/Writer.h>
#include <utility/Exceptions.h>

void Writer::write(std::string) {
    throw except::unexpected("Writer::write: This code should be unreachable.");
}
