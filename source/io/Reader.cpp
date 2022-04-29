#include <io/Reader.h>
#include <utility/Exceptions.h>

void Reader::read(std::string) {
    throw except::unexpected("Error in Reader::read: This code should be unreachable.");
}
