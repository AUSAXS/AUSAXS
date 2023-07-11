#include <em/detail/header/HeaderFactory.h>
#include <em/detail/header/MapHeader.h>
#include <em/detail/header/MRCHeader.h>
#include <em/detail/header/RECHeader.h>
#include <io/ExistingFile.h>
#include <utility/Exceptions.h>

std::unique_ptr<em::detail::header::MapHeader> em::detail::factory::create_header(const io::ExistingFile& path) {
    if (em::detail::header::MRCHeader::is_mrc(path)) {
        return std::make_unique<em::detail::header::MRCHeader>();
    } else if (em::detail::header::RECHeader::is_rec(path)) {
        return std::make_unique<em::detail::header::RECHeader>();
    } else {
        throw except::io_error("Unsupported map data format: \"" + path.extension() + "\"");
    }
}