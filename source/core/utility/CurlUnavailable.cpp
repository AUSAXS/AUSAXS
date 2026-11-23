// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <utility/Curl.h>

#include <utility/Console.h>

bool ausaxs::curl::download(const std::string& url, const io::File& /*unused*/, bool /*unused*/) {
    console::print_warning("curl::download: Downloads are unavailable in the WebAssembly build: \"" + url + "\".");
    return false;
}
