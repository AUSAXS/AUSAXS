// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/IOFwd.h>

#include <string>

namespace ausaxs::curl {
    /**
     * @brief Download the given URL as a file. 
     *
     * Redirects are followed, and an HTTP error status counts as a failure. The file is written under a temporary name and only moved into 
     * place once complete, so a failed or interrupted download never leaves a partial file at @a path.
     *
     * @param show_progress Print the progress of the download as it runs. Meant for large files.
     * @return True if the download was successful, false otherwise.
     */ 
    bool download(const std::string& url, const io::File& path, bool show_progress = false);
}