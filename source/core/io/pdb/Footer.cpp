// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <io/pdb/Footer.h>
#include <utility/StringUtils.h>

using namespace ausaxs;
using namespace ausaxs::io::pdb;

RecordType Footer::get_type() const {return RecordType::FOOTER;}

void Footer::parse_pdb(const std::string& s) {add(s);}

std::string Footer::as_pdb() const {return get();}

void Footer::add(const std::string& s) {contents.push_back(s);}

void Footer::remove(const std::string& type) {
    auto t = type + std::string(6 - type.size(), ' ');
    std::vector<std::string> new_contents;
    new_contents.reserve(contents.size());
    for (unsigned int i = 0; i < contents.size(); i++) {
        if (contents[i].substr(0, 6) != t) {
            new_contents.push_back(contents[i]);
        }
    }
    contents = std::move(new_contents);
}

std::string Footer::get() const {return utility::join(contents, "\n") + (size() == 0 ? "" : "\n");}

std::size_t Footer::size() const {return contents.size();}

bool Footer::operator==(const Footer& rhs) const = default;