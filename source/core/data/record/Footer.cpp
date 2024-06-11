/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <data/record/Footer.h>
#include <utility/StringUtils.h>

using namespace data::record;

Footer::Footer() noexcept = default;

Footer::~Footer() = default;

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