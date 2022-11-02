#include <data/Footer.h>

Record::RecordType Footer::get_type() const {return RecordType::FOOTER;}

void Footer::parse_pdb(const std::string s) {add(s);}

std::string Footer::as_pdb() const {return get();}

void Footer::add(const std::string s) {contents.push_back(s);}

void Footer::remove(std::string type) {
    std::vector<std::string> new_contents;
    new_contents.reserve(contents.size());
    for (unsigned int i = 0; i < contents.size(); ++i) {
        if (Record::get_type(contents[i].substr(0, 6)) != Record::get_type(type)) {
            new_contents.push_back(contents[i]);
        }
    }
    contents = std::move(new_contents);
}

std::string Footer::get() const {return utility::join(contents, "\n") + (size() == 0 ? "" : "\n");}

Footer& Footer::operator=(const Footer& footer) {
    contents = footer.contents;
    return *this;
}

size_t Footer::size() const {return contents.size();}