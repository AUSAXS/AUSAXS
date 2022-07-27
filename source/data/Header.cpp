#include <data/Header.h>

void Header::parse_pdb(std::string s) {add(s);}

std::string Header::as_pdb() const {return get();}

void Header::add(const std::string s) {contents.push_back(s);}

std::string Header::get() const {return utility::join(contents, "\n") + (size() == 0 ? "" : "\n");};

#include <iostream>
void Header::remove(std::string type) {
    std::vector<std::string> new_contents;
    new_contents.reserve(contents.size());
    for (unsigned int i = 0; i < contents.size(); i++) {
        if (contents[i].substr(0, 6) != type) {
            new_contents.push_back(contents[i]);
        }
    }
    std::cout << "Old size: " << contents.size() << ", new size: " << new_contents.size() << std::endl;
    contents = std::move(new_contents);
}

Header& Header::operator=(const Header& header) {
    contents = header.contents;
    return *this;
}

size_t Header::size() const {return contents.size();}