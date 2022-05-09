#include <utility/Multiset.h>
#include <utility/Exceptions.h>

Multiset::Multiset(const std::vector<Dataset>& data) : data(data) {}

Multiset::Multiset(const Dataset& data) : data({data}) {}

Multiset::Multiset(const Dataset& data1, const Dataset& data2) : data({data1, data2}) {}

const Dataset& Multiset::operator[](unsigned int i) const {
    return data[i];
}

Dataset& Multiset::operator[](unsigned int i) {
    return data[i];
}

const Dataset& Multiset::get_data(std::string name) const {
    if (names.count(name) == 0) {throw except::unknown_argument("Error in Multiset::get_data: No dataset named \"" + name + "\".");}
    return data[names.at(name)];
}

Dataset& Multiset::get_data(std::string name) {
    if (names.count(name) == 0) {throw except::unknown_argument("Error in Multiset::get_data: No dataset named \"" + name + "\".");}
    return data[names.at(name)];
}

const Dataset& Multiset::get_data(unsigned int i) const {
    return data[i];
}

Dataset& Multiset::get_data(unsigned int i) {
    return data[i];
}

size_t Multiset::size() const {return data.size();}

void Multiset::push_back(const Dataset& new_data) {
    data.push_back(new_data);
}

void Multiset::push_back(const Dataset&& new_data) {
    data.push_back(std::move(new_data));
}