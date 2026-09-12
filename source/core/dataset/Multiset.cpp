// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <dataset/Multiset.h>

#include <utility/Exceptions.h>

#include <algorithm>
#include <filesystem>

using namespace ausaxs;

Multiset::Multiset(const std::vector<Dataset2D>& data) : data(data) {}

Multiset::Multiset(const Dataset2D& data) : data({data}) {}

Multiset::Multiset(const Dataset2D& data1, const Dataset2D& data2) : data({data1, data2}) {}

const Dataset2D& Multiset::operator[](int i) const {
    return data[i];
}

Dataset2D& Multiset::operator[](int i) {
    return data[i];
}

const Dataset2D& Multiset::get_data(const std::string& name) const {
    if (!names.contains(name)) {throw except::unknown_argument("Multiset::get_data: No dataset named \"" + name + "\".");}
    return data[names.at(name)];
}

Dataset2D& Multiset::get_data(const std::string& name) {
    if (!names.contains(name)) {throw except::unknown_argument("Multiset::get_data: No dataset named \"" + name + "\".");}
    return data[names.at(name)];
}

const Dataset2D& Multiset::get_data(int i) const {
    return data[i];
}

Dataset2D& Multiset::get_data(int i) {
    return data[i];
}

std::size_t Multiset::size() const {
    return data.size();
}

bool Multiset::empty() const {
    return data.empty();
}

void Multiset::push_back(const Dataset2D& new_data) {
    data.push_back(new_data);
}

void Multiset::push_back(const Dataset2D&& new_data) {
    data.push_back(new_data);
}

void Multiset::ylimits(double min, double max) {ylimits({min, max});}

void Multiset::ylimits(const Limit& limit) {
    std::ranges::for_each(*this, [&limit] (Dataset2D& data) {data.limit_y(limit);});
}

void Multiset::save(const io::File& path) const {
    for (int i = 0; i < static_cast<int>(size()); i++) {
        data[i].save(path.str() + "/" + std::to_string(i) + ".txt");
    }
}

void Multiset::read(const io::Folder& path) {
    for (const auto& file : std::filesystem::recursive_directory_iterator(path.path())) { // loop over all files in the directory
        if (file.path().extension() == ".txt") {
            Dataset2D set(file.path().string());
            data.push_back(set);
        }
    }
}

std::vector<Dataset2D>::const_iterator Multiset::begin() const {return data.begin();}
std::vector<Dataset2D>::const_iterator Multiset::end() const {return data.end();}

std::vector<Dataset2D>::iterator Multiset::begin() {return data.begin();}
std::vector<Dataset2D>::iterator Multiset::end() {return data.end();}