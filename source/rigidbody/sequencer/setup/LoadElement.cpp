// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/setup/LoadElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/RigidBody.h>
#include <rigidbody/BodySplitter.h>
#include <settings/GeneralSettings.h>

#include <iostream>
#include <algorithm>

using namespace ausaxs::rigidbody::sequencer;

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& paths, const std::vector<std::string>& body_names, const std::string& saxs_path) : owner(owner) {
    if (auto loc = paths[0].find("%"); loc != std::string::npos) {
        rigidbody = std::make_unique<RigidBody>(data::Molecule(load_wildcarded(paths[0])));
    } else {
        auto rel_paths = paths;
        std::transform(paths.begin(), paths.end(), rel_paths.begin(), [this] (const std::string& path) {return lookup_file(path).first;});
        rigidbody = std::make_unique<RigidBody>(data::Molecule(paths));
    }

    // add default names
    for (unsigned int i = 0; i < rigidbody->size_body(); ++i) {
        owner->_get_body_names().emplace("b" + std::to_string(i+1), i);
    }

    // add custom names
    if (!body_names.empty() && body_names.size() != rigidbody->size_body()) {throw std::runtime_error("LoadElement::LoadElement: The number of body names does not match the number of bodies.");}
    if (!body_names.empty()) {
        for (unsigned int i = 0; i < rigidbody->size_body(); ++i) {
            owner->_get_body_names().emplace(body_names[i], i);
        }
    }
    owner->_set_active_body(rigidbody.get());

    if (!saxs_path.empty()) {
        auto path = lookup_file(saxs_path);
        owner->_set_saxs_path(path.first);
    }

    if (settings::general::verbose) {
        std::cout << "\tLoaded " << rigidbody->size_body() << " bodies from " << paths.size() << " files." << std::endl;
    }
}

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<int>& splits, const std::vector<std::string>& body_names, const std::string& saxs_path) 
    : owner(owner) 
{
    if (auto loc = path.find("%"); loc != std::string::npos) {
        rigidbody = std::make_unique<RigidBody>(data::Molecule(load_wildcarded(path)));
    } else {
        rigidbody = std::make_unique<RigidBody>(rigidbody::BodySplitter::split(lookup_file(path).first, splits));
    }

    // add default names
    for (unsigned int i = 0; i < rigidbody->size_body(); ++i) {
        owner->_get_body_names().emplace("b" + std::to_string(i+1), i);
    }

    // add custom names
    if (!body_names.empty() && body_names.size() != rigidbody->size_body()) {throw std::runtime_error("LoadElement::LoadElement: The number of body names does not match the number of bodies.");}
    if (!body_names.empty()) {
        for (unsigned int i = 0; i < rigidbody->size_body(); ++i) {
            owner->_get_body_names().emplace(body_names[i], i);
        }
    }
    owner->_set_active_body(rigidbody.get());

    if (!saxs_path.empty()) {
        auto path = lookup_file(saxs_path);
        owner->_set_saxs_path(path.first);
    }

    if (settings::general::verbose) {
        std::cout << "\tLoaded " << rigidbody->size_body() << " bodies from \"" << path << "\"." << std::endl;
    }
}

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<std::string>& body_names, const std::string& saxs_path) : owner(owner) {
    rigidbody = std::make_unique<RigidBody>(rigidbody::BodySplitter::split(lookup_file(path).first));

    // add default names
    for (unsigned int i = 0; i < rigidbody->size_body(); ++i) {
        owner->_get_body_names().emplace("b" + std::to_string(i+1), i);
    }

    // add custom names
    if (!body_names.empty() && body_names.size() != rigidbody->size_body()) {throw std::runtime_error("LoadElement::LoadElement: The number of body names does not match the number of bodies.");}
    if (!body_names.empty()) {
        for (unsigned int i = 0; i < rigidbody->size_body(); ++i) {
            owner->_get_body_names().emplace(body_names[i], i);
        }
    }
    owner->_set_active_body(rigidbody.get());

    if (!saxs_path.empty()) {
        auto path = lookup_file(saxs_path);
        owner->_set_saxs_path(path.first);
    }

    if (settings::general::verbose) {
        std::cout << "\tLoaded " << rigidbody->size_body() << " bodies from \"" << path << "\"." << std::endl;
    }
}

std::pair<std::string, bool> LoadElement::lookup_file(const std::string& path) {
    io::File file(path);
    if (file.exists()) {return {file, true};}

    auto config_folder = owner->_get_config_folder();
    io::File relative(config_folder, file.stem(), file.extension());
    if (relative.exists()) {return {relative, true};}

    return {file, false};
}

std::vector<std::string> LoadElement::load_wildcarded(const std::string& path) {
    static auto zero_pad_string = [] (int val, unsigned int pad) -> std::string {
        std::string s = std::to_string(val);
        if (s.size() < pad) {
            s.insert(0, pad - s.size(), '0');
        }
        return s;
    };
    std::vector<std::string> wildcarded_files;

    auto loc = path.find("%");
    int start = loc, end = loc;
    while (path[++end] == '%') {if (100 < end - start) {throw std::runtime_error("LoadElement::LoadElement: The maximum number of consecutive '%' characters is 100.");}}
    int counter = 0;

    // file numbered zero may or may not exist
    io::File file = path.substr(0, start) + zero_pad_string(counter, end - start) + path.substr(end);
    if (auto res = lookup_file(file); res.second) { // check filename padded with zeros
        wildcarded_files.push_back(res.first);
    } else if (1 < end - start) {   // check filename without padding
        file = path.substr(0, start) + std::to_string(counter) + path.substr(end);
        if (res = lookup_file(file); res.second) {wildcarded_files.push_back(res.first);}
    }
    ++counter;

    // check for wildcarded files with padding
    file = path.substr(0, start) + zero_pad_string(counter, end - start) + path.substr(end);
    auto res = lookup_file(file);
    while (res.second) {
        wildcarded_files.push_back(res.first);
        file = path.substr(0, start) + zero_pad_string(++counter, end - start) + path.substr(end);
        res = lookup_file(file);
    }

    // check for wildcarded files without padding
    if (wildcarded_files.size() <= 1) { // only check if no files were found with padding
        file = path.substr(0, start) + std::to_string(counter) + path.substr(end);
        res = lookup_file(file);
        while (res.second) {
            wildcarded_files.push_back(file.path());
            file = path.substr(0, start) + std::to_string(++counter) + path.substr(end);
            res = lookup_file(file);
        }
    }

    if (wildcarded_files.empty()) {throw std::runtime_error("LoadElement::LoadElement: No files found matching the wildcarded path.");}
    return wildcarded_files;
}

void LoadElement::run() {
    owner->_get_rigidbody() = rigidbody.get();
}