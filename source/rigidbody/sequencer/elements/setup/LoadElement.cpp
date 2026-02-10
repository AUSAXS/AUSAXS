// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/sequencer/elements/setup/LoadElement.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>
#include <rigidbody/constraints/ConstraintManager.h>
#include <rigidbody/Rigidbody.h>
#include <rigidbody/BodySplitter.h>
#include <data/Molecule.h>
#include <settings/GeneralSettings.h>

#include <algorithm>

using namespace ausaxs::rigidbody::sequencer;

LoadElement::~LoadElement() = default;

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& paths, const std::vector<std::string>& body_names) : owner(owner) {
    if (auto loc = paths[0].find("%"); loc != std::string::npos) {
        rigidbody = std::make_unique<Rigidbody>(data::Molecule(load_wildcarded(paths[0])));
    } else {
        auto rel_paths = paths;
        std::transform(paths.begin(), paths.end(), rel_paths.begin(), [this] (const std::string& path) {return lookup_file(path).first;});
        rigidbody = std::make_unique<Rigidbody>(data::Molecule(rel_paths));
    }

    // add default names
    for (unsigned int i = 0; i < rigidbody->molecule.size_body(); ++i) {
        owner->setup()._get_body_names().emplace("b" + std::to_string(i+1), detail::to_index(i));
    }

    // add custom names
    if (!body_names.empty()) {
        if (body_names.size() != rigidbody->molecule.size_body()) {
            throw std::runtime_error("LoadElement::LoadElement: The number of body names does not match the number of bodies.");
        }
        for (unsigned int i = 0; i < rigidbody->molecule.size_body(); ++i) {
            owner->setup()._get_body_names().emplace(body_names[i], detail::to_index(i));
        }
    }
    owner->setup()._set_active_body(rigidbody.get());

    if (settings::general::verbose) {
        std::cout << "\tLoaded " << rigidbody->molecule.size_body() << " bodies from " << paths.size() << " files." << std::endl;
    }
}

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<int>& splits, const std::vector<std::string>& body_names) 
    : owner(owner) 
{
    if (auto loc = path.find("%"); loc != std::string::npos) {
        rigidbody = std::make_unique<Rigidbody>(data::Molecule(load_wildcarded(path)));
    } else {
        rigidbody = std::make_unique<Rigidbody>(rigidbody::BodySplitter::split(lookup_file(path).first, splits));
    }

    // add default names
    for (unsigned int i = 0; i < rigidbody->molecule.size_body(); ++i) {
        owner->setup()._get_body_names().emplace("b" + std::to_string(i+1), detail::to_index(i));
    }

    // add custom names
    if (!body_names.empty()) {
        if (body_names.size() != rigidbody->molecule.size_body()) {
            throw std::runtime_error("LoadElement::LoadElement: The number of body names does not match the number of bodies.");
        }
        for (unsigned int i = 0; i < rigidbody->molecule.size_body(); ++i) {
            owner->setup()._get_body_names().emplace(body_names[i], detail::to_index(i));
        }
    }
    owner->setup()._set_active_body(rigidbody.get());

    if (settings::general::verbose) {
        std::cout << "\tLoaded " << rigidbody->molecule.size_body() << " bodies from \"" << path << "\"." << std::endl;
    }
}

LoadElement::LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<std::string>& body_names) : owner(owner) {
    rigidbody = std::make_unique<Rigidbody>(rigidbody::BodySplitter::split(lookup_file(path).first));

    // add default names
    for (unsigned int i = 0; i < rigidbody->molecule.size_body(); ++i) {
        owner->setup()._get_body_names().emplace("b" + std::to_string(i+1), detail::to_index(i));
    }

    // add custom names
    if (!body_names.empty()) {
        if (body_names.size() != rigidbody->molecule.size_body()) {
            throw std::runtime_error("LoadElement::LoadElement: The number of body names does not match the number of bodies.");
        }
        for (unsigned int i = 0; i < rigidbody->molecule.size_body(); ++i) {
            owner->setup()._get_body_names().emplace(body_names[i], detail::to_index(i));
        }
    }
    owner->setup()._set_active_body(rigidbody.get());

    if (settings::general::verbose) {
        std::cout << "\tLoaded " << rigidbody->molecule.size_body() << " bodies from \"" << path << "\"." << std::endl;
    }
}

/**
 * @brief Looks up a file relative to both the current working directory and the configuration file directory.
 * @return A pair containing the resolved file path and a boolean indicating whether the file exists.
 */
std::pair<std::string, bool> LoadElement::lookup_file(const std::string& path) {
    io::File file(path);
    if (file.exists()) {return {file, true};}

    auto config_folder = owner->setup()._get_config_folder();
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
    owner->_get_sequencer()->_set_rigidbody(rigidbody.get());
}