// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/BodyNameRegistry.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>

#include <algorithm>
#include <cassert>
#include <map>
#include <stdexcept>

using namespace ausaxs::rigidbody::sequencer::detail;

BodyNameRegistry::BodyLabels BodyNameRegistry::labels(int body) const {
    unsigned int index = static_cast<unsigned int>(to_index(body));
    auto d = defaults.find(index);
    if (d == defaults.end()) {
        throw std::runtime_error("BodyNameRegistry::labels: body " + std::to_string(body) + " is not registered.");
    }
    auto a = aliases.find(index);
    return {.default_name = d->second, .alias = a != aliases.end() ? a->second : std::string{}};
}

void BodyNameRegistry::add_body(int body, const BodyLabels& inherited) {
    assert(std::none_of(names.begin(), names.end(), [body](const auto& entry) {return from_index(entry.second).body == body;}) && "BodyNameRegistry::add_body: body index already exists in the registry.");
    assert(!inherited.default_name.empty() && "BodyNameRegistry::add_body: inherited labels have no default name.");

    // the counter is deliberately left alone: the inherited name was already drawn from it, and advancing it here would burn a number for no body
    add_entity(static_cast<unsigned int>(to_index(body)), inherited.default_name, inherited.alias);
}

void BodyNameRegistry::add_body(int body, const std::string& alias) {
    assert(std::none_of(names.begin(), names.end(), [body](const auto& entry) {return from_index(entry.second).body == body;}) && "BodyNameRegistry::add_body: body index already exists in the registry.");

    // the default name is minted from a monotonic counter rather than from the body index: remove() shifts surviving bodies down while deliberately keeping
    // their names, and elements that create bodies (split, copy) register the freed-up indices afterwards, so an index-derived name would regenerate a string
    // a surviving body already holds. Skip past any number a custom alias has claimed, and never rewind, so no name is ever handed out twice.
    std::string default_name;
    do {default_name = "b" + std::to_string(++bodies_registered);} while (names.contains(default_name));
    add_entity(static_cast<unsigned int>(to_index(body)), default_name, alias);
}

void BodyNameRegistry::add_replica(int body, int isymmetry, int replica) {
    // anchored on the base body's permanent default name for the same reason add_body avoids the index: the index is not stable across body removal/creation,
    // so an index-derived tag could collide with an unrelated body's replicas
    auto base = defaults.find(static_cast<unsigned int>(to_index(body)));
    if (base == defaults.end()) {
        throw std::runtime_error("BodyNameRegistry::add_replica: body " + std::to_string(body) + " is not registered.");
    }

    std::string tag_prefix = base->second + "s" + std::to_string(isymmetry + 1);
    unsigned int index = static_cast<unsigned int>(to_index(body, isymmetry, replica));
    add_entity(index, tag_prefix + "r" + std::to_string(replica), replica == 1 ? tag_prefix : std::string{});
}

void BodyNameRegistry::add_entity(unsigned int index, const std::string& default_name, const std::string& alias) {
    // both slots are validated before either is written, so a rejected registration leaves the registry untouched. Silently dropping a duplicate (as
    // unordered_map::emplace does) would leave one of the two entities addressable by nobody, which is far harder to diagnose than an outright failure.
    bool has_alias = !alias.empty() && alias != default_name;
    if (names.contains(default_name)) {
        throw std::runtime_error("BodyNameRegistry::add_entity: the name \"" + default_name + "\" is already in use.");
    }
    if (has_alias && names.contains(alias)) {
        throw std::runtime_error("BodyNameRegistry::add_entity: the name \"" + alias + "\" is already in use.");
    }

    names.emplace(default_name, index);
    defaults.emplace(index, default_name);
    if (has_alias) {
        names.emplace(alias, index);
        aliases.emplace(index, alias);
    }
}

void BodyNameRegistry::rename(std::string_view old_name, std::string_view new_name) {
    auto it = names.find(std::string{old_name});
    assert(it != names.end() && "BodyNameRegistry::rename: unknown body name.");
    unsigned int index = it->second;

    if (defaults.contains(index)) { // old_name refers to a tracked body; only its alias slot changes
        if (auto a = aliases.find(index); a != aliases.end()) {names.erase(a->second);}
        aliases[index] = std::string{new_name};
        names[std::string{new_name}] = index;
        return;
    }

    // not a tracked body name (e.g. a derived symmetry tag) - plain key swap
    names.erase(it);
    names.emplace(std::string{new_name}, index);
}

void BodyNameRegistry::remove(std::vector<int> body_indices) {
    std::sort(body_indices.begin(), body_indices.end());

    auto shift_of = [&] (int body) -> int {
        auto pos = std::lower_bound(body_indices.begin(), body_indices.end(), body);
        return static_cast<int>(std::distance(body_indices.begin(), pos));
    };
    auto is_erased = [&] (int body) {
        auto pos = std::lower_bound(body_indices.begin(), body_indices.end(), body);
        return pos != body_indices.end() && *pos == body;
    };

    for (auto it = names.begin(); it != names.end(); ) {
        auto sel = from_index(it->second);
        if (is_erased(sel.body)) {it = names.erase(it); continue;}
        if (auto shift = shift_of(sel.body); shift != 0) {it->second = to_index(sel.body - shift, sel.symmetry, sel.replica);}
        ++it;
    }

    auto reindex = [&] (std::unordered_map<unsigned int, std::string>& side) {
        std::unordered_map<unsigned int, std::string> updated;
        for (auto& [index, value] : side) {
            auto sel = from_index(index);
            if (is_erased(sel.body)) {continue;}
            auto shift = shift_of(sel.body);
            unsigned int new_index = shift == 0 ? index : static_cast<unsigned int>(to_index(sel.body - shift, sel.symmetry, sel.replica));
            updated.emplace(new_index, std::move(value));
        }
        side = std::move(updated);
    };
    reindex(defaults);
    reindex(aliases);
}

bool BodyNameRegistry::contains(std::string_view name) const {
    return names.contains(std::string{name});
}

unsigned int BodyNameRegistry::at(std::string_view name) const {
    return names.at(std::string{name});
}

std::string BodyNameRegistry::name(unsigned int index) const {
    if (auto alias = aliases.find(index); alias != aliases.end()) {return alias->second;}
    return defaults.at(index);
}

BodySymmetrySelector BodyNameRegistry::resolve(std::string_view name) const {
    auto it = names.find(std::string{name});
    if (it == names.end()) {
        std::string known;
        for (auto& [n, index] : names) {known += n + " ";}
        throw std::runtime_error("BodyNameRegistry::resolve: Unknown body name \"" + std::string{name} + "\". Known body names are: " + known);
    }
    return from_index(it->second);
}

int BodyNameRegistry::resolve_body(std::string_view name) const {
    auto sel = resolve(name);
    if (sel.symmetry != -1 || sel.replica != 0) {
        throw std::runtime_error("BodyNameRegistry::resolve_body: \"" + std::string{name} + "\" refers to a symmetry replica, not a base body.");
    }
    return sel.body;
}

std::unordered_map<std::string, unsigned int>::const_iterator BodyNameRegistry::begin() const {
    return names.begin();
}

std::unordered_map<std::string, unsigned int>::const_iterator BodyNameRegistry::end() const {
    return names.end();
}

std::vector<BodyNameRegistry::Group> BodyNameRegistry::group_by_index() const {
    std::map<unsigned int, std::vector<std::string>> grouped; // ordered so the result comes out in a stable, index-sorted order
    for (auto& entry : names) {grouped[entry.second].push_back(entry.first);}

    std::vector<Group> result;
    result.reserve(grouped.size());
    for (auto& [index, group_names] : grouped) {
        Group group;
        if (auto d = defaults.find(index); d != defaults.end()) {
            group.default_name = d->second;
            for (auto& name : group_names) {
                if (name != group.default_name) {group.others.push_back(std::move(name));}
            }
        } else {
            group.others = std::move(group_names);
        }
        std::sort(group.others.begin(), group.others.end());
        result.push_back(std::move(group));
    }
    return result;
}

std::vector<std::string> BodyNameRegistry::base_body_names() const {
    std::map<int, std::string> by_body; // ordered by body index, so the result comes out contiguous 0..n-1
    for (const auto& [index, default_name] : defaults) {
        auto sel = from_index(index);
        if (sel.symmetry != -1 || sel.replica != 0) {continue;} // base bodies only, no symmetry replicas
        auto alias = aliases.find(index);
        by_body[sel.body] = (alias != aliases.end()) ? alias->second : default_name;
    }

    std::vector<std::string> result;
    result.reserve(by_body.size());
    for (auto& [body, name] : by_body) {result.push_back(std::move(name));}
    return result;
}
