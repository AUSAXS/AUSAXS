// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/detail/BodyNameRegistry.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>

#include <algorithm>
#include <cassert>
#include <stdexcept>

using namespace ausaxs::rigidbody::sequencer::detail;

void BodyNameRegistry::add_body(int body, const std::string& alias) {
    assert(!has_body(body) && "BodyNameRegistry::add_body: body index already exists in the registry.");
    std::string default_name;
    do {default_name = "b" + std::to_string(++bodies_registered);} while (lookup.contains(default_name));
    add_entity(to_index(body), {.default_name = std::move(default_name), .alias = alias});
}

void BodyNameRegistry::add_body(int body, Entry inherited) {
    assert(!has_body(body) && "BodyNameRegistry::add_body: body index already exists in the registry.");
    assert(!inherited.default_name.empty() && "BodyNameRegistry::add_body: inherited entry has no default name.");

    // the counter is deliberately left alone: the inherited name was already drawn from it, and advancing it here would burn a number for no body
    add_entity(to_index(body), std::move(inherited));
}

void BodyNameRegistry::add_replica(int body, int isymmetry, int replica) {
    auto base = entries.find(to_index(body));
    if (base == entries.end()) {
        throw std::runtime_error("BodyNameRegistry::add_replica: body " + std::to_string(body) + " is not registered.");
    }

    std::string tag_prefix = base->second.default_name + "s" + std::to_string(isymmetry + 1);
    add_entity(to_index(body, isymmetry, replica), {
        .default_name = tag_prefix + "r" + std::to_string(replica),
        .alias = replica == 1 ? tag_prefix : std::string{} // the first replica doubles as the symmetry itself, and gets the shorter name for it
    });
}

void BodyNameRegistry::add_entity(int index, Entry entry) {
    // both names are validated before either is written, so a rejected registration leaves the registry untouched
    if (entry.alias == entry.default_name) {entry.alias.clear();} // not a distinct name, so not an alias
    if (lookup.contains(entry.default_name)) {
        throw std::runtime_error("BodyNameRegistry::add_entity: the name \"" + entry.default_name + "\" is already in use.");
    }
    if (!entry.alias.empty() && lookup.contains(entry.alias)) {
        throw std::runtime_error("BodyNameRegistry::add_entity: the name \"" + entry.alias + "\" is already in use.");
    }

    lookup.emplace(entry.default_name, index);
    if (!entry.alias.empty()) {lookup.emplace(entry.alias, index);}
    entries.emplace(index, std::move(entry));
}

void BodyNameRegistry::rename(std::string_view old_name, std::string_view new_name) {
    auto it = lookup.find(std::string{old_name});
    if (it == lookup.end()) {
        throw std::runtime_error("BodyNameRegistry::rename: unknown body name \"" + std::string{old_name} + "\".");
    }
    int index = it->second;

    if (auto taken = lookup.find(std::string{new_name}); taken != lookup.end() && taken->second != index) {
        throw std::runtime_error("BodyNameRegistry::rename: the name \"" + std::string{new_name} + "\" is already in use.");
    }

    Entry& entry = entries.at(index);
    if (!entry.alias.empty()) {lookup.erase(entry.alias);}
    if (new_name == entry.default_name) { // renaming an entity back to its permanent name just drops the alias
        entry.alias.clear();
        return;
    }
    entry.alias = std::string{new_name};
    lookup.emplace(entry.alias, index);
}

void BodyNameRegistry::remove(std::vector<int> body_indices) {
    std::sort(body_indices.begin(), body_indices.end());

    // every surviving body shifts down by the number of erased bodies preceding it
    auto shift_of = [&] (int body) {
        return static_cast<int>(std::distance(body_indices.begin(), std::lower_bound(body_indices.begin(), body_indices.end(), body)));
    };
    auto is_erased = [&] (int body) {
        auto pos = std::lower_bound(body_indices.begin(), body_indices.end(), body);
        return pos != body_indices.end() && *pos == body;
    };

    std::map<int, Entry> kept;
    for (auto& [index, entry] : entries) {
        auto sel = from_index(index);
        if (is_erased(sel.body)) {continue;} // drops the body's replicas along with it, since they encode the same body
        kept.emplace(to_index(sel.body - shift_of(sel.body), sel.symmetry, sel.replica), std::move(entry));
    }
    entries = std::move(kept);

    // the names themselves are untouched by the shift - a body keeps its identity across the erase - so the lookup only needs its indices refreshed
    rebuild_lookup();
}

void BodyNameRegistry::rebuild_lookup() {
    lookup.clear();
    for (const auto& [index, entry] : entries) {
        lookup.emplace(entry.default_name, index);
        if (!entry.alias.empty()) {lookup.emplace(entry.alias, index);}
    }
}

bool BodyNameRegistry::has_body(int body) const {
    auto it = entries.lower_bound(to_index(body)); // the base body sorts before its own replicas, which sort before the next body
    return it != entries.end() && from_index(it->first).body == body;
}

bool BodyNameRegistry::contains(std::string_view name) const {
    return lookup.contains(std::string{name});
}

BodySymmetrySelector BodyNameRegistry::resolve(std::string_view name) const {
    auto it = lookup.find(std::string{name});
    if (it == lookup.end()) {
        std::string known;
        for (const auto& [index, entry] : entries) { // listed in index order, so the suggestion reads in the same order as the bodies
            known += entry.display_name() + " ";
        }
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

const BodyNameRegistry::Entry& BodyNameRegistry::entry(int index) const {
    return entries.at(index);
}

const std::map<int, BodyNameRegistry::Entry>& BodyNameRegistry::all() const {
    return entries;
}

std::vector<std::string> BodyNameRegistry::base_body_names() const {
    std::vector<std::string> result;
    for (const auto& [index, entry] : entries) { // already ordered by encoded index, so the base bodies come out in body order
        auto sel = from_index(index);
        if (sel.symmetry == -1 && sel.replica == 0) {result.push_back(entry.display_name());}
    }
    return result;
}
