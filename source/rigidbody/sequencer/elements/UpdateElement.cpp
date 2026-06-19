// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/UpdateElement.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <data/Molecule.h>
#include <data/symmetry/MoleculeSymmetryFacade.h>
#include <data/detail/SimpleBody.h>
#include <utility/MultiThreading.h>
#include <utility/Console.h>

#include <string>
#include <utility>
#include <mutex>

using namespace ausaxs::rigidbody::sequencer;

UpdateElement::UpdateElement(observer_ptr<LoopElement> owner) : LoopElementCallback(owner) {
    // warn (once, at parse time) if nobody is going to read what we publish
    if (!live_consumer_connected) {
        console::print_warning("update: no live consumer is connected, disabling updates.");
    }

    // start each sequence with an empty buffer so stale frames from a previous run aren't served
    lock();
    x.clear(); y.clear(); z.clear();
    version = 0;
    unlock();
}

UpdateElement::~UpdateElement() = default;

void UpdateElement::run() {
    // nothing is reading the live structure, so don't spend cycles building and publishing it
    if (!live_consumer_connected) {return;}

    auto structure = owner->_get_molecule()->symmetry().explicit_structure();
    utility::multi_threading::get_global_pool()->detach_task(
        [atoms = std::move(structure.atoms)]() {
            lock();
            std::size_t n = atoms.size();
            x.resize(n); y.resize(n); z.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                x[i] = atoms[i].x();
                y[i] = atoms[i].y();
                z[i] = atoms[i].z();
            }
            ++version;
            unlock();
        }
    );
}

std::mutex mutex;
void UpdateElement::lock() {
    mutex.lock();
}

void UpdateElement::unlock() {
    mutex.unlock();
}

std::vector<std::string> UpdateElement::_valid_arguments() {
    return {"structure"};
}

std::unique_ptr<GenericElement> UpdateElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    if (!args.named.empty()) {throw except::parse_error("update", "Unexpected named argument.");}
    if (args.inlined.size() != 1) {throw except::parse_error("update", "Expected a single inline argument, e.g. \"update structure\".");}
    if (std::string(args.inlined[0]) != "structure") {
        throw except::parse_error("update", "Unsupported update target \"" + std::string(args.inlined[0]) + "\"; only \"structure\" is supported.");
    }
    return std::make_unique<UpdateElement>(owner);
}
