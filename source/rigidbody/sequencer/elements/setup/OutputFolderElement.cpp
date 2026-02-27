// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/OutputFolderElement.h>
#include <rigidbody/sequencer/elements/setup/SetupElement.h>
#include <rigidbody/sequencer/detail/parse_error.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <settings/GeneralSettings.h>
#include <utility/Logging.h>

using namespace ausaxs::rigidbody::sequencer;

OutputFolderElement::OutputFolderElement(observer_ptr<Sequencer> owner, const io::Folder& folder, Mode mode) : owner(owner) {
    std::string prefix = "";
    switch (mode) {
        case Mode::RELATIVE_TERMINAL:
            break;
        case Mode::RELATIVE_CONFIG:
            prefix = owner->setup()._get_config_folder() + "/";
            break;
    }

    settings::general::output = prefix + folder.path() + (folder.path().back() == '/' ? "" : "/");
    if (settings::general::verbose) {
        logging::log("OutputFolderElement: Setting output folder to \"" + settings::general::output + "\".");
    }
}

void OutputFolderElement::run() {}

std::unique_ptr<GenericElement> OutputFolderElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    enum class Args {path, mode};
    static std::unordered_map<Args, std::vector<std::string>> valid_args = {
        {Args::path, {"path", "folder", "anonymous"}},
        {Args::mode, {"mode", "relative"}}
    };
    auto path = args.get<std::string>(valid_args[Args::path]);
    auto mode = args.get<std::string>(valid_args[Args::mode], "relative_terminal");

    if (args.named.size() > 2) {throw except::parse_error("output", "Too many named arguments. Expected at most 2, but got " + std::to_string(args.named.size()) + ".");}
    if (args.inlined.size() == 1) {
        if (path.found) {throw except::parse_error("output", "Unexpected named argument \"" + args.named.begin()->first + "\".");}
        path.value = args.inlined[0];
        path.found = true;
    } else if (args.inlined.size() > 1) {
        throw except::parse_error("output", "Too many inline arguments. Expected at most 1, but got " + std::to_string(args.inlined.size()) + ".");
    }
    if (!path.found) {throw except::parse_error("output", "Missing required argument \"path\".");}

    if (mode.value == "relative" || mode.value == "relative_terminal") {
        return std::make_unique<OutputFolderElement>(owner->_get_sequencer(), io::Folder(path.value), OutputFolderElement::Mode::RELATIVE_TERMINAL);
    } else if (mode.value == "relative_config") {
        return std::make_unique<OutputFolderElement>(owner->_get_sequencer(), io::Folder(path.value), OutputFolderElement::Mode::RELATIVE_CONFIG);
    } else {
        throw except::parse_error("output", "Invalid argument for \"mode\": \"" + mode.value + "\". Expected one of {absolute, relative, relative_config}.");
    }
}