// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/OutputFolderElement.h>
#include <rigidbody/sequencer/elements/setup/SetupElement.h>
#include <rigidbody/sequencer/detail/ArgumentHelper.h>
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

    // the output folder is later concatenated with file names, so it must end in a separator.
    // '\' counts as one too, since Windows paths are written with it
    auto path = folder.path();
    bool terminated = !path.empty() && (path.back() == '/' || path.back() == '\\');
    settings::general::output = prefix + path + (terminated ? "" : "/");
    if (settings::general::verbose) {
        logging::log("OutputFolderElement: Setting output folder to \"" + settings::general::output + "\".");
    }
}

void OutputFolderElement::run() {}

namespace {
    enum class Args {path, mode};
    std::unordered_map<Args, std::vector<std::string>> args_map = {
        {Args::path, {"path", "folder", "anonymous"}},
        {Args::mode, {"mode", "relative"}}
    };
}

std::vector<std::string> OutputFolderElement::_valid_arguments() {
    static auto map = detail::get_arg_names<Args>(args_map);
    return map;
}

InlineSignature OutputFolderElement::_valid_inline_arguments() {
    return {.names = {"path"}, .min = 0, .max = 1};
}

std::unique_ptr<GenericElement> OutputFolderElement::_parse(observer_ptr<LoopElement> owner, ParsedArgs&& args) {
    auto path = args.get<std::string>(args_map[Args::path]);
    auto mode = args.get<std::string>(args_map[Args::mode], "relative_terminal");

    if (!args.inlined.empty()) {
        path.value = args.inlined[0];
        path.found = true;
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