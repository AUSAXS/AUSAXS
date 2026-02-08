// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#include <rigidbody/sequencer/elements/setup/OutputFolderElement.h>
#include <rigidbody/sequencer/elements/setup/SetupElement.h>
#include <rigidbody/sequencer/Sequencer.h>
#include <rigidbody/Rigidbody.h>
#include <settings/GeneralSettings.h>

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
        std::cout << "OutputFolderElement: Setting output folder to \"" << settings::general::output << "\"" << std::endl;
    }
}

void OutputFolderElement::run() {}