// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/detail/ParsedArgs.h>
#include <rigidbody/RigidbodyFwd.h>
#include <utility/observer_ptr.h>
#include <io/ExistingFile.h>

#include <string>
#include <vector>
#include <memory>

namespace ausaxs::rigidbody::sequencer {
    class LoadElement : public GenericElement {
        public:
            /**
             * @brief Load multiple bodies from multiple files. One body is loaded from each file.
             */
            LoadElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& paths, const std::vector<std::string>& body_names = {});

            /**
             * @brief Load multiple bodies from a single file, separated at the designated indices. 
             */
            LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<int>& split, const std::vector<std::string>& body_names = {});

            /**
             * @brief Load multiple bodies from a single file, separated by the chainID. 
             */
            LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<std::string>& body_names = {});

            ~LoadElement() override;

            void run() override;

            static std::vector<std::string> _valid_arguments();
            static std::unique_ptr<GenericElement> _parse(observer_ptr<LoopElement> owner, ParsedArgs&& args);

        protected:
            observer_ptr<Sequencer> owner;
            std::unique_ptr<rigidbody::Rigidbody> rigidbody;

            // the actual files loaded (after wildcard expansion and path lookup), in body order;
            // retained so subclasses (e.g. LoadElementWrapper) can re-read them for preview metadata
            std::vector<std::string> resolved_paths;

            std::vector<std::string> load_wildcarded(const std::string& path);
            std::pair<std::string, bool> lookup_file(const std::string& path);
    };
}