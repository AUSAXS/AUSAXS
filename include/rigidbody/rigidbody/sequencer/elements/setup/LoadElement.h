// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/SequencerFwd.h>
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
            LoadElement(observer_ptr<Sequencer> owner, const std::vector<std::string>& paths, const std::vector<std::string>& body_names = {}, const std::string& saxs_path = {});

            /**
             * @brief Load multiple bodies from a single file, separated at the designated indices. 
             */
            LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<int>& split, const std::vector<std::string>& body_names = {}, const std::string& saxs_path = {});

            /**
             * @brief Load multiple bodies from a single file, separated by the chainID. 
             */
            LoadElement(observer_ptr<Sequencer> owner, const std::string& path, const std::vector<std::string>& body_names = {}, const std::string& saxs_path = {});

            ~LoadElement() override;

            void run() override;
        
        private:
            observer_ptr<Sequencer> owner;
            std::unique_ptr<rigidbody::Rigidbody> rigidbody;

            std::vector<std::string> load_wildcarded(const std::string& path);
            std::pair<std::string, bool> lookup_file(const std::string& path);
    };
}