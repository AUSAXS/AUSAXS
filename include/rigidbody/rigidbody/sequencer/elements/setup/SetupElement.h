// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <io/ExistingFile.h>
#include <io/Folder.h>
#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <rigidbody/sequencer/detail/BodyNameRegistry.h>
#include <rigidbody/sequencer/elements/GenericElement.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>
#include <string_view>
#include <utility/observer_ptr.h>

#include <memory>
#include <string>
#include <vector>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Holds state shared by setup-phase sequence elements.
     */
    class SetupElement {
        public:
            SetupElement(observer_ptr<Sequencer> sequencer);
            SetupElement(observer_ptr<Sequencer> sequencer, io::ExistingFile saxs);

            /**
             * @brief Get the name identifiers of all loaded bodies.
             */
            detail::BodyNameRegistry& _body_name_registry();

            /**
             * @brief Resolve a name to the (body, symmetry, replica) selector it refers to. Accepts any known name, including a symmetry replica's tag.
             */
            detail::BodySymmetrySelector _get_body_index(std::string_view name) const;

            /**
             * @brief Resolve a name to a base body's index. Throws if the name refers to a symmetry replica.
             */
            int _get_body(std::string_view name) const;

            /**
             * @brief Set the currently active body for the setup.
             */
            void _set_active_body(observer_ptr<Rigidbody> body);

            /**
             * @brief Get the location of the configuration folder.
             *        This may be empty if no configuration file was loaded. 
             */
            std::string _get_config_folder() const;

            /**
             * @brief Set the location of the configuration folder.
             *        This is used to resolve relative paths in the configuration file.
             */
            void _set_config_folder(const io::Folder& folder);

            /**
             * @brief Set the location of the SAXS measurement data.
             */
            void _set_saxs_path(const io::ExistingFile& saxs);

            /**
             * @brief Get the path to the SAXS measurement data.
             */
            const io::ExistingFile& _get_saxs_path() const;

            /**
             * @brief Get the elements that are part of this setup.
             */
            std::vector<std::unique_ptr<GenericElement>>& _get_elements();

        private:
            observer_ptr<Sequencer> owner;
            detail::BodyNameRegistry body_names;
            observer_ptr<Rigidbody> active_body;
            io::Folder config_folder;
            io::ExistingFile saxs_path;
            std::vector<std::unique_ptr<GenericElement>> elements;
    };
}