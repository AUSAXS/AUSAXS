// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <data/DataFwd.h>
#include <io/ExistingFile.h>
#include <io/Folder.h>
#include <rigidbody/RigidbodyFwd.h>
#include <rigidbody/controller/IController.h>
#include <rigidbody/sequencer/detail/BodyNameRegistry.h>
#include <rigidbody/sequencer/elements/LoopElement.h>
#include <rigidbody/sequencer/elements/setup/BodySymmetrySelector.h>

#include <memory>
#include <string>
#include <string_view>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief The top-level driver of a rigidbody optimization.
     *
     * The optimization is described declaratively by a sequence script. A Sequencer is itself a
     * LoopElement, the root of the parsed element tree, and additionally owns the state shared by
     * the setup-phase elements: the body name registry, the SAXS data path, and the configuration
     * folder. The various `_get_*` accessors give nested elements access to that state. execute()
     * runs the parsed sequence and returns the resulting fit.
     */
    class Sequencer : public LoopElement {
        public:
            Sequencer();
            Sequencer(io::ExistingFile saxs);
            ~Sequencer() override;

            /**
             * @brief Execute the sequencer.
             * @return The result of the fit.
             */
            std::shared_ptr<fitter::FitResult> execute() override;

            /**
             * @brief Disabled. Sequencer must always be invoked via execute() to ensure proper setup.
             *        Calling run() directly skips rigidbody and controller initialization.
             */
            void run() override;

            /**
             * @brief Get the Rigidbody object.
             */
            observer_ptr<Rigidbody> _get_rigidbody() const override;
            void _set_rigidbody(observer_ptr<Rigidbody> rigidbody);

            /**
             * @brief Get the molecule object.
             * 
             * This is a convenience method to access the molecule from the Rigidbody.
             */
            observer_ptr<data::Molecule> _get_molecule() const override;

            /**
             * @brief Get the top Sequencer object.
             */
            observer_ptr<const Sequencer> _get_sequencer() const override;
            observer_ptr<Sequencer> _get_sequencer() override;

            /**
             * @brief Get the current rigidbody controller.
             */
            observer_ptr<controller::IController> _get_controller() const;

            /**
             * @brief Get the best configuration found so far.
             */
            observer_ptr<rigidbody::detail::MoleculeTransformParametersAbsolute> _get_best_conf() const override;

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

        private:
            observer_ptr<Rigidbody> rigidbody;
            detail::BodyNameRegistry body_names;
            io::Folder config_folder;
            io::ExistingFile saxs_path;
    };
}