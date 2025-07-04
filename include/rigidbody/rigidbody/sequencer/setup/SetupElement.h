// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/sequencer/SequencerFwd.h>

#include <rigidbody/sequencer/GenericElement.h>
#include <rigidbody/sequencer/LoopElementCallback.h>
#include <utility/observer_ptr.h>
#include <io/ExistingFile.h>
#include <io/Folder.h>

#include <string>
#include <vector>
#include <unordered_map>
#include <functional>

namespace ausaxs::rigidbody::sequencer {
    /**
     * @brief Set up the optimization problem.
     *        Once any method from the LoopElementCallback is called, the setup is considered complete.
     */
    class SetupElement : public LoopElementCallback {
        public:
            SetupElement(observer_ptr<Sequencer> owner);
            SetupElement(observer_ptr<Sequencer> owner, io::ExistingFile saxs);
            virtual ~SetupElement() = default;

            /**
             * @brief Set the overlap function for evaluating the overlap penalty. 
             *        The function is multiplied onto the distance histogram, with the sum of the product being the chi2 penalty. 
             *
             * @param func The monotonically decreasing weight function. 
             */
            SetupElement& set_overlap_function(std::function<double(double)> func);

            /**
             * @brief Load a body from a file. 
             * 
             * @param path Path to the file.
             * @param body_names Optional names for the bodies contained in the file.
             */
            SetupElement& load(const std::vector<std::string>& path, const std::vector<std::string>& body_names = {});

            /**
             * @brief Load a SAXS file. 
             */
            SetupElement& load(const io::ExistingFile& saxs);

            /**
             * @brief Load an existing rigidbody. 
             */
            SetupElement& load_existing(observer_ptr<RigidBody> rigidbody);

            /**
             * @brief Create a distance constraint between the two bodies at the specified atoms.
             * 
             * @param ibody1 Index of the first body.
             * @param ibody2 Index of the second body.
             * @param iatom1 Index of the first atom.
             * @param iatom2 Index of the second atom. 
             */
            SetupElement& distance_constraint(unsigned int ibody1, unsigned int ibody2, unsigned int iatom1, unsigned int iatom2);

            /**
             * @brief Create a distance constraint between the two bodies at the specified atoms.
             * 
             * @param body1 Name of the first body.
             * @param body2 Name of the second body.
             * @param iatom1 Index of the first atom.
             * @param iatom2 Index of the second atom. 
             */
            SetupElement& distance_constraint(const std::string& body1, const std::string& body2, unsigned int iatom1, unsigned int iatom2);

            /**
             * @brief Create a distance constraint between the two bodies at the closest atomic pair. 
             * 
             * @param body1 Index of the first body.
             * @param body2 Index of the second body.
             */
            SetupElement& distance_constraint_closest(unsigned int ibody1, unsigned int ibody2);

            /**
             * @brief Create a distance constraint between the two bodies at the closest atomic pair. 
             * 
             * @param body1 Name of the first body.
             * @param body2 Name of the second body.
             */
            SetupElement& distance_constraint_closest(const std::string& ibody1, const std::string& ibody2);

            /**
             * @brief Create a distance constraint between the two bodies at the center of masses. 
             * 
             * @param body1 Index of the first body.
             * @param body2 Index of the second body.
             */
            SetupElement& distance_constraint_center_mass(unsigned int ibody1, unsigned int ibody2);

            /**
             * @brief Create a distance constraint between the two bodies at the center of masses. 
             * 
             * @param body1 Name of the first body.
             * @param body2 Name of the second body.
             */
            SetupElement& distance_constraint_center_mass(const std::string& body1, const std::string& body2);

            /**
             * @brief Create a fixed constraint for the currently active body. 
             */
            SetupElement& fixed_constraint();

            /**
             * @brief Automatically create sequential distance constraints between all bodies. 
             *        Given the three bodies A, B, C, the constraints will be A-B, and B-C.
             */
            SetupElement& generate_linear_constraints();

            /**
             * @brief Automatically create distance constraints between all bodies that are close to each other.
             */
            SetupElement& generate_volumetric_constraints();

            /**
             * @brief Get the name identifiers of all loaded bodies.
             */
            std::unordered_map<std::string, unsigned int>& _get_body_names();

            /**
             * @brief Set the currently active body for the setup.
             */
            void _set_active_body(observer_ptr<RigidBody> body);

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

        protected:
            io::ExistingFile saxs_path;
            std::vector<std::unique_ptr<GenericElement>> elements;

        private:
            std::unordered_map<std::string, unsigned int> body_names;
            observer_ptr<RigidBody> active_body = nullptr;
            io::Folder config_folder;
    };
}