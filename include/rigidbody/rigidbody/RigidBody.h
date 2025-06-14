// SPDX-License-Identifier: LGPL-3.0-or-later
// Author: Kristian Lytje

#pragma once

#include <rigidbody/detail/RigidbodyInternalFwd.h>
#include <rigidbody/sequencer/SequencerFwd.h>
#include <fitter/FitterFwd.h>
#include <grid/GridFwd.h>
#include <data/DataFwd.h>
#include <data/Molecule.h>

#include <memory>

namespace ausaxs::rigidbody {
	class RigidBody : public data::Molecule {
		friend rigidbody::sequencer::Sequencer;
		public:
			template <typename... Args, typename = decltype(Molecule(std::declval<Args>()...))>
			RigidBody(Args&&... args) : Molecule(std::forward<Args>(args)...) {initialize();}
			virtual ~RigidBody() override;

			/**
			 * @brief Apply a calibration to this rigid body. 
			 * 
			 * This will fix the solvent scattering density to the fitted value.
			 */
			void apply_calibration(std::unique_ptr<fitter::FitResult> calibration);

			/**
			 * @brief Update the given fitter with the current rigid body parameters.
			 */
			void update_fitter();

			/**
			 * @brief Get the constraint manager for this rigid body.
			 */
			std::shared_ptr<constraints::ConstraintManager> get_constraint_manager() const;

			/**
			 * @brief Get the fitter for this rigid body. Note that this is a ConstrainedFitter, and will thus include chi2 contributions from the constraints. 
			 */
			std::shared_ptr<fitter::SmartFitter> get_fitter() const;

			/**
			 * @brief Create a new fitter for this rigid body. This fitter will not include any constraints.
			 */
			std::unique_ptr<fitter::SmartFitter> get_unconstrained_fitter(const io::ExistingFile& saxs) const;

			void set_constraint_manager(std::shared_ptr<rigidbody::constraints::ConstraintManager> constraints);

			void set_body_select_manager(std::shared_ptr<rigidbody::selection::BodySelectStrategy> body_selector);

			void set_transform_manager(std::shared_ptr<rigidbody::transform::TransformStrategy> transform);

			void set_parameter_manager(std::shared_ptr<rigidbody::parameter::ParameterGenerationStrategy> parameters);

		private:
			std::shared_ptr<constraints::ConstraintManager> constraints = nullptr;
			std::shared_ptr<fitter::FitResult> calibration = nullptr;
			std::shared_ptr<selection::BodySelectStrategy> body_selector;
			std::shared_ptr<transform::TransformStrategy> transform;
			std::shared_ptr<parameter::ParameterGenerationStrategy> parameter_generator;
			std::shared_ptr<fitter::SmartFitter> fitter;

			/**
			 * @brief Perform an optimization step.
			 * 
			 * @return True if a better configuration was found, false otherwise.
			 */
			bool optimize_step(detail::BestConf& best);

			/**
			 * @brief Prepare the fitter for this rigidbody.
			 */
			void prepare_fitter(const io::ExistingFile& measurement_path); 

			/**
			 * @brief Small initialization function.
			 */
			void initialize();
	};
}