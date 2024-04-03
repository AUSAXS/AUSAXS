#pragma once

#include <rigidbody/detail/RigidbodyInternalFwd.h>
#include <hydrate/GridFwd.h>
#include <fitter/FitterFwd.h>
#include <data/Molecule.h>

#include <memory>

namespace rigidbody {
	namespace detail {
		struct BestConf {
			BestConf();
			BestConf(std::shared_ptr<grid::Grid> grid, std::vector<data::record::Water> waters, double chi2) noexcept;
			~BestConf();
			std::shared_ptr<grid::Grid> grid;
			std::vector<data::record::Water> waters;
			double chi2;	
		};
	}

	class RigidBody : public data::Molecule {
		public:
			RigidBody(data::Molecule&& protein);

			RigidBody(const data::Molecule& protein);

			virtual ~RigidBody();

			/**
			 * @brief Perform a rigid-body optimization for this structure. 
			 */
			std::shared_ptr<fitter::Fit> optimize(const io::ExistingFile& measurement_path);

			std::shared_ptr<fitter::Fit> optimize_sequence(const io::ExistingFile& measurement_path);

			/**
			 * @brief Apply a calibration to this rigid body. 
			 * 
			 * This will fix the solvent scattering density to the fitted value.
			 */
			void apply_calibration(std::shared_ptr<fitter::Fit> calibration);

			/**
			 * @brief Update the given fitter with the current rigid body parameters.
			 */
			void update_fitter(std::shared_ptr<fitter::LinearFitter> fitter);

			/**
			 * @brief Get the constraint manager for this rigid body.
			 */
			std::shared_ptr<constraints::ConstraintManager> get_constraint_manager() const;

		protected:
			std::shared_ptr<constraints::ConstraintManager> constraints = nullptr;
			std::shared_ptr<fitter::Fit> calibration = nullptr;
			std::unique_ptr<selection::BodySelectStrategy> body_selector;
			std::unique_ptr<transform::TransformStrategy> transform;
			std::unique_ptr<parameter::ParameterGenerationStrategy> parameter_generator;
			std::shared_ptr<fitter::LinearFitter> fitter;

			/**
			 * @brief Perform an optimization step.
			 * 
			 * @return True if a better configuration was found, false otherwise.
			 */
			bool optimize_step(detail::BestConf& best);

			/**
			 * @brief Prepare the fitter for this rigidbody.
			 */
			void prepare_fitter(const std::string& measurement_path); 

			/**
			 * @brief Small initialization function.
			 */
			void initialize();
	};
}