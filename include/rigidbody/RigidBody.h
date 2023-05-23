#pragma once

#include <data/Protein.h>

#include <memory>

namespace grid {class Grid;}
namespace fitter {
	class Fit;
	class LinearFitter;
}
namespace rigidbody {
	namespace detail {
		struct BestConf {
			BestConf();
			BestConf(std::shared_ptr<grid::Grid> grid, std::vector<Water> waters, double chi2) noexcept;
			~BestConf();
			std::shared_ptr<grid::Grid> grid;
			std::vector<Water> waters;
			double chi2;	
		};

	}

	namespace selection {class BodySelectStrategy;}
	class ConstraintManager;
	class TransformStrategy;
	class ParameterGenerationStrategy;
	class RigidBody : public Protein {
		public:
			RigidBody(Protein&& protein);

			RigidBody(const Protein& protein);

			virtual ~RigidBody();

			/**
			 * @brief Perform a rigid-body optimization for this structure. 
			 */
			std::shared_ptr<fitter::Fit> optimize(const std::string& measurement_path);

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

			std::shared_ptr<ConstraintManager> constraints;
		protected:
			std::shared_ptr<fitter::Fit> calibration = nullptr;
			std::unique_ptr<rigidbody::selection::BodySelectStrategy> body_selector;
			std::unique_ptr<TransformStrategy> transform;
			std::unique_ptr<ParameterGenerationStrategy> parameter_generator;
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