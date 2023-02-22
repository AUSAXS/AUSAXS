#pragma once

#include <rigidbody/selection/BodySelectStrategy.h>

#include <random>

namespace rigidbody {
	/**
	 * @brief RandomSelect
	 * 
	 * This selection strategy randomly selects a Body. 
	 */
	class SequentialSelect : public BodySelectStrategy {
		public: 
			/**
			 * @brief Constructor.
			 */
			SequentialSelect(const RigidBody* rigidbody);

			/**
			 * @brief Destructor.
			 */
			~SequentialSelect() override;

			/**
			 * @brief Get the index of the next body to be transformed. 
			 */
			unsigned int next() override;

		private:
			std::mt19937 generator;                          // The random number generator. 
			std::uniform_int_distribution<int> distribution; // The random number distribution. 
	};	
}