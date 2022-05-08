#pragma once

// forwards declaration
class Protein;

#include <data/Atom.h>
#include <data/Body.h>
#include <data/StateManager.h>
#include <histogram/ScatteringHistogram.h>
#include <histogram/Histogram.h>

namespace hist {
		/**
		 * @brief A compact vector representation of the coordinates and weight of all atoms in a body. 
		 *        The idea is that by only extracting the absolute necessities for the distance calculation, more values can be stored
		 *        in the cache at any given time. This is meant as a helper class to DistanceCalculator.
		 */
	struct CompactCoordinates {
		struct Data {
			Data() {}
			Data(const Vector3& v, float w) : x(v.x()), y(v.y()), z(v.z()), w(w) {}
			float x, y, z, w;
		};

		CompactCoordinates() {}

		/**
		 * @brief Extract the necessary coordinates and weights from a body. 
		 */
		CompactCoordinates(const Body& body);

		/**
		 * @brief Extract the necessary coordinates and weights from a vector of hydration atoms. 
		 */
		CompactCoordinates(const vector<Hetatom>& atoms);

		size_t size;
		vector<Data> data;
	};

	/**
	 * @brief Simple data containers defined for clarity.  
	 */
	typedef Histogram PartialHistogram;
	typedef Histogram HydrationHistogram;

	/**
	 * @brief We also define the MasterHistogram type, which is identical to a PartialHistogram. 
	 *        We do this to make += and -= well-defined operations. 
	 */
	class MasterHistogram : public Histogram {
		public: 
			MasterHistogram() {}

			/**
			 * @brief Create a new Master Histogram. 
			 * @param p The current histogram. 
			 * @param p_base The constant, unchanging part of the histogram. 
			 */
			MasterHistogram(const vector<double>& p_base, const Axis& axis);

			/**
			 * @brief Add a PartialHistogram to the MasterHistogram. 
			 */
			MasterHistogram& operator+=(const PartialHistogram& rhs);

			/**
			 * @brief Subtract a PartialHistogram from the MasterHistogram. We have to use a lambda since the standard std::minus would
			 *        reverse the order of the entries.
			 */
			MasterHistogram& operator-=(const PartialHistogram& rhs);

			// The base part of the histogram which will never change. This contains all internal distances between atoms in each individual body.
			Histogram base;
	};

	/**
	 * The basic idea is that we have a bunch of partial histograms (contained in @a partials), which combined represents the total scattering histogram. 
	 * As an example, if we had 4 bodies, it would look something like this:
	 * 4       x
	 * 3     x
	 * 2   x
	 * 1 x
	 *   1 2 3 4
	 * The self-correlation partials are marked with an 'x'. They are constant and are thus precalculated when this class is initialized. 
	 * The upper and lower triangle are symmetric, and we can thus just calculate one of them and double the result. After all partials are initially
	 * generated, this class recalculates them whenever a body has changed. If body 2 is moved, the partials (1, 2), (2, 3), and (2, 4) must be recalculated. 
	 * 
	 * This is further complicated by the presence of the hydration layer. Since this does not belong to any individual body, it can be viewed as 
	 * a simple extension to the above example, so we now have {1, 2, 3, 4, H}. 
	 */

	/**
	 * @brief A smart distance calculator which efficiently calculates the scattering histogram.
	 */
	class PartialHistogramManager {
		public:
			PartialHistogramManager(Protein* protein); 

			/**
			 * @brief Initialize this object. The internal distances between atoms in each body is constant and cannot change. 
			 *        They are unaffected by both rotations and translations, and so we precalculate them. 
			 */
			void initialize();

			/**
			 * @brief Calculate only the total scattering histogram. 
			 */
			Histogram calculate();

			/**
			 * @brief Calculate all contributions to the scattering histogram. 
			 */
			ScatteringHistogram calculate_all();

			/**
			 * @brief Get a signalling object for signalling a change of state. 
			 *        Each body is supposed to hold one of these, and trigger it when they change state. 
			 */
			std::shared_ptr<StateManager::BoundSignaller> get_probe(unsigned int i);

			/**
			 * @brief Signal that the hydration layer was modified. 
			 *        This is supposed to be used only by the Protein class, which has direct access to this object. Thus a signalling object is unnecessary. 
			 */
			void signal_modified_hydration_layer();

			const StateManager& get_state_manager() const;

			StateManager& get_state_manager();

		private:
			const size_t size;                            // number of managed bodies
			StateManager statemanager;                    // a helper which keeps track of state changes in each body
			vector<CompactCoordinates> coords_p;          // a compact representation of the relevant data from the managed bodies
			CompactCoordinates coords_h;                  // a compact representation of the hydration data
			Protein* protein;                             // pointer to the parent Protein

			// histogram data
			MasterHistogram master;                       // the current total histogram
			vector<vector<PartialHistogram>> partials_pp; // the partial histograms
			vector<HydrationHistogram> partials_hp;       // the partial hydration-atom histograms
			HydrationHistogram partials_hh;               // the partial histogram for the hydration layer

			/**
			 * @brief Calculate the atom-atom distances between body @a index and all others. 
			 */
			void calc_pp(unsigned int index);

			/**
			 * @brief Calculate the atom-atom distances between body @a n and @a m. 
			 */
			void calc_pp(unsigned int n, unsigned int m);

			/**
			 * @brief Calculate the hydration-atom distances between the hydration layer and body @a index.
			 */
			void calc_hp(unsigned int index);

			/**
			 * @brief Calculate the hydration-hydration distances. 
			 */
			void calc_hh();

			/**
			 * @brief Calculate the self-correlation of a body. 
			 */
			void calc_self_correlation(unsigned int index);
	};
}