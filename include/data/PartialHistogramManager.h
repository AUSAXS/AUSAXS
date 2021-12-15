#include "data/Atom.h"
#include "data/Body.h"
#include "data/StateManager.h"
#include "ScatteringHistogram.h"

/**
 * @brief A compact vector representation of the coordinates and weight of all atoms in a body. 
 *        The idea is that by only extracting the absolute necessities for the distance calculation, more values can be stored
 *        in the cache at any given time. This is meant as a helper class to DistanceCalculator.
 */
struct CompactCoordinates {
    /**
     * @brief Extract the necessary coordinates and weights from a body. 
     */
    CompactCoordinates(const Body& body) : size(body.protein_atoms.size()) {
        for (size_t i = 0; i < size; i++) {
            const Atom& a = body.protein_atoms[i]; 
            data[4*i] = a.coords.x;
            data[4*i+1] = a.coords.y;
            data[4*i+2] = a.coords.z;
            data[4*i+3] = a.effective_charge*a.occupancy;
        }
    }

    /**
     * @brief Extract the necessary coordinates and weights from a vector of hydration atoms. 
     */
    CompactCoordinates(const vector<Hetatom>& atoms) : size(atoms.size()) {
        for (size_t i = 0; i < size; i++) {
            const Hetatom& a = atoms[i]; 
            data[4*i] = a.coords.x;
            data[4*i+1] = a.coords.y;
            data[4*i+2] = a.coords.z;
            data[4*i+3] = a.effective_charge*a.occupancy;
        }
    }

    size_t size;
    vector<float> data;
};

/**
 * @brief A simple data container for the hydration histogram. 
 *        This is only defined for consistency. 
 */
class HydrationHistogram {
    public: 
        HydrationHistogram() {}
        HydrationHistogram(const vector<double>& p_hh) : p_hh(p_hh) {}

        vector<double> p_hh;
};

/**
 * @brief A simple data container for the partial histograms. 
 */
class PartialHistogram {
    public:
        PartialHistogram() {}
        PartialHistogram(const vector<double>& p_pp, const vector<double>& p_hp) : p_pp(p_pp), p_hp(p_hp) {}

        vector<double> p_pp;
        vector<double> p_hp;
};

/**
 * @brief We also define the MasterHistogram type, which is identical to a PartialHistogram. 
 *        We do this to make += and -= well-defined operations. 
 */
class MasterHistogram {
    public: 
        MasterHistogram() {}

        /**
         * @brief Create a new Master Histogram. 
         * @param p The current histogram. 
         * @param p_base The constant, unchanging part of the histogram. 
         */
        MasterHistogram(const vector<double>& p_base, const vector<int>& axes) : p(p_base), axes(axes), p_base(p_base) {}

        /**
         * @brief Add a PartialHistogram to the MasterHistogram. 
         */
        MasterHistogram& operator+=(const PartialHistogram& rhs) {
            std::transform(rhs.p_pp.begin(), rhs.p_pp.end(), p.begin(), p.begin(), std::plus<double>());
            std::transform(rhs.p_hp.begin(), rhs.p_hp.end(), p.begin(), p.begin(), std::plus<double>());
        }

        /**
         * @brief Add a Hydrationhistogram to the MasterHistogram. 
         */
        MasterHistogram& operator+=(const HydrationHistogram& rhs) {
            std::transform(rhs.p_hh.begin(), rhs.p_hh.end(), p.begin(), p.begin(), std::plus<double>());
        }

        /**
         * @brief Subtract a PartialHistogram from the MasterHistogram. We have to use a lambda since the standard std::minus would
         *        reverse the order of the entries.
         */
        MasterHistogram& operator-=(const PartialHistogram& rhs) {
            std::transform(rhs.p_pp.begin(), rhs.p_pp.end(), p.begin(), p.begin(), [] (const double& r, const double& l) {return l-r;});
            std::transform(rhs.p_hp.begin(), rhs.p_hp.end(), p.begin(), p.begin(), [] (const double& r, const double& l) {return l-r;});
        }

        /**
         * @brief Subtract a HydrationHistogram from the MasterHistogram. We have to use a lambda since the standard std::minus would
         *        reverse the order of the entries.
         */
        MasterHistogram& operator-=(const HydrationHistogram& rhs) {
            std::transform(rhs.p_hh.begin(), rhs.p_hh.end(), p.begin(), p.begin(), [] (const double& r, const double& l) {return l-r;});
        }

        vector<double> p; // The master histogram itself.
        vector<int> axes; // The axes used for the histogram. 

    private:
        // The base part of the histogram which will never change. This contains all internal distances between atoms in each individual body.
        vector<double> p_base; 
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
        PartialHistogramManager(const vector<Body>& bodies, const vector<Hetatom>& hydration_atoms, std::shared_ptr<StateManager> statemanager) 
            : size(bodies.size()), statemanager(statemanager), coords(size), bodies(bodies), hydration_atoms(hydration_atoms), partials(size) {}

        /**
         * @brief Initialize this object. The internal distances between atoms in each body is constant and cannot change. 
         *        They are unaffected by both rotations and translations, and so we precalculate them. 
         */
        void initialize() {
            // generous sizes - 1000Å should be enough for just about any structure
            double width = setting::axes::scattering_intensity_plot_binned_width;
            vector<int> axes = {int(1000/width), 0, 1000}; 
            vector<double> p_base(axes[0], 0);

            for (int n = 0; n < size; n++) {
                // create more efficient access to the necessary variables
                CompactCoordinates current(bodies[n]);

                // calculate internal distances between atoms
                for (size_t i = 1; i < current.size; i++) {
                    for (size_t j = i; j < current.size; j++) {
                        float weight = current.data[4*i+3]*current.data[4*j+3];
                        float dx = current.data[4*i] - current.data[4*j];
                        float dy = current.data[4*i+1] - current.data[4*j+1];
                        float dz = current.data[4*i+2] - current.data[4*j+2];
                        float dist = sqrt(dx*dx + dy*dy + dz*dz);
                        p_base[dist/width] += 2*weight;
                    }
                }

                // calculate self-correlation
                for (size_t i = 0; i < current.size; i++) {p_base[0] += current.data[4*i+3]*current.data[4*i+3];}

                // store the coordinates for later
                coords[n] = current;
            }
            master = MasterHistogram(p_base, axes);
        }

        /**
         * @brief Calculate the total scattering histogram. 
         */
        ScatteringHistogram calculate() {
            if (statemanager->get_modified_hydration()) {
                calc_hh();

                for (int i = 0; i < size; i++) {
                    calc_hp(i);
                }
            }

            const vector<bool> modified_state = statemanager->get_modified_bodies();
            for (int i = 0; i < size; i++) {
                if (modified_state[i]) {
                    calc_hp(i);
                }
            }
        }

        /**
         * @brief Calculate the atom-atom distances between body @a i and all others. 
         */
        vector<double> calc_pp(const int& i) {

        }

        /**
         * @brief Calculate the hydration-atom distances between the hydration layer and body @a i.
         */
        vector<double> calc_hp(const int& i) {}

        /**
         * @brief Calculate the hydration-hydration distances. 
         */
        vector<double> calc_hh() {
            master -= partial_hydration; // subtract the previous hydration histogram

            double width = setting::axes::scattering_intensity_plot_binned_width;
            const vector<int>& axes = master.axes; 
            vector<double> p_hh(axes[0], 0);

            // calculate internal distances for the hydration layer
            CompactCoordinates coords(hydration_atoms);
            for (size_t i = 1; i < hydration_atoms.size(); i++) {
                for (size_t j = i; j < hydration_atoms.size(); j++) {
                    float weight = coords.data[4*i+3]*coords.data[4*j+3];
                    float dx = coords.data[4*i] - coords.data[4*j];
                    float dy = coords.data[4*i+1] - coords.data[4*j+1];
                    float dz = coords.data[4*i+2] - coords.data[4*j+2];
                    float dist = sqrt(dx*dx + dy*dy + dz*dz);
                    p_hh[dist/width] += 2*weight;
                }
            }

            // calculate self-correlation
            for (size_t i = 0; i < hydration_atoms.size(); i++) {p_hh[0] += coords.data[4*i+3]*coords.data[4*i+3];}

            partial_hydration.p_hh = std::move(p_hh);
            master += partial_hydration; // add the new hydration histogram
        }

    private:
        const size_t size; // number of managed bodies
        std::shared_ptr<StateManager> statemanager; // a helper which keeps track of state changes in each body
        vector<CompactCoordinates> coords; // a compact representation of the relevant data from the managed bodies
        const vector<Body>& bodies; // reference access to the managed bodies
        const vector<Hetatom>& hydration_atoms; // reference to the hydration layer 

        // histogram data
        MasterHistogram master; // the current total histogram
        vector<PartialHistogram> partials; // the partial histograms
        HydrationHistogram partial_hydration; // the partial histogram for the hydration layer
};