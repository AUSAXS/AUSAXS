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
    CompactCoordinates(const Body& body) : size(body.protein_atoms.size()) {
        for (size_t i = 0; i < body.protein_atoms.size(); i++) {
            const Atom& a = body.protein_atoms[i]; 
            data[4*i] = a.coords.x;
            data[4*i+1] = a.coords.y;
            data[4*i+2] = a.coords.z;
            data[4*i+3] = a.effective_charge*a.occupancy;
        }
    }

    const size_t size;
    vector<float> data;
};

/**
 * @brief A smart distance calculator which efficiently calculates the scattering histogram.
 */
class DistanceCalculator {
    public:
        DistanceCalculator(const vector<Body>& bodies, const vector<Hetatom>& hydration_atoms, std::shared_ptr<StateManager> statemanager) 
            : statemanager(statemanager), coords(bodies.size()), bodies(bodies), hydration_atoms(hydration_atoms) {}

        /**
         * @brief Initialize this object. The internal distances between atoms in each body is constant and cannot change. 
         *        They are unaffected by both rotations and translations, and so we precalculate them. 
         */
        void initialize() {
            for (int i = 0; i < bodies.size(); i++) {
                coords[i] = CompactCoordinates(bodies[i]);
            }

            // generous sizes - 1000Å should be enough for just about any structure
            double width = setting::axes::scattering_intensity_plot_binned_width;
            vector<int> axes = {int(1000/width), 0, 1000}; 
            vector<double> p_pp(axes[0], 0);
            vector<double> p_hh(axes[0], 0);
            vector<double> p_hp(axes[0], 0);
            vector<double> p_tot(axes[0], 0);

            for (const auto& body : bodies) {

            }
        }

        ScatteringHistogram calculate() {
            if (statemanager->get_modified_hydration()) {
                calc_hh();

                for (int i = 0; i < bodies.size(); i++) {
                    calc_hp(i);
                }
            }

            const vector<bool> modified_state = statemanager->get_modified_bodies();
            for (int i = 0; i < bodies.size(); i++) {
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
        vector<double> calc_hh() {}

    private:
        std::shared_ptr<StateManager> statemanager;
        vector<CompactCoordinates> coords; 
        const vector<Body>& bodies;
        const vector<Hetatom>& hydration_atoms;
};