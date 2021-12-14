#include <vector>

using std::vector;

/**
 * @brief A state manager which keeps track of changes in each body. 
 *        This is meant to be used in conjunction with DistanceCalculator, such that it only recalculates what is absolutely necessary. 
 */
class StateManager {
public:
    StateManager(const int& size) : size(size), _modified_protein_atoms(size, false), _modified_hydration(false) {}

    /**
     * @brief Mark that the protein atoms of a body was modified.
     * @param i index of the body. 
     */
    void modified_protein_atoms(const int& i) {
        _modified_protein_atoms[i] = true;
    }

    /**
     * @brief Mark that the hydration atoms of a body was modified.
     * @param i index of the body. 
     */
    void modified_hydration_layer() {
        _modified_hydration = true;
    }

    /**
     * @brief Reset all marks to false.
     */
    void reset() {
        _modified_protein_atoms = vector<bool>(size, false);
        _modified_hydration = false;
    }

private:
    const int size;
    vector<bool> _modified_protein_atoms;
    bool _modified_hydration;
};