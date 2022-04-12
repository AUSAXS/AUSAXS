#pragma once

#include <vector>
#include <memory>

/**
 * @brief A state manager which keeps track of changes in each body. 
 *        This is meant to be used in conjunction with DistanceCalculator, such that it only recalculates what is absolutely necessary. 
 */
class StateManager {
  public:
    /**
     * @brief A small probe for signalling changes which can be dispatched to other classes. 
     */
    class Signaller {
        public: 
            Signaller(unsigned int id, StateManager* const owner);

            /**
             * @brief Signal that the state of this object has changed. 
             */
            virtual void state_change() const;

        private: 
            StateManager* const owner;
            int id;
    };

    /**
     * @brief Dummy version of a Signaller object. This can be used to initialize an instance of Signaller. 
     */
    class UnboundSignaller : public Signaller {
        public: 
            UnboundSignaller();

            /**
             * @brief Does nothing.
             */
            void state_change() const override;
    };

    StateManager(unsigned int size);

    /**
     * @brief Mark that the protein atoms of all bodies were modified. 
     */
    void modified_all();

    /**
     * @brief Mark that the protein atoms of a body was modified.
     * @param i index of the body. 
     */
    void modified(const int i);

    /**
     * @brief Mark that the hydration atoms of a body was modified.
     * @param i index of the body. 
     */
    void modified_hydration_layer();

    /**
     * @brief Reset all marks to false.
     */
    void reset();

    /**
     * @brief Get a pointer to the @a ith probe so it can be dispatched to other classes.
     */
    std::shared_ptr<Signaller> get_probe(unsigned int i);

    /**
     * @brief Get a boolean vector which denotes if the state of a given body was changed. 
     */
    std::vector<bool> get_modified_bodies() const;

    /**
     * @brief Check if a given body has been marked as modified.
     */
    [[nodiscard]] bool is_modified(unsigned int i);

    /**
     * @brief Returns true if the hydration layer has been modified, false otherwise. 
     */
    bool get_modified_hydration() const;

  private:
    const int size;
    std::vector<bool> _modified;
    bool _modified_hydration;
    std::vector<std::shared_ptr<Signaller>> probes;
};