#pragma once

#include "data/Protein.h"

/**
 * @brief \class TransformationStrategy. 
 * 
 * This super-class defines the interface for the body transformation strategies for the rigid-body optimization. 
 * More specifically its implementations essentially specifies how other connected bodies are affected by a transformation. 
 */
class TransformationStrategy {
    public:
        /**
         * @brief Construtor. 
         */
        TransformationStrategy(const Protein& protein) : protein(protein) {}

        /**
         * @brief Destructor.
         */
        virtual ~TransformationStrategy() = default;

        /**
         * @brief Rotate a body. 
         */
        virtual void rotate(Body& body) = 0;

        /**
         * @brief Translate a body. 
         */
        virtual void translate(Body& body) = 0;

    protected: 
        const Protein& protein; // A reference to the protein to be optimized. We need this to access its constituent bodies. 
};