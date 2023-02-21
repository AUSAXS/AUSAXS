#pragma once

#include <algorithm>

#include <data/Protein.h>
#include <rigidbody/Constraint.h>

struct BodySplitter {
    /**
     * @brief Constructor. 
     * 
     * @param input The path to the input data file. 
     */
    static Protein split(const std::string& input, std::vector<int> splits);

    /**
     * @brief Creates a single constraint between each body and the next in the sequence. 
     * 
     * @param protein The protein to be constrained. 
     */
    static std::vector<rigidbody::Constraint> sequential_constraints(const Protein& protein);
};