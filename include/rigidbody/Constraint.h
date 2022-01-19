#pragma once

#include <memory>

#include "data/Atom.h"
#include "Exceptions.h"

/**
 * @brief \class Constraint. 
 * 
 * This class is the glue that keeps separate bodies together during the optimization. Each constraint is between two individual atoms of two different bodies, and works by 
 * adding a new term to the chi-square which is a function of the distance between the two atoms. Thus the optimizer is penalized depending on how much it separates the two 
 * constrained atoms. 
 */
class Constraint {
  public: 
    /**
     * @brief Constructor. 
     * 
     * Create a new constraint between a pair of atoms. 
     * 
     * @param atom1 
     * @param atom2 
     * @param body1 
     * @param body2 
     */
    Constraint(const Atom* const atom1, const Atom* const atom2, const Body* const body1, const Body* const body2) 
        : atom1(atom1), atom2(atom2), body1(body1), body2(body2) {

        // we only want to allow constraints between the backbone C-alpha structure
        if (atom1->element != "C" || atom2->element != "C") {
            throw except::invalid_argument("Error in Constraint::Constraint: Constraints only makes sense between the carbon-atoms of the backbone!");
        }

        // constraints within the same body doesn't make sense
        if (body1->uid == body2->uid) {
            throw except::invalid_argument("Error in Constraint::Constraint: Cannot create a constraint between atoms in the same body!");
        }

        // set the base radius and perform a sanity check
        r_base = atom1->distance(*atom2);
        if (r_base > 4) {
            throw except::invalid_argument("Error in Constraint::Constraint: The atoms being constrained are too far apart!");
        }

        // as long as uid_counter does not change after object creation, this should be unique
        uid = atom1->uid*Atom::uid_counter + atom2->uid;
    }

    /**
     * @brief Evaluate this constraint for the current positions. 
     */
    double evaluate() {return transform(r_base - atom1->distance(*atom2));}

    /**
     * @brief Check if a constraint is identical to this object. 
     * 
     * @param constraint The constraint to be checked for equality. 
     */
    bool operator==(const Constraint& constraint) const {
        return atom1 == constraint.atom1 && atom2 == constraint.atom2;
    }

    size_t uid;              // Unique identifier for this constraint. 
    double r_base;           // The normal distance between the two atoms. 
    const Atom* const atom1; // The first atom. 
    const Atom* const atom2; // The second atom. 
    const Body* const body1; // The first body.
    const Body* const body2; // The second body.

    /**
     * @brief Transforms a distance into a proper constraint for least-squares fitting. 
     * 
     * @param offset The radial offset between the new and original positions. 
     */
    static double transform(const double offset) {
        return offset*offset*offset*offset;
    }
};