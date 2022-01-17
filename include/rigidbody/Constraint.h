#pragma once

#include <memory>

#include "data/Atom.h"
#include "Exceptions.h"

/**
 * @brief \class Constraint. 
 * 
 * This class acts like a constraint on the distance between a pair of atoms. 
 */
class Constraint {
  public: 
    /**
     * @brief Constructor.
     * 
     * Create a new constraint between a pair of atoms. The 
     * 
     * @param atom1 
     * @param atom2 
     */
    Constraint(const std::shared_ptr<Atom> const atom1, const std::shared_ptr<Atom> const atom2) : atom1(atom1), atom2(atom2) {
        // we only want to allow constraints between the backbone C-alpha structure
        if (atom1->element != "C" || atom2->element != "C") {
            throw except::invalid_argument("Error in Constraint::Constraint: Constraints only makes sense between the carbon-atoms of the backbone!");
        }

        // set the base radius and perform a sanity check
        r_base = atom1->distance(*atom2);
        if (r_base > 4) {
            throw except::invalid_argument("Error in Constraint::Constraint: The atoms being constrained are too far apart!");
        }
    }

    /**
     * @brief Evaluate this constraint for the current positions. 
     */
    double evaluate() {return transform(r_base - atom1->distance(*atom2));}

  private:
    double r_base;                           // The normal distance between the two atoms. 
    const std::shared_ptr<Atom> const atom1; // The first atom. 
    const std::shared_ptr<Atom> const atom2; // The second atom. 

    /**
     * @brief Transforms a distance into a proper constraint for least-squares fitting. 
     * 
     * @param offset The radial offset between the new and original positions. 
     */
    static double transform(const double offset) {
        return offset*offset*offset*offset;
    }
};