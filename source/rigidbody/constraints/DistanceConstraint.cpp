/*
This software is distributed under the GNU Lesser General Public License v3.0. 
For more information, please refer to the LICENSE file in the project root.
*/

#include <rigidbody/constraints/DistanceConstraint.h>
#include <settings/RigidBodySettings.h>
#include <constants/Constants.h>
#include <data/Molecule.h>
#include <data/Body.h>

using namespace ausaxs;
using namespace ausaxs::rigidbody::constraints;
using namespace ausaxs::data;

DistanceConstraint::DistanceConstraint(data::Molecule* protein, unsigned int ibody1, unsigned int ibody2, unsigned int iatom1, unsigned int iatom2) 
    : protein(protein), ibody1(ibody1), ibody2(ibody2), iatom1(iatom1), iatom2(iatom2) {
    const Body& body1 = protein->get_body(ibody1);
    const Body& body2 = protein->get_body(ibody2);
    const AtomFF& atom1 = body1.get_atom(iatom1);
    const AtomFF& atom2 = body2.get_atom(iatom2);

    // we only want to allow constraints between the backbone C-alpha structure
    if (form_factor::to_atom_type(atom1.form_factor_type()) != constants::atom_t::C || form_factor::to_atom_type(atom2.form_factor_type()) != constants::atom_t::C) {
        throw except::invalid_argument("DistanceConstraint::DistanceConstraint: Constraints only makes sense between the carbon-atoms of the backbone!");
    }

    // constraints within the same body doesn't make sense
    if (body1.get_uid() == body2.get_uid()) {
        throw except::invalid_argument("DistanceConstraint::DistanceConstraint: Cannot create a constraint between atoms in the same body!");
    }

    // set the base radius and perform a sanity check
    r_base = atom1.coordinates().distance(atom2.coordinates());
    if (r_base > 4) {
        throw except::invalid_argument(
            "DistanceConstraint::DistanceConstraint: The atoms being constrained are too far apart!\n"
            "Atom 1: " + form_factor::to_string(atom1.form_factor_type()) + " in body " + std::to_string(ibody1) + " at " + atom1.coordinates().to_string() + "\n"
            "Atom 2: " + form_factor::to_string(atom2.form_factor_type()) + " in body " + std::to_string(ibody2) + " at " + atom2.coordinates().to_string() + "\n"
        );
    }
}

DistanceConstraint::DistanceConstraint(data::Molecule* protein, unsigned int ibody1, unsigned int ibody2, bool center_mass)
    : protein(protein), ibody1(ibody1), ibody2(ibody2) {
    const Body& body1 = protein->get_body(ibody1);
    const Body& body2 = protein->get_body(ibody2);

    if (center_mass) {
        // find the atoms closest to the center of mass of the two bodies
        double min_distance = std::numeric_limits<double>::max();
        double target_distance = body1.get_cm().distance2(body2.get_cm());
        for (unsigned int i = 0; i < body1.get_atoms().size(); i++) {
            if (form_factor::to_atom_type(body1.get_atom(i).form_factor_type()) != constants::atom_t::C) {
                continue;
            }

            for (unsigned int j = 0; j < body2.get_atoms().size(); j++) {
                if (form_factor::to_atom_type(body2.get_atom(j).form_factor_type()) != constants::atom_t::C) {
                    continue;
                }

                double distance = body1.get_atom(i).coordinates().distance2(body2.get_atom(j).coordinates());
                if (std::abs(distance - target_distance) < min_distance) {
                    min_distance = std::abs(distance - target_distance);
                    iatom1 = i;
                    iatom2 = j;
                }
            }
        }
        r_base = min_distance;
    } else {
        // find the closest atoms in the two bodies
        double min_distance = std::numeric_limits<double>::max();
        for (unsigned int i = 0; i < body1.get_atoms().size(); i++) {
            if (form_factor::to_atom_type(body1.get_atom(i).form_factor_type()) != constants::atom_t::C) {
                continue;
            }

            for (unsigned int j = 0; j < body2.get_atoms().size(); j++) {
                if (form_factor::to_atom_type(body2.get_atom(j).form_factor_type()) != constants::atom_t::C) {
                    continue;
                }

                double distance = body1.get_atom(i).coordinates().distance2(body2.get_atom(j).coordinates());
                if (distance < min_distance) {
                    min_distance = distance;
                    iatom1 = i;
                    iatom2 = j;
                }
            }
        }
        r_base = min_distance;
        if (r_base > 4) {
            auto& atom1 = body1.get_atom(iatom1);
            auto& atom2 = body2.get_atom(iatom2);
            throw except::invalid_argument(
                "DistanceConstraint::DistanceConstraint: The atoms being constrained are too far apart!\n"
                "Atom 1: " + form_factor::to_string(atom1.form_factor_type()) + " in body " + std::to_string(ibody1) + " at " + atom1.coordinates().to_string() + "\n"
                "Atom 2: " + form_factor::to_string(atom2.form_factor_type()) + " in body " + std::to_string(ibody2) + " at " + atom2.coordinates().to_string() + "\n"
            );
        }
    }
    if (iatom1 == iatom2) {throw except::invalid_argument("DistanceConstraint::DistanceConstraint: Could not find atoms to constrain.");}
}

DistanceConstraint::DistanceConstraint(data::Molecule* protein, const AtomFF& atom1, const AtomFF& atom2) : protein(protein) {
    // we only want to allow constraints between the backbone C-alpha structure
    if (form_factor::to_atom_type(atom1.form_factor_type()) != constants::atom_t::C || form_factor::to_atom_type(atom2.form_factor_type()) != constants::atom_t::C) {
        throw except::invalid_argument("DistanceConstraint::DistanceConstraint: Constraints only makes sense between the carbon-atoms of the backbone!");
    }

    auto[loc1, loc2] = find_host_bodies(atom1, atom2);

    ibody1 = loc1.body;
    ibody2 = loc2.body;
    iatom1 = loc1.atom;
    iatom2 = loc2.atom;

    // constraints within the same body doesn't make sense
    if (ibody1 == ibody2) {throw except::invalid_argument("DistanceConstraint::DistanceConstraint: Cannot create a constraint between atoms in the same body!");}

    // set the base radius and perform a sanity check
    r_base = atom1.coordinates().distance(atom2.coordinates());
    if (r_base > settings::rigidbody::bond_distance) {
        throw except::invalid_argument("DistanceConstraint::DistanceConstraint: The atoms being constrained are too far apart!");
    }
}

std::pair<DistanceConstraint::AtomLoc, DistanceConstraint::AtomLoc> DistanceConstraint::find_host_bodies(const AtomFF& atom1, const AtomFF& atom2) const {
    int ibody1 = -1, ibody2 = -1;
    int iatom1 = -1, iatom2 = -1;
    for (unsigned int ibody = 0; ibody < protein->size_body(); ++ibody) {
        const Body& body = protein->get_body(ibody);
        for (unsigned int iatom = 0; iatom < body.size_atom(); ++iatom) {
            auto& atom = body.get_atom(iatom);
            if (atom1 == atom) {
                ibody1 = ibody;
                iatom1 = iatom;
                break; // a1 and a2 *must* be from different bodies, so we break
            } else if (atom2 == atom) {
                ibody2 = ibody;
                iatom2 = iatom;
                break; // same
            }
        }
    }

    // check that both b1 and b2 were found
    if (ibody1 == -1 || ibody2 == -1) {
        throw except::invalid_argument("DistanceConstraint::find_host_bodies: Could not determine host bodies for the two atoms.");
    }

    return std::make_pair(AtomLoc{.body=ibody1, .atom=iatom1}, AtomLoc{.body=ibody2, .atom=iatom2});
}

double DistanceConstraint::evaluate() const {
    const AtomFF& atom1 = protein->get_body(ibody1).get_atom(iatom1);
    const AtomFF& atom2 = protein->get_body(ibody2).get_atom(iatom2);
    return transform(r_base - atom1.coordinates().distance(atom2.coordinates()));
}

bool DistanceConstraint::operator==(const DistanceConstraint& constraint) const {
    return ibody1 == constraint.ibody1 
        && ibody2 == constraint.ibody2 
        && iatom1 == constraint.iatom1 
        && iatom2 == constraint.iatom2;
}

double DistanceConstraint::transform(double offset) {
    return offset*offset*offset*offset*10;
}

const AtomFF& DistanceConstraint::get_atom1() const {
    return protein->get_body(ibody1).get_atom(iatom1);
}

const AtomFF& DistanceConstraint::get_atom2() const {
    return protein->get_body(ibody2).get_atom(iatom2);
}

const Body& DistanceConstraint::get_body1() const {
    return protein->get_body(ibody1);
}

const Body& DistanceConstraint::get_body2() const {
    return protein->get_body(ibody2);
}

Body& DistanceConstraint::get_body1() {
    return protein->get_body(ibody1);
}

Body& DistanceConstraint::get_body2() {
    return protein->get_body(ibody2);
}

std::string DistanceConstraint::to_string() const {
    std::stringstream ss;
    ss << "Constraint between (" << form_factor::to_string(protein->get_body(ibody1).get_atom(iatom1).form_factor_type()) << ") and "
                             "(" << form_factor::to_string(protein->get_body(ibody2).get_atom(iatom2).form_factor_type()) << ")";
    return ss.str();
}