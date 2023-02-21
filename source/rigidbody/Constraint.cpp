#include <rigidbody/Constraint.h>
#include <data/Protein.h>

using namespace rigidbody;

Constraint::Constraint(const Protein* protein, unsigned int ibody1, unsigned int ibody2, unsigned int iatom1, unsigned int iatom2) 
    : protein(protein), ibody1(ibody1), ibody2(ibody2), iatom1(iatom1), iatom2(iatom2) {
    const Body& body1 = protein->body(ibody1);
    const Body& body2 = protein->body(ibody2);
    const Atom& atom1 = body1.atoms(iatom1);
    const Atom& atom2 = body2.atoms(iatom2);

    // we only want to allow constraints between the backbone C-alpha structure
    if (atom1.element != constants::symbols::carbon || atom2.element != constants::symbols::carbon) {
        throw except::invalid_argument("Constraint::Constraint: Constraints only makes sense between the carbon-atoms of the backbone!");
    }

    // constraints within the same body doesn't make sense
    if (body1.uid == body2.uid) {throw except::invalid_argument("Constraint::Constraint: Cannot create a constraint between atoms in the same body!");}

    // set the base radius and perform a sanity check
    r_base = atom1.distance(atom2);
    if (r_base > 4) {throw except::invalid_argument("Constraint::Constraint: The atoms being constrained are too far apart!");}

    uid = uid_counter++;
}

Constraint::Constraint(const Protein* protein, const Atom& atom1, const Atom& atom2) : protein(protein) {
    // we only want to allow constraints between the backbone C-alpha structure
    if (atom1.element != constants::symbols::carbon || atom2.element != constants::symbols::carbon) {
        throw except::invalid_argument("Constraint::Constraint: Constraints only makes sense between the carbon-atoms of the backbone!");
    }

    auto[loc1, loc2] = find_host_bodies(atom1, atom2);

    ibody1 = loc1.body;
    ibody2 = loc2.body;
    iatom1 = loc1.atom;
    iatom2 = loc2.atom;

    // constraints within the same body doesn't make sense
    if (ibody1 == ibody2) {throw except::invalid_argument("Constraint::Constraint: Cannot create a constraint between atoms in the same body!");}

    // set the base radius and perform a sanity check
    r_base = atom1.distance(atom2);
    if (r_base > 4) {throw except::invalid_argument("Constraint::Constraint: The atoms being constrained are too far apart!");}

    uid = uid_counter++;
}

std::pair<rigidbody::detail::AtomLoc, rigidbody::detail::AtomLoc> Constraint::find_host_bodies(const Atom& atom1, const Atom& atom2) const {
    int ibody1 = -1, ibody2 = -1;
    int iatom1 = -1, iatom2 = -1;
    for (unsigned int ibody = 0; ibody < protein->bodies.size(); ibody++) {
        const Body& body = protein->body(ibody);
        for (unsigned int iatom = 0; iatom < body.atoms().size(); iatom++) {
            if (atom1 == body.atoms(iatom)) {
                ibody1 = ibody;
                iatom1 = iatom;
                break; // a1 and a2 *must* be from different bodies, so we break
            } else if (atom2 == body.atoms(iatom)) {
                ibody2 = ibody;
                iatom2 = iatom;
                break; // same
            }
        }
    }

    // check that both b1 and b2 were found
    if (ibody1 == -1 || ibody2 == -1) {
        throw except::invalid_argument("RigidBody::create_constraint: Could not determine host bodies for the two atoms.");
    }

    return std::make_pair(detail::AtomLoc{.body=ibody1, .atom=iatom1}, detail::AtomLoc{.body=ibody2, .atom=iatom2});
}

double Constraint::evaluate() const {
    const Atom& atom1 = protein->body(ibody1).atoms(iatom1);
    const Atom& atom2 = protein->body(ibody2).atoms(iatom2);
    return transform(r_base - atom1.distance(atom2));
}

bool Constraint::operator==(const Constraint& constraint) const {
    return ibody1 == constraint.ibody1 
        && ibody2 == constraint.ibody2 
        && iatom1 == constraint.iatom1 
        && iatom2 == constraint.iatom2;
}

double rigidbody::Constraint::transform(double offset) {
    return offset*offset*offset*offset;
}

const Atom& Constraint::get_atom1() const {
    return protein->body(ibody1).atoms(iatom1);
}

const Atom& Constraint::get_atom2() const {
    return protein->body(ibody2).atoms(iatom2);
}

const Body& Constraint::get_body1() const {
    return protein->body(ibody1);
}

const Body& Constraint::get_body2() const {
    return protein->body(ibody2);
}