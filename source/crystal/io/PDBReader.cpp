#include <crystal/io/PDBReader.h>

#include <data/Protein.h>
#include <utility/Axis3D.h>

std::pair<Basis3D, std::vector<Vector3<double>>> crystal::io::PDBReader::read(const std::string& input) const {
    Protein protein(input);
    auto prot_atoms = protein.atoms();

    double expansion = 5;

    Axis3D axis;
    for (const auto& atom : prot_atoms) {
        auto position = atom.get_coordinates();
        if (position.x() < axis.x.min) {axis.x.min = position.x();}
        if (position.x() > axis.x.max) {axis.x.max = position.x();}
        if (position.y() < axis.y.min) {axis.y.min = position.y();}
        if (position.y() > axis.y.max) {axis.y.max = position.y();}
        if (position.z() < axis.z.min) {axis.z.min = position.z();}
        if (position.z() > axis.z.max) {axis.z.max = position.z();}
    }

    // center the protein in the middle of the box
    double factor = (expansion-1)/2;
    protein.translate({-axis.x.min + factor*axis.x.span(), -axis.y.min + factor*axis.y.span(), -axis.z.min + factor*axis.z.span()});
    std::vector<Vector3<double>> atoms;
    for (const auto& atom : prot_atoms) {
        atoms.push_back(atom.get_coordinates());
    }

    // expand the box by 2 times
    axis.x.min -= axis.x.span()*expansion;
    axis.x.max += axis.x.span()*expansion;
    axis.y.min -= axis.y.span()*expansion;
    axis.y.max += axis.y.span()*expansion;
    axis.z.min -= axis.z.span()*expansion;
    axis.z.max += axis.z.span()*expansion;

    Basis3D basis({axis.x.span(), 0, 0}, {0, axis.y.span(), 0}, {0, 0, axis.z.span()});
    return std::make_pair(std::move(basis), std::move(atoms));
}