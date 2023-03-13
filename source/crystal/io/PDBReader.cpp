#include <crystal/io/PDBReader.h>

#include <data/Protein.h>
#include <utility/Axis3D.h>

std::pair<Basis3D, std::vector<Vector3<double>>> crystal::io::PDBReader::read(const std::string& input) const {
    Protein protein(input);
    auto prot_atoms = protein.atoms();

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
    axis.x.min -= axis.x.span()*5;
    axis.x.max += axis.x.span()*5;
    axis.y.min -= axis.y.span()*5;
    axis.y.max += axis.y.span()*5;
    axis.z.min -= axis.z.span()*5;
    axis.z.max += axis.z.span()*5;

    protein.translate({-axis.x.min, -axis.y.min, -axis.z.min});
    std::vector<Vector3<double>> atoms;
    for (const auto& atom : prot_atoms) {
        atoms.push_back(atom.get_coordinates());
    }

    Basis3D basis({axis.x.span(), 0, 0}, {0, axis.y.span(), 0}, {0, 0, axis.z.span()});
    return std::make_pair(std::move(basis), std::move(atoms));
}