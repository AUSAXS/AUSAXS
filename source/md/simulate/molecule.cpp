#include <md/simulate/molecule.h>

using namespace gmx;
using namespace gmx::simulate;

Molecule::Molecule(MoleculeOptions& options) : options(options) {}

SimulateMoleculeOutput Molecule::simulate() {return SimulateMoleculeOutput();}

std::tuple<GROFile, TOPFile> Molecule::setup() {return std::make_tuple(GROFile(), TOPFile());}

GROFile Molecule::minimize() {return GROFile();}

std::tuple<GROFile, NDXFile> Molecule::thermalize() {return std::make_tuple(GROFile(), NDXFile());}