// includes
#include <vector>
#include <string>

// my own includes
#include "source/pdbml_reader.h"
#include "source/Atom.cpp"
#include "source/Protein.cpp"
#include "source/tools.cpp"

int main(int argc, char const *argv[])
{
    Protein protein("temp.xml");
    auto[dp, dh] = protein.calc_distances();
    save(dp, "temp/distances_protein.txt");
    save(dh, "temp/distances_hydration.txt");

    std::cout << "Protein atoms:" << std::endl;
    for (Atom* a : protein.protein) {
        a->print();
    }

    std::cout << "Hydration atoms:" << std::endl;
    for (Atom* a : protein.hydration) {
        a->print();
    }
    return 0;
}