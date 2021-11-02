// includes
#include <vector>
#include <string>

// my own includes
#include "source/pdbml_reader.h"
#include "source/Atom.cpp"

int main(int argc, char const *argv[])
{
    std::string file = "temp.xml";
    pdbml_reader reader(file);
    std::vector<Atom*> atoms = reader.read();
    for (Atom* a : atoms) {
        a->print();
        std::cout << a->to_pdbml() << std::endl;
    }

    return 0;
}