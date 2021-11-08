// includes
#include <vector>
#include <string>

// my own includes
#include "Protein.cpp"

int main(int argc, char const *argv[])
{
    if (argc != 3) {
        std::cout << "Missing input!" << std::endl;
        exit(1);
    }
    Protein protein(argv[1]);
    protein.generate_new_hydration();
    protein.save(argv[2]);
    return 0;
}