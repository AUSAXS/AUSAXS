#include <dsrpdb/Protein.h>
#include <dsrpdb/geometry.h>
#include <fstream>
#include <vector>
#include <cassert>

int main(int, char *[]){
  std::ifstream input("check_bonds.pdb");
  dsrpdb::Protein helix(input);
  //helix.dump(std::cout);
  std::vector<dsrpdb::Point> points;
  std::vector<dsrpdb::Protein::Bond> bonds;
  dsrpdb::backbone_coordinates_and_bonds(helix, std::back_inserter(points), std::back_inserter(bonds));
  assert(bonds.size() == points.size()-2);
  return EXIT_SUCCESS;
}
