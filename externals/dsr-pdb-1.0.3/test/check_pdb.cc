/* Copyright 2004
Stanford University

This file is part of the DSR PDB Library.

The DSR PDB Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The DSR PDB Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the DSR PDB Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#include "dsrpdb/PDB.h"
#include <fstream>
#include <cassert>

int main(int argc, char *argv[]){
  //dsr::Residue res= dsr::Residue(dsr::Residue::VAL);
  //res.write(std::cout);
  assert(argc==3);
  std::ifstream in(argv[1]);
  dsrpdb::PDB p(in);
  //p.write(std::cout);
  std::ofstream of(argv[2]);
  
  std::cout << "There are " << p.number_of_models() << " models." << std::endl;
  
  for (unsigned int i=0; i< p.number_of_models(); ++i){
    const dsrpdb::Model &m= p.model(i);
    std::cout << "Model " << i << " has " << m.number_of_chains() << " chains" << std::endl;
    for (unsigned int j=0; j< m.number_of_chains(); ++j){
      std::cout << "Chain " << j << " has " << m.chain(j).number_of_residues() << " residues" << std::endl;
    }
  }

  p.write(of);
  
  return EXIT_SUCCESS;
}
