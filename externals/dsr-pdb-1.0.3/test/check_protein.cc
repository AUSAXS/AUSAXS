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

#include "dsrpdb/Protein.h"
#include "dsrpdb/geometry.h"
#include "dsrpdb/iterator.h"
#include <fstream>
#include <cassert>
#include <iterator>

int main(int argc, char *argv[]){
  //dsr::Residue res= dsr::Residue(dsr::Residue::VAL);
  //res.write(std::cout);
  assert(argc==3);
  std::ifstream in(argv[1]);
  dsrpdb::Protein p(in);
  //p.write(std::cout);
  std::ofstream of(argv[2]);
  
  std::cout << "There are " << p.number_of_residues() << " residues." << std::endl;
  
  p.write_pdb(of);
  //p.dump(std::cout);

  // Demonstrate the geometry functions
  std::vector<dsrpdb::Point> atms;
  dsrpdb::coordinates(p.atoms_begin(), p.atoms_end(), std::back_inserter(atms));
  std::cout << "There are " << atms.size() << " atoms." << std::endl;
  assert(std::distance(p.atoms_begin(), p.atoms_end()) == atms.size());
  
  std::vector<dsrpdb::Protein::Bond> bds(p.bonds_begin(), p.bonds_end());
  std::cout << "There are " << bds.size() << " bonds." << std::endl;
  

  std::vector<dsrpdb::Point> bb;
  dsrpdb::backbone_coordinates(p.atoms_begin(), p.atoms_end(), std::back_inserter(bb));
  assert(bb.size() < atms.size());
  {
    std::vector<std::pair<int,int> > bonds;
    std::vector<dsrpdb::Point> points;
    dsrpdb::coordinates_and_bonds(p, std::back_inserter(points), std::back_inserter(bonds));
    assert(bonds.size() == bds.size());
    assert(points.size() == atms.size());
  }
  {
    std::vector<dsrpdb::Point> points;
    dsrpdb::coordinates(p.atoms_begin(), p.atoms_end(), std::back_inserter(points));
    assert(points.size() == atms.size());
  }
  {
    std::vector<dsrpdb::Point> points;
    std::copy(backbone_coordinates_begin(p),
	      backbone_coordinates_end(p),
	      std::back_inserter(points));
    assert(points.size() == bb.size());
  }
  {
    std::vector<std::pair<int,int> > bonds;
    std::vector<dsrpdb::Point> points;
    dsrpdb::simplified_coordinates_and_bonds(p, 
					     std::back_inserter(points), 
					     std::back_inserter(bonds));
    //std::cerr << bonds.size()  << " " << atms.size() << " " << points.size() << std::endl;
    assert(bonds.size() == points.size()-1);
    assert(points.size() == 4*p.number_of_residues());
  }

  {
    for (dsrpdb::Protein::Atoms_iterator it= p.atoms_begin(); it != p.atoms_end(); ++it){
      dsrpdb::Atom::Index ind= it->second.index();
      dsrpdb::Residue::Atom_label al= it->first;
      dsrpdb::Residue& res= p.residue_containing_atom(ind);
      assert(&res >= &*p.residues_begin() 
	     && &res < &*p.residues_begin() + p.number_of_residues());
      assert(res.atom_label(ind) != dsrpdb::Residue::AL_INVALID);
      assert(res.atom_label(ind) == al);
    }
  }

  {
    for (dsrpdb::Protein::Bonds_iterator it= p.bonds_begin(); it != p.bonds_end(); ++it){
      dsrpdb::Atom::Index ind0= it->first;
      dsrpdb::Atom::Index ind1= it->second;
      assert(p.atom(ind0) != p.atom(dsrpdb::Atom::Index()));
      assert(p.atom(ind1) != p.atom(dsrpdb::Atom::Index()));
      int ir0= p.residue_containing_atom(ind0).index();
      int ir1= p.residue_containing_atom(ind1).index();
      int diff = ir1 - ir0;
      assert(diff==0 || diff ==1);
    }
  }

  return EXIT_SUCCESS;
}
