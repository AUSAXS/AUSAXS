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
#include <dsrpdb/align.h>
#include <dsrpdb/align.h>
#include <dsrpdb/distance.h>
#include <dsrpdb/Transform.h>

int main(int argc, char *argv[]){
  using namespace dsrpdb;
  //dsr::Residue res= dsr::Residue(dsr::Residue::VAL);
  //res.write(std::cout);
  assert(argc==2);
  std::ifstream in(argv[1]);
  dsrpdb::PDB p(in);
  
  std::cout << "There are " << p.number_of_models() << " models." << std::endl;
  
  for (unsigned int i=0; i< p.number_of_models(); ++i){
    const dsrpdb::Model &m= p.model(i);
    std::cout << "Model " << i << " has " << m.number_of_chains() << " chains" << std::endl;
    for (unsigned int j=0; j< m.number_of_chains(); ++j){
      dsrpdb::Protein p = m.chain(j);
      double rot[3][3]={{0.309976, 0.851651, -0.422618},
			{-0.741526, 0.494764, 0.453154},
			{0.595025, 0.172916, 0.78488}},
	trans[3]={4,5,6};
      dsrpdb::Transform tr(rot, trans);

      for (dsrpdb::Protein::Residues_iterator rit= p.residues_begin();
	   rit != p.residues_end(); ++rit){
	for (dsrpdb::Residue::Atoms_iterator ait= rit->atoms_begin(); 
	     ait != rit->atoms_end(); ++ait){
	  ait->second.set_cartesian_coords(tr(ait->second.cartesian_coords()));
	}
      }

      std::vector<dsrpdb::Point> pa, pp;
      dsrpdb::backbone_coordinates(m.chain(j).atoms_begin(), m.chain(j).atoms_end(),
				   std::back_inserter(pa));
      dsrpdb::backbone_coordinates(p.atoms_begin(), p.atoms_end(),
				   std::back_inserter(pp));

      dsrpdb::Transform trp= dsrpdb::compute_transform_taking_first_to_second( pa, pp);

      double mdiff=tr.error(trp);
      assert(mdiff<.1);

      dsrpdb::align_second_protein_to_first(  m.chain(j), p);

      double err= no_align_cRMS(m.chain(j), p);
      double ca_err= cRMS(m.chain(j), p);
      assert(err < 1e-5);
      std::cout << err << " " << ca_err << std::endl;
    }
  }

  
  return EXIT_SUCCESS;
}
