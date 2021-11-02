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
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <dsrpdb/Transform.h>

#ifdef PDB_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif


int main(int argc, char *argv[]){
  std::string input_file, output_file;
  bool print_help=false;
  bool verbose=false;
  bool dali=false;
  bool matrix=false;

#ifdef PDB_USE_BOOST_PROGRAM_OPTIONS
  boost::program_options::options_description o("Allowed options"), po, ao;
  o.add_options()
    ("help", boost::program_options::bool_switch(&print_help),
     "produce help message")
    ("dali,d", boost::program_options::bool_switch(&dali), 
     "Use a DALI transformation matrix as pasted from an email from DALI.")
    ("matrix,m", boost::program_options::bool_switch(&matrix),
     "Enter a single transformation matrix.")
    ("verbose,v", boost::program_options::bool_switch(&verbose),
     "print out any errors that occur during reading of the pdb file.");
   
  po.add_options()
    ("input-pdb", boost::program_options::value< std::string>(&input_file),
     "input file")
    ("output-pdb", boost::program_options::value< std::string>(&output_file),
     "output file");
  ao.add(o).add(po);

  boost::program_options::positional_options_description p;
  p.add("input-pdb", 1);
  p.add("output-pdb", 1);

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				options(ao).positional(p).run(), vm);
  boost::program_options::notify(vm); 

  if (input_file.empty() || output_file.empty() || print_help || (dali && matrix)) {
    std::cout << "This program transforms a pdb file by reading a transformation"
	      << " matrix. The matrix is either specified by pasting the transform lines"
	      << " from a DALI email, entering a 4x4 transformation matrix (with a"
	      << " coordinate ordering x,y,z,w) or as separate rotational and "
	      << " translational components.\n";
    std::cout << "usage: " << argv[0] << " input-pdb output-pdb\n\n";
    std::cout << o << "\n";
    return EXIT_FAILURE;
  }

#else
  if (argc != 3){
    std::cerr << "The program is built without boost/program_options.hpp.\n";
    std::cerr << "useage: " << argv[0] << " [-c] file.pdb output.pdb" << std::endl;
    return EXIT_FAILURE;
  }
  input_file = argv[1];
  output_file = argv[2];
#endif

  
  // std::cout << input_file << " " << output_template << " " << split_domain << " " << split_chains << std::endl;

  

  //= new char[strlen(argv[2]+1000)];
 
  double rot[3][3];
  double trans[3];
  if (dali) {
    // const char format[]="%d: %d %s %s %le %le %le %le"
    const char format[]="%le %le %le %le";
    for (unsigned int i=0; i< 3; ++i){
      char buf[1000];
      std::cin.getline(buf, 1000);
      if (sscanf(buf+27, format, &rot[i][0], &rot[i][1], &rot[i][2], &trans[i]) != 4) {
	std::cerr << "Error parsing Dali matrix.\n";
	return EXIT_FAILURE;
      };
      //std::cin >> n >> jc >> js >> js >> x >> y >> z >> t;
      //std::cout <<  x << " " << y << " " << z << " " << t << std::endl;
    }
  } else if (matrix) {
    const char format[]="%le %le %le %le";
    for (unsigned int i=0; i< 4; ++i){
      char buf[1000];
      std::cin.getline(buf, 1000);
      if (sscanf(buf, format, &rot[i][0], &rot[i][1], &rot[i][2], &trans[i]) != 4) {
	std::cerr << "Error parsing matrix, expected 4 floats.\n";
	return EXIT_FAILURE;
      };
    }
  } else {
    std::cout << "Rotation:\n";
    std::cout << "> " << std::flush;
    for (unsigned int i=0; i< 3; ++i){
      char buf[1000];
      std::cin.getline(buf, 1000);
      if (sscanf(buf, "%le %le %le", &rot[i][0], &rot[i][1], &rot[i][2]) != 3) {
	std::cerr << "Error parsing rotation matrix, expected 3 floats.\n";
	return EXIT_FAILURE;
      };
      if (i != 2) {
	std::cout << "> " << std::flush;
      }
    }
    std::cout << "Translation:\n";
    std::cout << "> "<< std::flush;
    char buf[1000];
    std::cin.getline(buf, 1000);
    if (sscanf(buf, "%le %le %le", &trans[0], &trans[1], &trans[2]) != 3) {
      std::cerr << "Error parsing translation, expected 3 floats.\n";
      return EXIT_FAILURE;
    };
  }
  
  dsrpdb::Transform t(rot, trans);
  std::cout  << t;
  std::ifstream in(input_file.c_str());
  if (!in) {
    std::cerr << "Error opening input file " << input_file << std::endl;
    return EXIT_FAILURE;
  }

  dsrpdb::PDB pdb(in, verbose);

  if (verbose) std::cout << "Input PDB has " << pdb.number_of_models() << " models." << std::endl;
 
  for (unsigned int i=0;i< pdb.number_of_models(); ++i){
    dsrpdb::Model &m= pdb.model(i);
    std::cout << "Model " << i << " has " << m.number_of_chains() << " chains."<< std::endl;
    for (unsigned int j=0; j< m.number_of_chains(); ++j){
      dsrpdb::Protein &p= m.chain(j);
      for (dsrpdb::Protein::Residues_iterator rit= p.residues_begin(); rit != p.residues_end(); ++rit){
	for (dsrpdb::Residue::Atoms_iterator ait= rit->atoms_begin(); ait != rit->atoms_end(); ++ait){
	  ait->second.set_cartesian_coords(t(ait->second.cartesian_coords()));
	}
      }
    }
  }

  std::ofstream out(output_file.c_str());
  if (!out) {
    std::cerr << "Error opening output file " << output_file << std::endl;
    return EXIT_FAILURE;
  }

  pdb.write(out);

  //delete[] buf;
  return EXIT_SUCCESS;
}
