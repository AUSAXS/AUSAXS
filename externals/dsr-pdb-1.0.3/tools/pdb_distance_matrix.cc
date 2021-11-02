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
#include <dsrpdb/distance.h>

#ifdef PDB_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif

#ifdef PDB_USE_MAGICK
#include <Magick++.h>


/*!
  \example pdb_write_distmat.cc

  This example shows how to write a distance matrix to a file.
*/

int main(int argc, char *argv[]){
  std::string input_file, output_file;
  bool print_help=false;
  bool verbose=false;
  bool interactive=false;
  bool all_atoms=false;
  double contact_map_threshold=std::numeric_limits<double>::infinity();

#ifdef PDB_USE_BOOST_PROGRAM_OPTIONS
  {
    boost::program_options::options_description o("Allowed options"), po, ao;
    o.add_options()
      ("help", boost::program_options::bool_switch(&print_help), "produce help message")
      ("verbose,v", boost::program_options::bool_switch(&verbose), "print verbose error messages")
      ("all-atoms,a", boost::program_options::bool_switch(&all_atoms),
       "Output the distances between all atoms, not just the C_alphas.")
      ("image,i",boost::program_options::value<std::string>(&output_file), "Output image file.")
      ("contact-threshold,c",boost::program_options::value<double>(&contact_map_threshold), "Output a contact map with the given threshold. (assuming the image argument is passed).")
      ("query,q", boost::program_options::bool_switch(&interactive), "Enter a mode to do interactive queries about edge lengths.");
    po.add_options()
      ("input-pdb", boost::program_options::value< std::string>(&input_file), "input file");

    ao.add(o).add(po);
    
    boost::program_options::positional_options_description p;
    p.add("input-pdb", 1);
    p.add("output-image", 1);
    
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				  options(ao).positional(p).run(), vm);
    boost::program_options::notify(vm);
    
    
    
    //boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    //boost::program_options::notify(vm);    
    
    if (input_file.empty() || print_help) {
      std::cout << "This program writes the distance matrix to an image file.\n";
      std::cout << "usage: " << argv[0] << " input-pdb\n\n";
      std::cout << o << "\n";
      return EXIT_FAILURE;
    }
  }

#else
  if (argc != 3){
    std::cerr << "The program is built without boost/program_options.hpp.\n";
    std::cerr << "useage: " << argv[0] << " [-c] file.pdb template%d[%c].pdb" << std::endl;
    return EXIT_FAILURE;
  }
  input_file = argv[1];
  output_file = argv[2];
#endif
 
  
  
  
  std::ifstream in(input_file.c_str());
  if (!in) {
    std::cerr << "Error opening input file " << input_file << std::endl;
    return EXIT_FAILURE;
  }

  dsrpdb::PDB pdb(in, verbose);

  if (verbose) std::cout << "Input PDB has " << pdb.number_of_models() << " models." << std::endl;
  if (pdb.number_of_models() != 1){
    std::cout << "Attempting to write multiple image files: assuming the output argument has a %d in it.\n";
  }
 
  for (unsigned int i=0;i< pdb.number_of_models(); ++i){
    dsrpdb::Model &m= pdb.model(i);
    dsrpdb::Matrix arr;
    if (all_atoms) {
      arr= distance_matrix(m);
    } else {
      arr= ca_distance_matrix(m);
    }

    //std::cout << arr.dim1() << " " << arr.dim2() << std::endl;

    if (interactive) {
      std::cout << "Entering interactive mode for model " << i << std::endl;
      std::cout << "There are " << arr.dim1() << " atoms." << std::endl;
      while (true) {
	std::cout << "> " << std::flush;
	char buf[100];
	std::cin.getline(buf, 100);
	if (strcmp(buf, "quit") ==0  || strcmp(buf,"exit") ==0){
	  break;
	} else {
	  int a,b;
	  if (sscanf(buf, "%d %d", &a, &b)==2){
	    if (a > arr.dim1() || b > arr.dim1() || a<0 || b <0){
	      std::cerr << "Index " << a << " or " << b << " out of range." << std::endl;	  
	    } else {
	      std::cout << arr[a][b] << std::endl;
	    }
	  } else if (sscanf(buf, "r %d %d", &a, &b) ==2) {
	    if (a > 0 && b > 0) {
	      dsrpdb::Residue::Index ia(a), ib(b);
	      if (m.chain(0).has_residue(ia) && m.chain(0).has_residue(ib)) {
		dsrpdb::Point pa= m.chain(0).residue(ia).atom(dsrpdb::Residue::AL_CA).cartesian_coords();
		dsrpdb::Point pb= m.chain(0).residue(ib).atom(dsrpdb::Residue::AL_CA).cartesian_coords();
		dsrpdb::Squared_distance sd;
		double d= std::sqrt(sd(pa,pb));
		std::cout << "The Ca distance is " << d << std::endl;
	      } else {
		std::cerr << "Invalid residues picked.\n";
	      }
	    } else {
	      std::cerr << "Invalid residues indices picked (must be greater than 0).\n";
	    }
	  } else {
	    std::cerr << "Error parsing input. Please enter a pair of indices or 'exit'." << std::endl;
	  }
	 
	}
      }
    }

    if (!output_file.empty()) {

      char buf[100000];
      if (pdb.number_of_models() != 1){
	sprintf(buf, output_file.c_str(), i);
      } else {
	sprintf(buf, output_file.c_str());
      }
      Magick::Geometry geom(arr.dim1(),arr.dim2());


      Magick::Image im(geom, "red");

      if (contact_map_threshold != std::numeric_limits<double>::infinity()){
	for (int j=0; j< arr.dim1(); ++j){
	  for (int k=0; k< arr.dim2(); ++k){
	    if (arr[j][k] < contact_map_threshold){
	      im.pixelColor(j, k, Magick::ColorGray(1));
	    } else {
	      im.pixelColor(j,k, Magick::ColorGray(0));
	    }
	  }
	}
      } else {
	double max=0;
	for (int i=0; i< arr.dim1(); ++i){
	  for (int j=0; j< arr.dim2(); ++j){
	    max=std::max(max, arr[i][j]);
	    //min=std::min(min, arr[i][j]);
	  }
	}
	//std::cout << "Maximum distance is " << max << std::endl;
	for (int i=0; i< arr.dim1(); ++i){
	  for (int j=0; j< arr.dim2(); ++j){
	    im.pixelColor(i,j, Magick::ColorGray(arr[i][j]/max));
	  }
	}
      }
	
      im.write(buf);
    }
  }
  return EXIT_SUCCESS;
  }

#else

int main(int, char*[]){
  bool this_program_requires_image_magick;
  std::cerr << "This program requires Image Magick++\n";
  return EXIT_FAILURE;
}

#endif
