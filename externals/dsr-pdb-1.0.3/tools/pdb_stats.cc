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

#ifdef PDB_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif

int main(int argc, char *argv[]){
  std::string input_file;
  std::string output_file;
  bool print_help=false;
  bool verbose=false;

#ifdef PDB_USE_BOOST_PROGRAM_OPTIONS
  {
    boost::program_options::options_description o("Allowed options"), po, ao;
    o.add_options()
      ("help", boost::program_options::bool_switch(&print_help), "produce help message")
      ("verbose,v", boost::program_options::bool_switch(&verbose), 
       "Print error messages from reading the pdb files.");
    po.add_options()
      ("input-pdb", boost::program_options::value< std::string >(&input_file),
       "The input and output files.");

    ao.add(o).add(po);
    
    boost::program_options::positional_options_description p;
    p.add("input-pdb", -1);

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				  options(ao).positional(p).run(), vm);
    boost::program_options::notify(vm);

    
    if (print_help) {
      std::cout << "Print some statistics of a pdb file.\n";
      std::cout << "usage: " << argv[0] << " pdb_file.pdb\n";
      std::cout << o << "\n";
      return EXIT_FAILURE;
    }
  }

#else
  input_file=argv[1];
#endif

  std::ifstream in(input_file.c_str());
  if (!in){
    std::cerr << "Error opening input file " << input_file << std::endl;
    return EXIT_FAILURE;
  }
  
  dsrpdb::PDB inpdb(in, verbose);
    
  std::cout << "PDB has " << inpdb.number_of_models() << " models." << std::endl;

  for (dsrpdb::PDB::Models_iterator mit= inpdb.models_begin(); mit != inpdb.models_end(); ++mit){
    dsrpdb::Model m= *mit;
    std::cout << "Model has " << m.number_of_chains() << " chains." << std::endl;
    for (unsigned int i=0; i< m.number_of_chains(); ++i){
      std::cout << m.chain(i).chain() << " ";
    }
    std::cout << std::endl;
  }
  
 
  
  return EXIT_SUCCESS;
}
