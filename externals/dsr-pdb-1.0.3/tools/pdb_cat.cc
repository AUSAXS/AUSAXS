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
  std::vector<std::string> input_files;
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
      ("input-pdbs", boost::program_options::value< std::vector<std::string> >(&input_files),
       "The input and output files.");

    ao.add(o).add(po);
    
    boost::program_options::positional_options_description p;
    p.add("input-pdbs", -1);

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				  options(ao).positional(p).run(), vm);
    boost::program_options::notify(vm);

    
    if (input_files.size() <2 || print_help) {
      std::cout << "Concatenate a bunch of pdb files into one pdb file with many models.\n";
      std::cout << "usage: " << argv[0] << " input-pdb-0 input-pdb-1 ... output-pdb\n";
      std::cout << o << "\n";
      return EXIT_FAILURE;
    }

    output_file=input_files.back();
    input_files.pop_back();
  }

#else
  for (int i=1; i< argc-1; ++i){
    input_files.push_back(argv[i]);
  }
  output_file = argv[argc-1];
#endif

  
  dsrpdb::PDB outpdb;

  for (unsigned int i=0; i < input_files.size(); ++i){
    std::ifstream in(input_files[i].c_str());
    if (!in){
      std::cerr << "Error opening input file " << input_files[i] << std::endl;
      return EXIT_FAILURE;
    }
    dsrpdb::PDB inpdb(in, verbose);
    
    for (dsrpdb::PDB::Models_iterator mit= inpdb.models_begin(); mit != inpdb.models_end(); ++mit){
      dsrpdb::Model m= *mit;
      m.set_index(outpdb.number_of_models());
      outpdb.new_model(m);
    }
  }
  
  if (output_file.empty()){
    outpdb.write(std::cout);
  } else {
    std::ofstream out(output_file.c_str());
    outpdb.write(out);
  }


  return EXIT_SUCCESS;
}
