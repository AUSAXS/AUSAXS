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

#include <dsrpdb/PDB.h>

#include <vector>
#include <iterator>
#include <fstream>
#include <dsrpdb/align.h>
#include <dsrpdb/distance.h>


#ifdef PDB_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif

// return the cRMS
double align_to_points(const std::vector<dsrpdb::Point> &points,
		       char mode, dsrpdb::Protein &p){
  std::vector<dsrpdb::Point> ppoints;
  switch (mode) {
  case 'b':
    dsrpdb::backbone_coordinates(p.atoms_begin(), p.atoms_end(),
				 std::back_inserter(ppoints));
    break;
  case 'c':
    dsrpdb::ca_coordinates(p.atoms_begin(), p.atoms_end(),
			   std::back_inserter(ppoints));
    break;
  case 'a':
    dsrpdb::coordinates(p.atoms_begin(), p.atoms_end(),
			std::back_inserter(ppoints));
    break;
  };

  assert(points.size() == ppoints.size());

  dsrpdb::Transform tr= dsrpdb::compute_transform_taking_first_to_second(ppoints, points);
  //std::cout << tr << std::endl;

  for (dsrpdb::Protein::Atoms_iterator it= p.atoms_begin(); it != p.atoms_end(); ++it){
    dsrpdb::Point np= tr(it->second.cartesian_coords());
    it->second.set_cartesian_coords(np);
  }

  /*
  {
   
   dsrpdb::Transform tr2= dsrpdb::compute_transform_taking_first_to_second(points, 
									   ppoints);
   std::cout << tr2 << std::endl;
    
  }

  {
    std::vector<dsrpdb::Point> npoints(ppoints.size());
    for (unsigned int i=0; i< points.size(); ++i){
      dsrpdb::Point tp = tr(ppoints[i]);
      npoints[i]=tp;
    }
    dsrpdb::Transform tr2= dsrpdb::compute_transform_taking_first_to_second(npoints,
									    points);
    std::cout << tr2 << std::endl;
    
  }

  {
    std::vector<dsrpdb::Point> npoints(ppoints.size());
    for (unsigned int i=0; i< points.size(); ++i){
      dsrpdb::Point tp = tr(points[i]);
      npoints[i]=tp;
    }
    dsrpdb::Transform tr2= dsrpdb::compute_transform_taking_first_to_second(npoints,
									    ppoints);
    std::cout << tr2 << std::endl;
    
    }*/
  
  double dist=0;

  dsrpdb::Squared_distance sd;

  for (unsigned int i=0; i< points.size(); ++i){
    dsrpdb::Point tp = tr(ppoints[i]);
    dist += std::sqrt(sd(tp, points[i]));
  }
  return dist/ points.size();
}



int main(int argc, char *argv[]){
  std::string base_file, input_file, output_file;
  bool print_help=false;
  bool crms=false;
  char mode='b';
  bool warn=false;
#ifdef PDB_USE_BOOST_PROGRAM_OPTIONS
  {
    boost::program_options::options_description o("Allowed options"), po, ao;
    o.add_options()
      ("aligned-pdb,o", boost::program_options::value< std::string>(&output_file),
       "Where to write the result of aligning input-pdb to base-points.")
      ("atoms,a", boost::program_options::value< char>(&mode),
       "Which atoms to use from the protein. Values are b (backbone), c (c-alpha), a (all).")
      ("crms,c", boost::program_options::bool_switch(&crms),
       "Output the cRMS between the two pdbs (after alignment).")
      ("verbose,v", boost::program_options::bool_switch(&warn),
       "Warn about errors parsing pdb file.")
      ("help", boost::program_options::bool_switch(&print_help),
       "Produce help message");
    po.add_options()("base-points", boost::program_options::value< std::string>(&base_file),
		     "Base points")
      ("input-pdb", boost::program_options::value< std::string>(&input_file),
       "input file");
    ao.add(o).add(po);
    
    boost::program_options::positional_options_description p;
    p.add("base-points", 1);
    p.add("input-pdb", 1);

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				  options(ao).positional(p).run(), vm);
    boost::program_options::notify(vm);


    if (base_file.empty() || input_file.empty() || print_help
	|| (mode != 'b' && mode != 'c' && mode != 'a')) {
      std::cout << "This program aligns aligns a protein with a set of points and"
	" can compute distances between them.\n";
      std::cout << "usage: " << argv[0] << " base-points input-pdb\n\n";
      std::cout << o << "\n";
      return EXIT_FAILURE;
    }
  }

#else
  if (argc != 3){
    std::cerr << "The program is built without boost/program_options.hpp.\n";
    std::cerr << "useage: " << argv[0] << " file0.pdb file1.pdb" << std::endl;
    return EXIT_FAILURE;
  }
  input_file = argv[1];
  output_file = argv[2];
#endif
  
  std::ifstream in(input_file.c_str());

  if (!in){
    std::cerr<< "Error opening input file: " << input_file << std::endl;
    return EXIT_FAILURE;
  }
  std::ifstream bin(base_file.c_str());

  if (!bin){
    std::cerr<< "Error opening input file: " << base_file << std::endl;
    return EXIT_FAILURE;
  }

  dsrpdb::PDB input(in, warn);
  
  std::vector<dsrpdb::Point> points;
  while (bin) {
    char buf[10000];
    bin.getline(buf, 10000);
    if (!bin) break;
    if (buf[0]=='#') continue;
    std::istringstream iss(buf);
    dsrpdb::Point pt;
    iss >> pt;
    points.push_back(pt);
  }
  
  std::cout << "Read " << points.size() << " points.\n";

  if ( !output_file.empty() || crms) {
    for (unsigned int i=0; i< input.number_of_models(); ++i){
      dsrpdb::Model &m= input.model(i);
      dsrpdb::Protein &p= m.chain(0);
      double cRMS= align_to_points(points, mode, p);
      if (crms) {
	std::cout << cRMS << std::endl;
      }
    }
  }

  if (!output_file.empty()) {
    std::ofstream out(output_file.c_str());
    if (!out){
      std::cerr << "Error opening output file " 
		<< output_file << std::endl;
      return EXIT_FAILURE;
    }
    input.write(out);
  }

 

  return EXIT_SUCCESS;

}

