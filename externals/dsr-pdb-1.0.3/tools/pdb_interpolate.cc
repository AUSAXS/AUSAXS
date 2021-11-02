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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <queue> 


struct atom_label_t {
  typedef boost::vertex_property_tag kind;
};
    
typedef int Atom_label;


//typedef boost::property<boost::edge_weight_t, double> Edge_property;
typedef boost::property<atom_label_t, Atom_label > Vertex_property;
typedef boost::adjacency_list<boost::vecS, boost::vecS,
			      boost::undirectedS, Vertex_property> Protein_graph;
typedef boost::property_map<Protein_graph, atom_label_t>::type Atom_label_map;
//typedef boost::property_map<Protein_graph, color_map_t>::type Color_map;
typedef boost::graph_traits<Protein_graph> Protein_graph_traits;
typedef Protein_graph_traits::vertex_iterator Vertex_iterator;
typedef Protein_graph_traits::vertex_descriptor Vertex_descriptor;
typedef Protein_graph_traits::edge_descriptor Edge_descriptor;



struct Spherical_point {
  Spherical_point(){
    for (unsigned int i=0; i< 3; ++i){ p_[i]= std::numeric_limits<double>::infinity();}
  }
  Spherical_point(double phi, double psi, double r){
    p_[0]=phi;
    p_[1]=psi;
    p_[2]=r;
  }

  double phi() const {
    return p_[0];
  }
  double psi() const {
    return p_[1];
  }
  double r() const {
    return p_[2];
  }

  operator bool() const {
    return p_[0] != std::numeric_limits<double>::infinity();
  }

  double p_[3];
};

template <class It, class Oit>
void spherical_average(It b0, It e0, It b1, It, double f, Oit out) {
  for (It c0=b0, c1=b1; c0 != e0; ++c0, ++c1){
    *out= Spherical_point(f*c0->phi()+ (1-f) * c1->phi(),
			  f*c0->psi()+ (1-f) * c1->psi(),
			  f*c0->r()+ (1-f) * c1->r());
    ++out;
  }
}


void protein_to_graph(const dsrpdb::Protein &p, Protein_graph &g){
  std::vector<Vertex_descriptor> index_vertex_map;
  Atom_label_map alm= boost::get(atom_label_t(), g);
  for (dsrpdb::Protein::Const_atoms_iterator ait= p.atoms_begin(); ait != p.atoms_end(); ++ait){
    Vertex_descriptor v= boost::add_vertex(g);
    boost::put(alm, v, ait->index());
    if (ait->index() > index_vertex_map.size()-1){
      index_vertex_map.resize(ait->index()+1);
    }
    index_vertex_map[ait->index()]= v;
  }

  for (dsrpdb::Protein::Bonds_iterator rit= p.bonds_begin(); rit != p.bonds_end(); ++rit){
    boost::add_edge(index_vertex_map[rit->first], index_vertex_map[rit->second], g);
  }
  
}

struct Compare_vertices {
  Compare_vertices(const Atom_label_map &alm): map_(alm){}
  bool operator()(Vertex_descriptor va, Vertex_descriptor vb) const {
    return map_[va] < map_[vb];
  }
};

void pred_map(const Protein_graph &g, std::vector<Vertex_descriptor> &preds){
  Compare_vertices comp(boost::get(atom_label_t(), g))
  std::priority_queue<Vertex_descriptor, std::vector<Vertex_descriptor>, 
    Compare_vertices> q(comp);

  q.insert(*std::min_element(boost::vertices(g).first,
			     boost::vertices(g).second,
			     q));
  while (!q.empty()){
    Vertex_descriptor c= q.top();
    q.pop();
    
  }
  
  
}

/*void frame(Vertex_descriptor v, const Protein_graph &g, const dsrpdb::Protein &p, 
	   dsrpdb::Point &center, dsrpdb::Vector &e0, dsrpdb::Vector &e1) {
  
	   }*/

Spherical_point cartesian_to_spherical(Vertex_descriptor v, const Protein_graph &g, const dsrpdb::Protein &p) {
  
  /*dsrpdb::Point
  std::pair<dsrpdb::Point, dsrpdb::Point> ps= parents(v, g);
  dsrpdb::Point tp= g.atom(boost::get(atom_label_t(), g)[v]).cartesian_coordinates();*/
  
}

dsrpdb::Point spherical_to_cartesian(const Spherical_point &pt, Vertex_descriptor v, const Protein_graph &g, const dsrpdb::Protein &p) {
  std::vector<Vertex_descriptor> preds;
}


// The point ordering/labeling are kind of bad
template <class Oit>
void cartesian_to_spherical(const dsrpdb::Protein &p, Oit out){
  if (p.number_of_residues()==0) return;

  Protein_graph g;
  protein_to_graph(p,g);
  std::vector<Vertex_descriptor> preds;
  pred_map(g, preds);
  for (Vertex_iterator it= boost::vertices(g).first; it != boost::vertices(g).second; ++it){
    *out= cartesian_to_spherical(*it, g, p);
    ++out;
  }
}

template <class It>
void spherical_to_cartesian(It b, It e, const Protein_graph &g, dsrpdb::Protein &p){
  It c= b;
  std::vector<Vertex_descriptor> preds;
  pred_map(g, preds);
  for (Vertex_iterator it= boost::vertices(g).first; it != boost::vertices(g).second; ++it, ++c){
    dsrpdb::Point pt= spherical_to_cartesian(*c, *it, g, p);
    it->set_cartesian_coordinates(pt);
  }
}


int main(int argc, char *argv[]){
  std::string input_file_0, input_file_1, output_file;
  bool print_help=false;
  bool verbose=false;

#ifdef PDB_USE_BOOST_PROGRAM_OPTIONS
  boost::program_options::options_description o("Allowed options"), po, ao;
  o.add_options()
    ("help", boost::program_options::bool_switch(&print_help), "produce help message");
   
  po.add_options()
    ("input0-pdb", boost::program_options::value< std::string>(&input_file_0), "input file")
    ("input1-pdb", boost::program_options::value< std::string>(&input_file_1), "input file")
    ("output-pdb", boost::program_options::value< std::string>(&output_file), "output file");
  ao.add(o).add(po);

  boost::program_options::positional_options_description p;
  p.add("input0-pdb", 1);
  p.add("input1-pdb", 1);
  p.add("output-pdb", 1);

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				options(ao).positional(p).run(), vm);
  boost::program_options::notify(vm);


 
  //boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  //boost::program_options::notify(vm);    

  if (input_file_0.empty() || input_file_1.empty() || output_file.empty() || print_help) {
    std::cout << "This program interpolates between two pdb files.\n";
    std::cout << "usage: " << argv[0] << " input0-pdb input1-pdb output-pdb\n\n";
    std::cout << o << "\n";
    return EXIT_FAILURE;
  }

#else
  if (argc != 4){
    std::cerr << "The program is built without boost/program_options.hpp.\n";
    std::cerr << "useage: " << argv[0] << " [-c] file.pdb output.pdb" << std::endl;
    return EXIT_FAILURE;
  }
  input_file_0 = argv[1];
  input_file_1 = argv[2];
  output_file = argv[3];
#endif


  std::ifstream in0(input_file_0.c_str());
  if (!in0) {
    std::cerr << "Error opening input file " << input_file_0 << std::endl;
    return EXIT_FAILURE;
  }
  std::ifstream in1(input_file_1.c_str());
  if (!in1) {
    std::cerr << "Error opening input file " << input_file_1 << std::endl;
    return EXIT_FAILURE;
  }

  dsrpdb::PDB pdb0(in0, verbose);
  dsrpdb::PDB pdb1(in1, verbose);

  if (pdb0.number_of_models() != pdb1.number_of_models()){
    std::cerr << "The two pdbs must have the same number of models.\n";
    return EXIT_FAILURE;
  }
 
  dsrpdb::PDB pdba=pdb0;

  for (unsigned int i=0;i< pdb0.number_of_models(); ++i){
    dsrpdb::Model &m0= pdb0.model(i);
    dsrpdb::Model &m1= pdb1.model(i);
    dsrpdb::Model &ma= pdba.model(i);
    if (m0.number_of_chains() != m1.number_of_chains()){
      std::cerr << "The two pdbs must have the same number of chains in each model.\n";
      return EXIT_FAILURE;
    }

    for (unsigned int j=0; j< m0.number_of_chains(); ++j){
      dsrpdb::Protein &p0= m0.chain(j);
      dsrpdb::Protein &p1= m1.chain(j);
      dsrpdb::Protein &pa= ma.chain(j);

      if (p0.number_of_residues() != p1.number_of_residues()){
	std::cerr << "The two pdbs must have the same number of residues in each chain.\n";
	return EXIT_FAILURE;
      }
      
      // make this better later
      std::vector<Spherical_point> points_0, points_1;

      cartesian_to_spherical(p0, std::back_inserter(points_0));
      cartesian_to_spherical(p1, std::back_inserter(points_1));

      std::vector<Spherical_point> average_points;

      spherical_average(points_0.begin(), points_0.end(),
			points_1.begin(), points_1.end(),
			.5, std::back_inserter(average_points));
      spherical_to_cartesian(average_points.begin(), average_points.end(), pa);
    }
  }

  std::ofstream out(output_file.c_str());
  if (!out) {
    std::cerr << "Error opening output file " << output_file << std::endl;
    return EXIT_FAILURE;
  }

  pdba.write(out);

  //delete[] buf;
  return EXIT_SUCCESS;
}
