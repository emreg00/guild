//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Author: Jeremy G. Siek
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#include <boost/config.hpp>
#include <iostream>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_utility.hpp>
#include "clustering_coefficient.hpp"
#include <map>


int main()
{
  using namespace boost;
  enum { A, B, C, D, E, F, N };
  const char* name = "ABCDEF";

  // A directed graph
  
  typedef adjacency_matrix<directedS> Graph;

  typedef boost::graph_traits < Graph >::vertex_descriptor Vertex;

  Graph g(N);
/*
  add_edge(B, C, g);
  add_edge(B, F, g);
  add_edge(C, A, g);
  add_edge(C, C, g);
  add_edge(D, E, g);
  add_edge(E, D, g);
  add_edge(F, A, g);
*/
  add_edge(A, B, g);
  add_edge(B, C, g);
  add_edge(C, D, g);
  add_edge(B, D, g);
  add_edge(D, E, g);

  std::cout << "vertex set: ";
  print_vertices(g, name);
  std::cout << std::endl;

  std::cout << "edge set: ";
  print_edges(g, name);
  std::cout << std::endl;

  std::cout << "out-edges: " << std::endl;
  print_graph(g, name);
  std::cout << std::endl;

  // An undirected graph

  typedef adjacency_matrix<undirectedS> UGraph;
  UGraph ug(N);
  //add_edge(B, C, ug);
  //add_edge(B, F, ug);
  //add_edge(C, A, ug);
  //add_edge(D, E, ug);
  //add_edge(F, A, ug);
  add_edge(A, B, ug);
  add_edge(B, C, ug);
  add_edge(C, D, ug);
  add_edge(B, D, ug);
  add_edge(D, E, ug);

  std::cout << "vertex set: ";
  print_vertices(ug, name);
  std::cout << std::endl;

  std::cout << "edge set: ";
  print_edges(ug, name);
  std::cout << std::endl;

  std::cout << "incident edges: " << std::endl;
  print_graph(ug, name);
  std::cout << std::endl;

  std::map<Vertex, float> vertexToCoefficient;
  associative_property_map< std::map<Vertex, float> > mapCoefficient(vertexToCoefficient);
  clustering_coefficient(g, mapCoefficient);

  std::map<Vertex, float>::iterator it, itEnd;
  for(it=vertexToCoefficient.begin(), itEnd=vertexToCoefficient.end(); it != itEnd; ++it)
	std::cout << it->first << ": " << it->second << std::endl;

  return 0;
}
