/*
    GUILD (Genes Underlying Inheritance Linked Disorders) implements several 
    graph based algorithms for scoring relevance of a node in the network in 
    terms of a phenotype using known associations in the node's neighborhood 
    for that phenotype. GUILD has been applied to the prioritization of genes 
    for several human disorders. 2011 - Emre Guney (Unviersitat Pompeu Fabra)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef BOOST_GRAPH_CLUSTERING_COEFFICIENT_HPP
#define BOOST_GRAPH_CLUSTERING_COEFFICIENT_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/unordered_set.hpp>

#ifdef BOOST_NO_TEMPLATED_ITERATOR_CONSTRUCTORS
#  include <iterator>
#endif

/*  Author: Emre Guney, 2009, Pompeu Fabra University

    This algorithm is to calculate clustering coefficients (c.c.) of each vertex in a graph
    The c.c. vertex v will be stored in coefficient[v].

    Algorithm: 
    Let G = (V,E) be a graph with vertices. Where (u,v) in E, for each v in V, algorithm calculates c.c.
    based on the following formula:
	c.c. = # of edges connecting neighbors of v / Combination(# of neighbors, 2) 
	Therefore,
	For undirected graphs [(u,v) and (v,u) considered identical]:
	c.c. = |(u,v)| / ( |u| * (|u| - 1) / 2 )

	For directed graphs:
	c.c. = |(u,v)| / ( |u| * (|u| - 1) )

    Reference:
    D. J. Watts and Steven Strogatz (June 1998). "Collective dynamics of 'small-world' networks". Nature 393: 440–442. doi:10.1038/30918
*/

namespace boost {
    
    template <typename VertexListGraph, typename CoefficientMap>
    void 
    clustering_coefficient(const VertexListGraph& G, CoefficientMap coefficient)
    {
	typedef graph_traits<VertexListGraph> GraphTraits;
	typedef typename GraphTraits::vertex_descriptor Vertex;
	typedef typename property_traits<CoefficientMap>::value_type size_type;

        function_requires< VertexListGraphConcept<VertexListGraph> >();
        function_requires< ReadWritePropertyMapConcept<CoefficientMap, Vertex> >();

	typename GraphTraits::vertex_iterator vt, vt_end;
	typename GraphTraits::adjacency_iterator ut, ut_end;
	typename GraphTraits::out_edge_iterator et, et_end;

	typedef typename GraphTraits::directed_category Directed;

	// For each vertex
	for(tie(vt, vt_end) = vertices(G); vt != vt_end; ++vt) 
	{
	    // Get neighbor vertices
	    unordered_set<Vertex> neighbor_vertices;
	    for(tie(ut, ut_end) = adjacent_vertices(*vt, G); ut != ut_end; ++ut) 
	    {
		if(*ut != *vt)
		    neighbor_vertices.insert(*ut);
	    }
	    // Calculate clustering coefficient applying above formula
	    typename unordered_set<Vertex>::const_iterator it, it_end;
	    unsigned int N=neighbor_vertices.size(), n=0;
	    size_type clustering_coefficient = 0.0;
	    // For each neighbor of the vertex
	    for(it=neighbor_vertices.begin(), it_end=neighbor_vertices.end(); it != it_end; ++it) 
	    {
		// Count number of connected neighbors
		for(tie(et, et_end)=out_edges(*it, G); et != et_end; ++et) 
		{
			if(neighbor_vertices.find(target(*et, G)) != it_end) {
			    n += 1;
			}
		}
	    }
	    // Taking this condition below the loop would improve performance at a cost of reduced readability
	    if(N == 1 || n==0) {
		clustering_coefficient = 0;
	    } else {
		clustering_coefficient = float(n) / (N * (N-1));
	    }
	    // Save the coefficient of current vertex 
	    put(coefficient, *vt, clustering_coefficient);  
	}
	return;
    }
}

#endif
