#ifndef BOOST_GRAPH_CLUSTERING_COEFFICIENT_HPP
#define BOOST_GRAPH_CLUSTERING_COEFFICIENT_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map.hpp>
#include <boost/unordered_set.hpp>

#ifdef BOOST_NO_TEMPLATED_ITERATOR_CONSTRUCTORS
#  include <iterator>
#endif

/* 
    This algorithm is to calculate clustering coefficients (c.c.) of each vertex in a graph
    The c.c. vertex v will be stored in coefficient[v].

    Algorithm: 
    Let G = (V,E) be a graph with vertices. Where (u,v) in E, for each v in V, algorithm calculates c.c.
    based on the following formula:
	For undirected graphs:
	c.c. = # of edges connecting neighbors of v / Combination(# of neighbors, 2) = |(u,v)| / ( |u| * (|u| - 1) / 2 )

	For directed graphs:
	c.c. = |(u,v)| / ( |u| * (|u| - 1) )

    Reference:
    D. J. Watts and Steven Strogatz (June 1998). "Collective dynamics of 'small-world' networks". Nature 393: 440–442. doi:10.1038/30918
*/

namespace boost {
    
    template <typename VertexListGraph, typename CoefficientMap, typename VertexIndexMap>
    void 
    clustering_coefficient(const VertexListGraph& G, CoefficientMap coefficient, VertexIndexMap index)
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
	    for(it=neighbor_vertices.begin(), it_end=neighbor_vertices.end(); it != it_end; ++it) 
	    {
		for(tie(et, et_end)=out_edges(*it, G); et != et_end; ++et) 
		{
		    //if(is_convertible<Directed*, directed_tag*>::value == false) {
		    //	if(index[*it] > index[target(*et, G)]) {
		    //	    if(neighbor_vertices.find(target(*et, G)) != it_end) {
		    //		n += 1;
		    //	    }
		    //	}
		    //} else {
			if(neighbor_vertices.find(target(*et, G)) != it_end) {
			    n += 1;
			}
		    //}
		}
	    }
	    if(N == 1 || n==0) {
		clustering_coefficient = 0;
	    } else {
		clustering_coefficient = float(n) / (N * (N-1));
	    }
	    // Clustering coefficient differs in case of directed/undirected graph
	    //if(is_convertible<Directed*, directed_tag*>::value == false) 
	    //	clustering_coefficient *= 2;

	    // Save the coefficient of current vertex 
	    put(coefficient, *vt, clustering_coefficient);  
	}
	return;
    }

    template <typename VertexListGraph, typename CoefficientMap>
    void 
    clustering_coefficient(const VertexListGraph& G, CoefficientMap coefficient)
    {
	return clustering_coefficient(G, coefficient, get(vertex_index, G));
    }
}

#endif
