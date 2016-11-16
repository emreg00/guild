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

#ifndef VERTEX_H_
#define VERTEX_H_

#include "hash.h"
#include "globals.h"
#include "Edge.h"

namespace fflow {

typedef std::unordered_map<int, Edge*> MapIntToEdge;
//typedef hash_map<int, Edge*> MapIntToEdge;
//typedef MapIntToEdge::value_type value_type;
//typedef MapIntToEdge::iterator iterator;

class Vertex
{
friend ostream& operator<<(ostream& output, const Vertex& v);
public:
	int id;
	Data data;
	int degree;
	//float cClustering;
	MapIntToEdge *pMapIdToEdge; // neighbor id to edge
//public:
	Vertex();
	Vertex(int id);
	Vertex(int vId, Data vData);
	virtual ~Vertex();
	bool containsNeighbor(int id);
	void addNeighbor(int vId);
	void addNeighbor(int vId, Data eData);
	void removeNeighbor(int vId);
};

}

#endif /*VERTEX_H_*/
