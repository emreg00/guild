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

#include "Vertex.h"

using namespace fflow;

Vertex::Vertex() {
	id = 0;
	degree = 0;
	//cClustering = 0;
	pMapIdToEdge = new MapIntToEdge;
}

Vertex::Vertex(int vId) {
	id = vId;
	degree = 0;
	//cClustering = 0;
	pMapIdToEdge = new MapIntToEdge;
}

Vertex::Vertex(int vId, Data vData) {
	id = vId;
	data = vData;
	degree = 0;
	//cClustering = 0;
	pMapIdToEdge = new MapIntToEdge;
}

Vertex::~Vertex() {
	// iterate over Edges and remove them
	MapIntToEdge::iterator it = pMapIdToEdge->begin(), itEnd = pMapIdToEdge->end();
	for(; it != itEnd; ++it)
		delete it->second;
	// removeMap
	delete pMapIdToEdge;
}

bool Vertex::containsNeighbor(int id) {
	MapIntToEdge::iterator it;
	it = pMapIdToEdge->find(id);
	if(it != pMapIdToEdge->end()) {
	    return true;
	}
	return false;
}

void Vertex::addNeighbor(int vId) {
	if(!containsNeighbor(vId)) {
		(*pMapIdToEdge)[vId] = new Edge(id, vId);
		degree++;
	} else {
		cerr << "Warning: Neighbor already contained (ignoring) " << endl;
	}
}

void Vertex::addNeighbor(int vId, Data eData) {
	if(!containsNeighbor(vId)) {
		(*pMapIdToEdge)[vId] = new Edge(id, vId, eData);
		degree++;
	} else {
			cerr << "Warning: Neighbor already contained (ignoring) " << endl;
	}
}

void Vertex::removeNeighbor(int vId){
	MapIntToEdge::iterator it;
	it = pMapIdToEdge->find(vId);
	if(it != pMapIdToEdge->end()) { //it->second != NULL) {
		delete it->second;
		//it->second = NULL;
		pMapIdToEdge->erase(it);
		degree--;
	} else {
		cerr << "Warning: trying to remove an unexisting neighbor" << endl;
	}
}

ostream& fflow::operator<<(ostream& output, const Vertex& v) {
    output << "(" <<  v.id << "{" << v.data.reserve << "}: ";
	MapIntToEdge::iterator it, itEnd;
	MapIntToEdge *map = v.pMapIdToEdge;
	for(it = map->begin(), itEnd = map->end(); it != itEnd; ++it)
		output << it->first << ", ";
	output << ")";
    return output; 
}

