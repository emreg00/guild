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

#include "Graph.h"

using namespace fflow;

Graph::Graph() {
	pMapIdToNode = new MapIntToVertex;
	idLast = 0;
    flagDirected = false;
    nNode = 0;
    nEdge = 0;
}

Graph::Graph(bool fDirected) {
	pMapIdToNode = new MapIntToVertex;
	idLast = 0;
    flagDirected = fDirected;
    nNode = 0;
    nEdge = 0;
}

Graph::~Graph() {
	MapIntToVertex::iterator it = pMapIdToNode->begin(), itEnd = pMapIdToNode->end();
	for(; it != itEnd; ++it)
		delete it->second;
	mapIdToName.clear();
	mapNameToId.clear();
}

bool Graph::containsVertex(string name) {
	MapStringToInt::iterator iterator;
	//iterator = mapNameToId.find(name.c_str());
	iterator = mapNameToId.find(name);
	if(iterator != mapNameToId.end()) {
	    return true;
	}
	return false;
}

bool Graph::containsVertex(int id) {
	MapIntToString::iterator iterator;
	iterator = mapIdToName.find(id);
	if(iterator != mapIdToName.end()) {
	    return true;
	}
	return false;
}

int Graph::addVertex() {
	idLast++;
	(*pMapIdToNode)[idLast] = new Vertex(idLast);
	nNode++;
    return idLast;
}

int Graph::addVertex(string name) {
	if(!containsVertex(name)) {
		idLast++;
		(*pMapIdToNode)[idLast] = new Vertex(idLast);
		//mapNameToId[name.c_str()] = idLast; 
		mapNameToId[name] = idLast; 
		mapIdToName[idLast] = name;
		nNode++;
		return idLast;
	}
	cerr << "Warning: Vertex already contained (ignoring) " << endl;
	//return mapNameToId[name.c_str()];
	return mapNameToId[name];
}

void Graph::addVertex(int id, Data data) {
	if(!containsVertex(id)) { 
		(*pMapIdToNode)[id] = new Vertex(id, data);
		nNode++;
	} else {
		cerr << "Warning: Vertex already contained (ignoring) " << endl;
	}
}

int Graph::addVertex(string name, Data data) {
	if(!containsVertex(name)) {
		idLast++;
	    (*pMapIdToNode)[idLast] = new Vertex(idLast, data); 
		//mapNameToId[name.c_str()] = idLast; 
		mapNameToId[name] = idLast; 
		mapIdToName[idLast] = name;
		nNode++;
		return idLast;
	}
	cerr << "Warning: Vertex already contained (ignoring) " << endl;
	//return mapNameToId[name.c_str()];
	return mapNameToId[name];
}

void Graph::addEdge(int idSource, int idTarget) {
	if(containsVertex(idSource) && containsVertex(idTarget)) {
		(*pMapIdToNode)[idSource]->addNeighbor(idTarget);
	    if(!flagDirected and idSource != idTarget)
		    (*pMapIdToNode)[idTarget]->addNeighbor(idSource);
	    nEdge++;
	} else {
		cerr << "Warning: At least one of the vertex is missing " << endl;
	}
}

void Graph::addEdge(int idSource, int idTarget, Data eData) {
	if(containsVertex(idSource) && containsVertex(idTarget)) {
		(*pMapIdToNode)[idSource]->addNeighbor(idTarget, eData);
	    if(!flagDirected and idSource != idTarget)
		    (*pMapIdToNode)[idTarget]->addNeighbor(idSource, eData);
	    nEdge++;
	} else {
		cerr << "Warning: At least one of the vertices is missing " << endl;
	}
}

void Graph::addEdge(string nameSource, string nameTarget) {
	if(containsVertex(nameSource) && containsVertex(nameTarget)) {
	    int idSource, idTarget;
	    //idSource = mapNameToId[nameSource.c_str()];
	    idSource = mapNameToId[nameSource];
	    //idTarget = mapNameToId[nameTarget.c_str()];
	    idTarget = mapNameToId[nameTarget];
	    (*pMapIdToNode)[idSource]->addNeighbor(idTarget);
	    if(!flagDirected and idSource != idTarget)
		    (*pMapIdToNode)[idTarget]->addNeighbor(idSource);
	    nEdge++;
	} else {
		cerr << "Warning: At least one of the vertices is missing " << endl;
	}
}

void Graph::addEdge(string nameSource, string nameTarget, Data eData) {
	if(containsVertex(nameSource) && containsVertex(nameTarget)) {
	    int idSource, idTarget;
	    //idSource = mapNameToId[nameSource.c_str()];
	    idSource = mapNameToId[nameSource];
            //idTarget = mapNameToId[nameTarget.c_str()];
            idTarget = mapNameToId[nameTarget];
	    (*pMapIdToNode)[idSource]->addNeighbor(idTarget, eData);
	    if(!flagDirected and idSource != idTarget)
		(*pMapIdToNode)[idTarget]->addNeighbor(idSource, eData);
	    nEdge++;
	} else {
		cerr << "Warning: At least one of the vertices is missing " << endl;
	}
}

void Graph::removeVertex(int vId) {
    Vertex *pVertex = (*pMapIdToNode)[vId];
    if(!flagDirected) {
        // remove neighbors edges to this edge
    	MapIntToEdge::iterator it = pVertex->pMapIdToEdge->begin(), itEnd = pVertex->pMapIdToEdge->end();
    	for(; it != itEnd; ++it) {
            removeEdgeDirected(it->first, vId); // (*pMapIdToNode)[it->first]->removeNeighbor(vId); // removeEdge(vId, it->first); 
            nEdge--; //!! should be erased if removeEdge is used. 
        }
    }
    // remove node from the graph
    delete pVertex; // does erase call destructor somehow maybe with delete??
    pMapIdToNode->erase(vId);
    string name = mapIdToName[vId];
	//mapNameToId.erase(name.c_str());
	mapNameToId.erase(name);
	mapIdToName.erase(vId);
	nNode--;
	//cout << " size: " << pMapIdToNode->size() << endl;
	//cout<< "removeVertex " << *this << endl;
}

void Graph::removeEdge(int vIdSource, int vIdTarget) {
	removeEdgeDirected(vIdSource, vIdTarget);
	if(!flagDirected and vIdSource != vIdTarget) {
		removeEdgeDirected(vIdTarget, vIdSource);
	}
	nEdge--;
	// Vertex *pVertex = (*pMapIdToNode)[vId];
	// if(pVertex->pMapIdToEdge->size() == 0) removeVertex(vId) 
}

// nEdge is not modified since this method supposed to be called only from removeEdge or removeVertex 
void Graph::removeEdgeDirected(int vIdSource, int vIdTarget) {
	Vertex *pVertex = (*pMapIdToNode)[vIdSource];
	pVertex->removeNeighbor(vIdTarget);
	// cout << *pVertex << " removed edge with: " << vIdTarget << endl;
}

/*
ostream& operator<<(ostream& output, const Graph& g) {
	//MapIntToVertex *pMap = g.pMapIdToNode;
	MapIntToVertex *pMap = g.getPMapIdToNode();
	MapIntToVertex::iterator it, itEnd;
	//MapIntToString map2 = g.mapIdToName;
        MapIntToString map2 = g.getMapIdToName();
	itEnd = pMap->end(); 
	for(it = pMap->begin(); it != itEnd; ++it)
		output << "[" << map2[it->first] << "]" << *(it->second) << " - ";
    output << endl;
    return output;  
}
*/

