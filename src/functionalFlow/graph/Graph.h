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

#ifndef GRAPH_H_
#define GRAPH_H_

#include "globals.h"
#include "hash.h"
#include "Vertex.h"
#include "Edge.h"

namespace fflow {

typedef std::tr1::unordered_map<int, Vertex*> MapIntToVertex;
//typedef hash_map<int, Vertex*> MapIntToVertex; // hash_map<char*, int, hash<char*>, eqstr> Map;
//typedef hash_map<string, int, hash<string>, isEqualString> MapStringToInt; 
typedef std::tr1::unordered_map<string, int> MapStringToInt;
//typedef hash_map<const char*, int, hash<const char*>, eqstr> MapStringToInt;
//typedef MapIntToVertex::value_type value_type;
//typedef MapIntToVertex::iterator iterator;
typedef std::tr1::unordered_map<int, string> MapIntToString;
//typedef hash_map<int, string> MapIntToString;

class Graph {
//friend ostream& operator<<(ostream& output, const Graph& g);
private:
void removeEdgeDirected(int vIdSource, int vIdTarget);
protected:
    MapIntToVertex *pMapIdToNode; // for mapping of node ids to Vertices
    MapStringToInt mapNameToId; // for mapping of node names to node ids
    MapIntToString mapIdToName; // for mapping of node ids to node names
    int idLast; // id assigned to last added node
    bool flagDirected;
    int nNode;
    int nEdge;    
public:
    Graph();
    Graph(bool fDirected);
    virtual ~Graph();
    //MapIntToVertex * getPMapIdToNode() const { return pMapIdToNode; }
    //MapIntToString getMapIdToName() const { return mapIdToName; }
    bool containsVertex(string name);
    bool containsVertex(int id);
    int addVertex();
    int addVertex(string name);
    void addVertex(int id, Data vData);
    int addVertex(string name, Data vData);
    void addEdge(int idSource, int idTarget);
    void addEdge(int idSource, int idTarget, Data eData);
    void addEdge(string nameSource, string nameTarget);
    void addEdge(string nameSource, string nameTarget, Data eData);
    void removeVertex(int vId);
    void removeEdge(int vIdSource, int vIdTarget); 
};

}

#endif /*GRAPH_H_*/
