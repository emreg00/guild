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

#ifndef FUNCTIONALFLOW_H_
#define FUNCTIONALFLOW_H_

#include "graph/globals.h"
//#include "hash.h"
//#include "Vertex.h"
//#include "Edge.h"
#include "graph/Graph.h"
#include <cmath>

using namespace fflow;

class FunctionalFlow: public fflow::Graph {
private:    
    bool flagDebug;
    bool flagRemoveSeedEffect;
    float seedScoreThreshold;
    std::unordered_map<int, void*> mapSource;
    //hash_map<int, void*> mapSource;
    string fileOutput;
public:
    FunctionalFlow(); 
    FunctionalFlow(string fileProtein, string fileInteraction, string fileOutput, float seedScoreThreshold = 0.01, bool flagRemoveSeedEffect = false, bool flagDebug = false);
    virtual ~FunctionalFlow();
    void loadProteins(string const &fileName);
    void loadInteractions(string const &fileName);
    void filterNetwork(int degreeNodeMin = -INT_INFINITY, int degreeNodeMax = INT_INFINITY, float capacityMin = -INFINITY, bool flagRemoveSelfEdges = true);
    void run(int nIteration);
    void updateNetwork();
    void printNetwork(string explanation = "");
    void printNodes(string explanationPrefix = "");
    void printEdges();
    void printNodesAndScores();
    void outputScores();
};

#endif // FUNCTIONALFLOW_H_

