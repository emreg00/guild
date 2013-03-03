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

#ifndef NETRANK_HPP_9054862562754458964905368042156
#define NETRANK_HPP_9054862562754458964905368042156

#include "graph/Graph.hpp"

class Netrank {
private:
    // MEMBERS
    Graph network;
    
    std::string outputFile;
    bool flagUseEdgeScore; // if true, page rank values of neighbors are scaled by edge weight 
    bool flagAccumulateToInitialNodeScore; // if true, accumulated shortest path score is added to node score

public:
    Netrank(); 
    Netrank(std::string fileNode, std::string fileEdge, std::string fileOutput, bool flagUseEdgeScore = true, bool flagAccumulateToInitialNodeScore = false);
    ~Netrank();
    Graph & getNetwork() { return network; };
    void run(unsigned int nIteration = 20, float dFactor = 0.85); 
};

#endif // NETRANK_HPP_

