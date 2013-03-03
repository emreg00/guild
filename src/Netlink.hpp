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

#ifndef NETLINK_HPP_9054862562754458964905368042156
#define NETLINK_HPP_9054862562754458964905368042156

#include "graph/Graph.hpp"
#include "ScoreNetwork.hpp"

class Netlink: public ScoreNetwork {
private:
    // MEMBERS
    float threshold;
    // inherited from ScoreNetwork
    //Graph network;
    //std::string outputFile;
    //bool flagAccumulateToInitialNodeScore; // if true, linker degree score is added to node score
    //bool flagVerbose;

public:
    Netlink(); 
    Netlink(std::string fileNode, std::string fileEdge, std::string fileOutput, float t, bool flagAccumulateToInitialNodeScore = false, bool flagVerbose = false);
    ~Netlink();
    //Graph & getNetwork() { return network; };
    //void run(int threshold); 
    void updateNodeScore(Vertex v);
};

#endif // NETLINK_HPP_

