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

#ifndef NETSCORE_HPP_9054862562754458964905368042156
#define NETSCORE_HPP_9054862562754458964905368042156

#include "ScoreNetwork.hpp"

#include <boost/unordered_map.hpp>

// multiplication of weight of edges on the path, iteration of arrival, number of equidistant arrivals
typedef boost::tuple<float, int, int> Message;
typedef boost::unordered_map<unsigned int,  Message> UIntToMessage;
//typedef std::map<unsigned int,  Message> UIntToMessage;
typedef boost::unordered_map<Vertex, UIntToMessage* > VertexToMessageMap;
//typedef std::map<Vertex, UIntToMessage* > VertexToMessageMap;

class Netscore: public ScoreNetwork {
public:
    Netscore(); 
    Netscore(std::string fileOutput, bool fUseEdgeScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false); 
    Netscore(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fUseEdgeScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false); 
    ~Netscore();
    //setNetwork(const Graph * g);

    // methods below were protected before, made public to be able to use netzcore & netscore in combination
    void initializeScoring();
    void finalizeScoring();
    void initializeRepetition();
    void finalizeIteration();
    void updateNodeScore(Vertex v);

private:
    UIntToMessage * getVertexMessageMap(Vertex const v) { return messageMaps[v]; };
    void createVertexMessageMap(Vertex const v) { messageMaps[v] = new UIntToMessage(); };
    float getVertexScoreInitial(Vertex const v) const { return mapScoreInitial.at(v); };
    void setVertexScoreInitial(Vertex const v, float vData) { mapScoreInitial[v] = vData; };
    float getVertexScoreInitialOriginal(Vertex const v) const { return mapScoreInitialOriginal.at(v); };
    void setVertexScoreInitialOriginal(Vertex const v, float vData) { mapScoreInitialOriginal[v] = vData; };

    // MEMBERS
    VertexToMessageMap messageMaps;
    VertexToFloat mapScoreInitial;
    VertexToFloat mapScoreInitialOriginal;
};

#endif // NETSCORE_HPP_

