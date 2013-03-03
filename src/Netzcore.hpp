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

#ifndef NETZCORE_HPP_9054862562754458964905368042156
#define NETZCORE_HPP_9054862562754458964905368042156

#include "ScoreNetwork.hpp"

//#include <boost/unordered_map.hpp>

#include <list>

//typedef boost::unordered_map<unsigned int,  Message> UIntToMessage;
///typedef boost::unordered_map<Vertex, UIntToMessage* > VertexToMessageMap;

template <class InputIterator>
std::pair<float, float> calculateMeanAndSigma(InputIterator first, InputIterator beyond);

class Netzcore: public ScoreNetwork {
public:
    Netzcore(); 
    Netzcore(std::string fileNode, std::string fileOutput, std::string prefixSampled, unsigned int nSample, bool fUseEdgeScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false); 
    Netzcore(std::string fileNode, std::string fileEdge, std::string fileOutput, std::string prefixSampled, unsigned int nSample, bool fUseEdgeScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false); 
    ~Netzcore();

    // methods below were protected before, made public to be able to use netzcore & netscore in combination
    void initializeScoring();
    void finalizeScoring();
    void initializeIteration();
    void finalizeIteration();
    void updateNodeScore(Vertex v);

private:
    //UIntToMessage * getVertexMessageMap(Vertex const v) { return messageMaps[v]; };
    //void createVertexMessageMap(Vertex const v) { messageMaps[v] = new UIntToMessage(); };

    float getVertexScoreInitial(Vertex const v) const { return mapScoreInitial.at(v); };
    void setVertexScoreInitial(Vertex const v, float vData) { mapScoreInitial[v] = vData; };

    void loadSampledGraphs();
    void printSampledGraphs();
    void updateSampledGraphScores();

    // MEMBERS
    //VertexToMessageMap messageMaps;
    VertexToFloat mapScoreInitial;
    std::list<Graph*> sampledGraphs;
    std::string nodeFileName;
    std::string prefixSampledGraphs;
    unsigned int nSampledGraphs;
    float mean;
    float sigma;
};

#endif // NETZCORE_HPP_

