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

#include "Netrank.hpp"

#include <iostream>

Netrank::Netrank() 
{
    flagAccumulateToInitialNodeScore = false;
}

Netrank::Netrank(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fUseEdgeScore, bool fAccumulateToInitialNodeScore)
{
    flagUseEdgeScore = fUseEdgeScore;
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore;
    network.loadNodes(fileNode); 
    network.loadEdges(fileEdge); 
    outputFile = fileOutput;
}

Netrank::~Netrank() {
}

void Netrank::run(unsigned int nIteration, float dFactor) 
{ 
    std::map<Vertex, float>::iterator vt, vtEnd;
    std::map<Vertex, float> mapRank;
    //getNetwork().calculatePageRank(mapRank, nIteration, dFactor);
    //getNetwork().calculatePageRankWithPriors(mapRank, nIteration, dFactor);
    getNetwork().calculatePageRankWithWeightedPriors(mapRank, nIteration, dFactor);
    float score = 0.0;
    for(vt = mapRank.begin(), vtEnd = mapRank.end(); vt != vtEnd; ++vt) 
    {
	score = 0.0;
	if(flagAccumulateToInitialNodeScore) {
	    score = getNetwork().getVertexScore(vt->first);
	}
	score += vt->second;
	getNetwork().setVertexScore(vt->first, score);
    }
    getNetwork().scaleVertexScores(SCALE_BETWEEN_ZERO_AND_ONE);
    getNetwork().outputScores(outputFile);
    return;
}


