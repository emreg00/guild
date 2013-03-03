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

#include "Netshort.hpp"

#include <iostream>

Netshort::Netshort() 
{
    flagAccumulateToInitialNodeScore = false;
}

Netshort::Netshort(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fAccumulateToInitialNodeScore)
{
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore;
    network.loadNodes(fileNode); 
    network.loadEdges(fileEdge, true); // flag causes edge scores to be inverted (1/weight)
    outputFile = fileOutput;
}

Netshort::~Netshort() {
}

void Netshort::run() 
{ 
    VertexIterator it, itEnd;
    std::map<Vertex, float>::iterator vt, vtEnd;
    std::map<Vertex, Vertex> mapPredecessor;
    std::map<Vertex, float> mapDistance;
    float score = 0;
    //std::cout << flagAccumulateToInitialNodeScore << std::endl;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	getNetwork().calculateShortestPath(*it, mapDistance, mapPredecessor); 
	score = 0;
	//std::cout << score << std::endl;
	for(vt = mapDistance.begin(), vtEnd = mapDistance.end(); vt != vtEnd; ++vt) 
	{
	    //std::cout << vt->second << " ";
	    //score += vt->second;
	    if(vt->second != 0)
		score += 1/vt->second;
	}
	//getNetwork().setVertexScoreUpdated(*it, score);
	//score = 1 / score; // # 1000/score #! depends on the initial weights
	if(flagAccumulateToInitialNodeScore) 
	{
	    score += getNetwork().getVertexScore(*it);
	}
	getNetwork().setVertexScore(*it, score);
	//std::cout << std::endl << getNetwork().getVertexIndex(*it)  << " " << getNetwork().getVertexName(*it) << " " << 10000/score << std::endl;
    }
    getNetwork().scaleVertexScores(SCALE_BETWEEN_ZERO_AND_ONE);
    getNetwork().outputScores(outputFile);
    return;
}



/////
float Netshort::testRun() {
    std::map<Vertex, Vertex> mapPredecessor;
    std::map<Vertex, float> mapDistance;
    std::map<Vertex, float>::iterator vt, vtEnd;
    getNetwork().calculateShortestPath(getNetwork().getVertex("25604"), mapDistance, mapPredecessor); 
    //mapDistance = getNetwork().calculateShortestPath(getNetwork().getVertex("25604")); 
    float score = 0;
    for(vt = mapDistance.begin(), vtEnd = mapDistance.end(); vt != vtEnd; ++vt) 
    {
	//std::cout << getNetwork().getVertexIndex(vt->first) << std::endl;
	//std::cout << getNetwork().getVertexName(vt->first) << " " << vt->second << std::endl;
	score += vt->second;
    }
    return 1000/score;
}


