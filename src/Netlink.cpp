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

#include "Netlink.hpp"

#include <iostream>
#include <cstdlib>

Netlink::Netlink() : ScoreNetwork()
{
    //flagAccumulateToInitialNodeScore = false;
    //flagVerbose = false;
    threshold = 0.0;
}

Netlink::Netlink(std::string fileNode, std::string fileEdge, std::string fileOutput, float t, bool fAccumulateToInitialNodeScore, bool fVerbose) : ScoreNetwork(fileNode, fileEdge, fileOutput, false, fAccumulateToInitialNodeScore, false, fVerbose)
{
    //flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore;
    //flagVerbose = fVerbose;
    //getNetwork().loadNodes(fileNode); 
    //getNetwork().loadEdges(fileEdge); 
    //outputFile = fileOutput;
    threshold = t;
}

Netlink::~Netlink() {
}

void Netlink::updateNodeScore(Vertex v)
{ 
    AdjVertexIterator ut, utEnd;
    float score = 0.0;
    unsigned int i = 0;
    if(flagVerbose) {
	std::cout << "Checking node " << getNetwork().getVertexName(v) << std::endl;
    }
    for(boost::tie(ut, utEnd) = getNetwork().getAdjacentVertexIteratorOfVertex(v); ut != utEnd; ++ut) {
	if(v != *ut) { // Skip self edges
	    if(flagVerbose) {
		std::cout << "- Score from neighbor " << getNetwork().getVertexName(*ut) << ": " << getNetwork().getVertexScore(*ut) << std::endl;
	    }
	    score += getNetwork().getVertexScore(*ut);
	    i += 1;
	}
    }
    // Scale by number of neighbors
    //if(i!=0) score /= i;
    if(score >= threshold) {
	score = 1.0;
    } else {
	score = 0.0;
    }
    if(flagAccumulateToInitialNodeScore) {
	score += getNetwork().getVertexScore(v);
    }
    setVertexScoreUpdated(v, score);
    return;
}

/*
void Netlink::run(float threshold)
{ 
    VertexIterator it, itEnd;
    AdjVertexIterator ut, utEnd;
    float score = 0.0;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	if(flagVerbose) {
	    std::cout << "Checking node " << getNetwork().getVertexName(*it) << std::endl;
	}
	score = 0.0;
	for(boost::tie(ut, utEnd) = getNetwork().getAdjacentVertexIteratorOfVertex(*it); ut != utEnd; ++ut) {
	    if(*it != *ut) { // Skip self edges
		if(flagVerbose) {
		    std::cout << "- Score from neighbor " << getNetwork().getVertexName(*ut) << ": " << getNetwork().getVertexScore(*ut) << std::endl;
		}
		score += getNetwork().getVertexScore(*ut);
	    }
	}
	if(score >= threshold) {
	    score = 1.0;
	}
	if(flagAccumulateToInitialNodeScore) {
	    score += getNetwork().getVertexScore(*it);
	}
	getNetwork().setVertexScore(*it, score);
    }
    getNetwork().outputScores(outputFile);
    return;
}
*/

