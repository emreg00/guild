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

#include "Netrandom.hpp"

#include <iostream>

#include <ctime>
#include <cstdlib>

Netrandom::Netrandom() 
{
    flagAccumulateToInitialNodeScore = false;
}

Netrandom::Netrandom(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fAccumulateToInitialNodeScore)
{
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore;
    network.loadNodes(fileNode); 
    network.loadEdges(fileEdge); 
    outputFile = fileOutput;
}

Netrandom::~Netrandom() {
}

void Netrandom::run() 
{ 
    VertexIterator it, itEnd;
    srand ( time(NULL) );
    float score = 0.0;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	score = 0.0;
	if(flagAccumulateToInitialNodeScore) {
	    score = getNetwork().getVertexScore(*it);
	}
	score += float(rand())/RAND_MAX;
	getNetwork().setVertexScore(*it, score);
    }
    getNetwork().outputScores(outputFile);
    return;
}


