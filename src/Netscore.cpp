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

#include "Netscore.hpp"

#include <iostream>

Netscore::Netscore() : ScoreNetwork()
{
}

Netscore::Netscore(std::string fileOutput, bool fUseEdgeScore, bool fAccumulateToInitialNodeScore, bool fResetSeedScoresToInitial, bool fVerbose) : ScoreNetwork(fileOutput, fUseEdgeScore, fAccumulateToInitialNodeScore, fResetSeedScoresToInitial, fVerbose)
{
}

Netscore::Netscore(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fUseEdgeScore, bool fAccumulateToInitialNodeScore, bool fResetSeedScoresToInitial, bool fVerbose) : ScoreNetwork(fileNode, fileEdge, fileOutput, fUseEdgeScore, fAccumulateToInitialNodeScore, fResetSeedScoresToInitial, fVerbose)
{
    //std::cout << "Nscore fAccumulate: " << flagAccumulateToInitialNodeScore << std::endl;
    //std::cout << "Nscore file: " << fileOutput << std::endl;
    //if base constructor is not called, need to assign all like 
    //outputFile = fileOutput;
    //flagVerbose = flagVerbose;
    //this->flagAccumulateToInitialNodeScore = flagAccumulateToInitialNodeScore; 
    // Problem caused by pythonic assignment of instance variable with the same name as argument
}

Netscore::~Netscore()
{
    VertexToMessageMap::iterator it, itEnd;
    for(it = messageMaps.begin(), itEnd = messageMaps.end(); it != itEnd; ++it) 
    {
	delete it->second;
    }
}

void Netscore::initializeScoring() {
    scaleNodeScores(SCALE_BETWEEN_ZERO_AND_ONE);
    VertexIterator it, itEnd;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	createVertexMessageMap(*it);
	setVertexScoreUpdated(*it, 0.0);
	// Moved below to initializeRepetition
	//setVertexScoreInitial(*it, getNetwork().getVertexScore(*it)); // storing scaled (between 0 and 1) initial scores
	//getNetwork().setVertexScore(*it, 0.0);
	setVertexScoreInitialOriginal(*it, getNetwork().getVertexScore(*it)); // storing scaled (between 0 and 1) initial scores
    }
}

void Netscore::finalizeScoring()
{
    VertexIterator it, itEnd;
    // Add initial score if accumulation to initial score is desired 
    if(flagAccumulateToInitialNodeScore == true) 
    {
	for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
	{
	    getNetwork().setVertexScore(*it, getNetwork().getVertexScore(*it) + getVertexScoreInitialOriginal(*it)); // before was adding only initial scores from the last repetition
	}
    }
    scaleNodeScores(SCALE_BETWEEN_INITIAL_MIN_AND_MAX_SCORE);
}

void Netscore::initializeRepetition()
{
    VertexIterator it, itEnd;
    UIntToMessage * pMapMessage;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	//std::cout << getNetwork().getVertexName(*it) << std::endl;
	//mapMessage = getNetwork().getVertexMessageMap(*it); 
	//mapMessage[getNetwork().getVertexIndex(*it)] = std::pair<float, int>(1.0, 0);
	pMapMessage = getVertexMessageMap(*it); 
	if(repeatCounter != 1) {
	    pMapMessage->clear();
	}
	(*pMapMessage)[getNetwork().getVertexIndex(*it)] = boost::make_tuple(1.0, 0, 1); //std::pair<float, int>(1.0, 0);
	//UIntToMessage::iterator mt, mtEnd;
	//for(mt=mapMessage.begin(), mtEnd=mapMessage.end(); mt != mtEnd; ++mt)
	//    std::cout << mt->first << " " << mt->second.second << std::endl;
	//mapMessage = getNetwork().getVertexMessageMap(*it); 
	//for(mt=mapMessage.begin(), mtEnd=mapMessage.end(); mt != mtEnd; ++mt)
	//    std::cout << mt->first << " " << mt->second.second << std::endl;
	setVertexScoreInitial(*it, getNetwork().getVertexScore(*it)); // storing scaled (between 0 and 1) scores from the last repetition's iteration or initial scores 
	//setVertexScoreUpdated(*it, 0.0); // moved back to initializeScoring()
	getNetwork().setVertexScore(*it, 0.0);
    }
}

void Netscore::finalizeIteration()
{
    scaleNodeScores(SCALE_BETWEEN_ZERO_AND_ONE);
}

void Netscore::updateNodeScore(Vertex v) 
{
    AdjVertexIterator vt, vtEnd;
    UIntToMessage::iterator it, itEnd, itSearch;
    UIntToMessage * pMapMessage = getVertexMessageMap(v), mapMessageNeighbor;
    float tempScore = 0.0;
    //boost::unordered_map<unsigned int, unsigned int> mapIterationToCount;
    unsigned int nMessagesRecieved = 0;
    //std::cout << "Nscore fverbose : " << flagVerbose << " fAccumulate: " << flagAccumulateToInitialNodeScore << std::endl; //<< " T F: " << true << " " << false <<std::endl;
    if(flagVerbose)
        std::cout << "-Checking node: " << getNetwork().getVertexName(v) << std::endl; //<<  " (iteration " << iterationCounter << ")" << std::endl;
    for(boost::tie(vt, vtEnd) = getNetwork().getAdjacentVertexIteratorOfVertex(v); vt != vtEnd; ++vt) 
    {
	if(v == *vt) { // Skip self edges
	    continue;
	}
	if(flagVerbose)
	    std::cout << "--Messages from neighbor: " << getNetwork().getVertexName(*vt) << std::endl; 
	mapMessageNeighbor = *(getVertexMessageMap(*vt)); 
	for(it=mapMessageNeighbor.begin(), itEnd=mapMessageNeighbor.end(); it != itEnd; ++it)
	{
	    //std::cout << boost::get<1>(it->second) << " " << iterationCounter << std::endl;
	    // Consider the message if only coming from a previous iteration
	    if(boost::get<1>(it->second) < iterationCounter)
	    {
		//tempScore = it->second.first;
		tempScore = boost::get<0>(it->second);
		if(flagUseEdgeScore) 
		{
		    tempScore *= getNetwork().getEdgeScore(v, *vt);
		}
		// Sum up messages from the same node arriving at the same time
		itSearch = pMapMessage->find(it->first);
		if(itSearch == pMapMessage->end())
		{
		    (*pMapMessage)[it->first] = boost::make_tuple(tempScore, iterationCounter, 1); //std::pair<float, int>(tempScore, iterationCounter); // iterationCounter == it->second.second + 1
		    if(flagVerbose)
			std::cout << "---message of " << getNetwork().getVertexName(it->first) << " scoring factor: " << tempScore << std::endl; 
		} 
		else 
		{
		    if(boost::get<1>(itSearch->second) == iterationCounter) 
		    {
			boost::get<0>(itSearch->second) += tempScore;
			boost::get<2>(itSearch->second) += 1;
			if(flagVerbose)
			    std::cout << "---message of " << getNetwork().getVertexName(it->first) << " scoring factor: " << tempScore << std::endl; 
		    }
		}
	    }
	}
    }
    tempScore = 0.0;
    nMessagesRecieved = 0;
    for(it=pMapMessage->begin(), itEnd=pMapMessage->end(); it != itEnd; ++it)
    {
	// Vertex score contains accumulated score over the iterations so only consider messages from that iteration
	// Initial scores of all nodes are used while calculating node score based on the messages in the message array
	if(boost::get<1>(it->second) == iterationCounter) 
	{
	    tempScore += boost::get<0>(it->second) * getVertexScoreInitial(getNetwork().getVertex(it->first));
	    // If scaling by iteration is desired
	    //if(iterationCounter != 0)
	    //	tempScore /= iterationCounter; 

	    //nMessagesRecieved += 1;
	    nMessagesRecieved += boost::get<2>(it->second);
	    if(flagVerbose)
		std::cout << "--Evaluating message of " << getNetwork().getVertexName(it->first) << ": accumulated score: " << tempScore << " (previous + " << boost::get<0>(it->second) << " * " << getVertexScoreInitial(getNetwork().getVertex(it->first)) <<" )" << std::endl; 
	}
    }
    if(nMessagesRecieved != 0) 
    {
        tempScore /= nMessagesRecieved; 
    }
    tempScore += getNetwork().getVertexScore(v);
    // Remove initial score if accumulation to initial score is not desired
    //if(flagAccumulateToInitialNodeScore == false && iterationCounter == 1) 
    //{
    //	tempScore -= getVertexScoreInitial(v);
    //}
    // Update scoreUpdated (for error calculation)
    setVertexScoreUpdated(v, tempScore);
    if(flagVerbose)
	std::cout << "Score calculated for " << getNetwork().getVertexName(v) << ": " << tempScore << std::endl; //<< " (iteration " << iterationCounter << ")" << std::endl;
}

