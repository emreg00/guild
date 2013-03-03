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

#include "Netzcore.hpp"

#include <iostream>
#include <sstream>

Netzcore::Netzcore() : ScoreNetwork()
{
    nSampledGraphs = 0;
    mean = 0;
    sigma = 0;
}

Netzcore::Netzcore(std::string fileNode, std::string fileOutput, std::string prefixSampled, unsigned int nSampled, bool fUseEdgeScore, bool fAccumulateToInitialNodeScore, bool fResetSeedScoresToInitial, bool fVerbose) : ScoreNetwork(fileOutput, fUseEdgeScore, fAccumulateToInitialNodeScore, fResetSeedScoresToInitial, fVerbose)
{ 
    prefixSampledGraphs = prefixSampled;
    nodeFileName = fileNode;
    nSampledGraphs = nSampled;
    mean = 0;
    sigma = 0;
}

Netzcore::Netzcore(std::string fileNode, std::string fileEdge, std::string fileOutput, std::string prefixSampled, unsigned int nSampled, bool fUseEdgeScore, bool fAccumulateToInitialNodeScore, bool fResetSeedScoresToInitial, bool fVerbose) : ScoreNetwork(fileNode, fileEdge, fileOutput, fUseEdgeScore, fAccumulateToInitialNodeScore, fResetSeedScoresToInitial, fVerbose)
{ 
    prefixSampledGraphs = prefixSampled;
    nodeFileName = fileNode;
    nSampledGraphs = nSampled;
    mean = 0;
    sigma = 0;
}

Netzcore::~Netzcore()
{
    std::list<Graph*>::iterator it, itEnd;
    for(it=sampledGraphs.begin(), itEnd=sampledGraphs.end(); it != itEnd; ++it)
    {
	delete *it;
    }
}

void Netzcore::loadSampledGraphs()
{
    std::ostringstream oss;
    for(unsigned int i = 1; i <= nSampledGraphs; i++)
    {	
	oss.str("");
	oss << prefixSampledGraphs << i;
	//std::cout << oss.str() << std::endl;
	Graph* pG = new Graph();
	pG->loadNodes(nodeFileName);
	pG->loadEdges(oss.str());
	sampledGraphs.push_back(pG);
    }
}

void Netzcore::printSampledGraphs()
{
    std::list<Graph*>::iterator it, itEnd;
    Graph* pG;
    for(it=sampledGraphs.begin(), itEnd=sampledGraphs.end(); it != itEnd; ++it)
    {	
	pG = *it;
	pG->print(true);
    }
}

void Netzcore::updateSampledGraphScores()
{
    std::list<Graph*>::iterator it, itEnd;
    Graph* pG;
    VertexIterator vt, vtEnd;
    for(it=sampledGraphs.begin(), itEnd=sampledGraphs.end(); it != itEnd; ++it)
    {	
	pG = *it;
	for(boost::tie(vt, vtEnd) = getNetwork().getVertexIterator(); vt != vtEnd; ++vt) {
	    pG->setVertexScore(pG->getVertex(getNetwork().getVertexName(*vt)), getNetwork().getVertexScore(*vt));
	}
    }
}

void Netzcore::initializeScoring() {
    VertexIterator it, itEnd;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	//createVertexMessageMap(*it);
	setVertexScoreInitial(*it, getNetwork().getVertexScore(*it));
	setVertexScoreUpdated(*it, 0.0);
    }
    loadSampledGraphs();
    //printSampledGraphs();
    //float vertex_scores[] = { 0.0, 1.0, 0.0, 1.0 };
    //std::pair<float, float> tempPair = calculateMeanAndSigma(vertex_scores, vertex_scores+sizeof(vertex_scores)/sizeof(float));
    //std::cout << tempPair.first << " " << tempPair.second << std::endl;
}

void Netzcore::finalizeScoring()
{
    VertexIterator it, itEnd;
    // Add initial score if accumulation to initial score is desired 
    if(flagAccumulateToInitialNodeScore == true) 
    {
	for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
	{
	    getNetwork().setVertexScore(*it, getNetwork().getVertexScore(*it) + getVertexScoreInitial(*it));
	}
    }
    scaleNodeScores(SCALE_BETWEEN_INITIAL_MIN_AND_MAX_SCORE);
}

void Netzcore::initializeIteration()
{
    std::list<Graph*>::iterator it, itEnd;
    Graph* pG;
    Vertex u;
    AdjVertexIterator vt, vtEnd;
    VertexIterator t, tEnd;
    float tempScore = 0.0, sumScore = 0.0;
    unsigned int i = 0;

    std::list<float> scores;
    for(it=sampledGraphs.begin(), itEnd=sampledGraphs.end(); it != itEnd; ++it) {	
	pG = *it;
	for(boost::tie(t, tEnd) = pG->getVertexIterator(); t != tEnd; ++t) {
	    u = *t;
	    if(flagVerbose)
		std::cout << "Checking neighbors of node: " << pG->getVertexName(u) << std::endl;
	    sumScore = 0.0;
	    tempScore = 0.0;
	    i = 0;
	    bool firstTime = true;
	    for(boost::tie(vt, vtEnd) = pG->getAdjacentVertexIteratorOfVertex(u); vt != vtEnd; ++vt) 
	    {
		if(u == *vt) { // Skip (one of the) self edges
		    continue; // skip self edges
		    if(!firstTime) {
			continue;
		    } 
		    else {
			firstTime = false;
		    }
		}
		tempScore = pG->getVertexScore(*vt);
		if(flagUseEdgeScore) 
		{
		    tempScore *= pG->getEdgeScore(u, *vt);
		}
		sumScore += tempScore;
		i += 1;
		if(flagVerbose)
		    std::cout << "--Score of random neighbor " << pG->getVertexName(*vt) << ": " << tempScore << std::endl; 
	    }
	    if(i == 0) { // Disconnected node 
		tempScore = sumScore;
	    } else {
		tempScore = sumScore / i;
	    }
	    if(flagVerbose)
		std::cout << "-Average score of random neighbors " << tempScore << std::endl; 
	    scores.push_back(tempScore);
	}
    }
    boost::tie(this->mean, this->sigma) = calculateMeanAndSigma(scores.begin(), scores.end());
    if(flagVerbose)
	std::cout << "Mean: " << this->mean << " sigma: " << this->sigma << std::endl; 
}

void Netzcore::finalizeIteration()
{
    scaleNodeScores(SCALE_BETWEEN_ZERO_AND_ONE);
    updateSampledGraphScores();
}

void Netzcore::updateNodeScore(Vertex v) 
{
    AdjVertexIterator vt, vtEnd;
    std::list<Graph*>::iterator it, itEnd;
    float tempScore = 0.0, sumScore = 0.0;
    std::pair<float, float> tempPair(0.0,0.0);
    unsigned int i = 0;
    //Vertex u;
    //Graph* pG;
    if(flagVerbose)
        std::cout << "Checking node: " << getNetwork().getVertexName(v) << std::endl;

    // Before local mean & sigma was considered
    /*
    std::list<float> scores;
    for(it=sampledGraphs.begin(), itEnd=sampledGraphs.end(); it != itEnd; ++it)
    {	
	pG = *it;
	u = pG->getVertex(getNetwork().getVertexName(v));
	//pG->print(false);
	sumScore = 0.0;
	tempScore = 0.0;
	i = 0;
	for(boost::tie(vt, vtEnd) = pG->getAdjacentVertexIteratorOfVertex(u); vt != vtEnd; ++vt) 
	{
	    if(u == *vt) { // Skip self edges
		continue;
	    }
	    tempScore = pG->getVertexScore(*vt);
            if(flagUseEdgeScore) 
	    {
		tempScore *= pG->getEdgeScore(u, *vt);
            }
	    sumScore += tempScore;
	    i += 1;
	    if(flagVerbose)
		std::cout << "--Score of random neighbor " << pG->getVertexName(*vt) << ": " << tempScore << std::endl; 
	}
	tempScore = sumScore / i;
	if(flagVerbose)
	    std::cout << "-Average score of random neighbors " << tempScore << std::endl; 
	scores.push_back(tempScore);
    }
    tempPair = calculateMeanAndSigma(scores.begin(), scores.end());
    */
    // Now considering global mean and sigma
    tempPair = std::make_pair(this->mean, this->sigma);

    sumScore = 0.0;
    tempScore = 0.0;
    i = 0;
    bool firstTime = true;
    for(boost::tie(vt, vtEnd) = getNetwork().getAdjacentVertexIteratorOfVertex(v); vt != vtEnd; ++vt) 
    {
	if(v == *vt) { // Skip (one of the) self edges
	    continue; // skip self edges
	    if(!firstTime) {
		continue;
	    } 
	    else {
		firstTime = false;
	    }
	}
	tempScore = getNetwork().getVertexScore(*vt);
	if(flagUseEdgeScore) 
	{
	    tempScore *= getNetwork().getEdgeScore(v, *vt);
	} 
	sumScore += tempScore;
	i += 1;
	if(flagVerbose)
	    std::cout << "--Score of neighbor " << getNetwork().getVertexName(*vt) << ": " << tempScore << std::endl; 
    }
    if(i == 0) { // Disconnected node 
	tempScore = sumScore;
    } else {
	tempScore = sumScore / i;
    }
    if(flagVerbose)
	std::cout << "-Average score of neighbors " << tempScore << std::endl; 
    //std::cout << "s: " << tempScore << " m: " << tempPair.first << " sig: " << tempPair.second << std::endl;
    tempScore = tempScore - tempPair.first;
    if(tempScore == 0 && tempPair.second == 0) {
	std::cout << "0/0 due to Zero variance!" << std::endl;	
	tempScore = 0.0;
    } else {
	tempScore = tempScore / tempPair.second;
    }
    /*
    if(tempPair.second == 0) //(isnan(tempScore))
    {
	std::cout << "Zero variance!" << std::endl;	
	//tempScore = getNetwork().getVertexScore(v);
	if(tempScore == 0) {
	    tempScore = 0.0;
	} else {
	    tempScore = tempScore / tempPair.second;
	}
    } 
    else 
    {
	tempScore = tempScore / tempPair.second;
    }
    */
    // Update scoreUpdated (for error calculation)
    setVertexScoreUpdated(v, tempScore);
    if(flagVerbose)
	std::cout << "Score calculated for " << getNetwork().getVertexName(v) << ": " << tempScore << std::endl; // << " after iteration " << iterationCounter << std::endl;
}

template <class InputIterator>
std::pair<float, float> calculateMeanAndSigma(InputIterator first, InputIterator beyond)
{
    float sum = 0.0, squareSum = 0.0, mean=0.0;
    unsigned int count = 0;
    InputIterator first2(first);
    while(first != beyond)
    {
	//std::cout << *first << ","; //std::endl;
	sum += *first;
	squareSum += sq(*first);
	count++;
	++first;
    }
    mean = sum/count;
    // Without Bessel's correction
    //return std::make_pair(mean, sqrt(squareSum/count-sq(mean)));
    // Calculate sample (unbiased) standard deviation 
    first = first2;
    sum = 0.0;
    while(first != beyond)
    {
	sum += sq(*first-mean);
	++first;
    }
    return std::make_pair(mean, sqrt(sum/(count-1)));
}


