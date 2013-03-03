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

#include "ScoreNetwork.hpp"
#include <iostream>
//#include <cstdlib> 
//#include <ctime>
//#include <sstream>

using namespace std;

ScoreNetwork::ScoreNetwork() 
{
    typeTransfer = IDENTITY;
    //minScore = INFINITY;
    //maxScore = -INFINITY;
    //nIteration = 0;
    repeatCounter = 0;
    iterationCounter = 0;
    flagVerbose = false;
    flagUseEdgeScore = true;
    flagAccumulateToInitialNodeScore = true;
    flagResetSeedScoresToInitial = false;
    outputFile = "";
    pNetwork = 0;
}

ScoreNetwork::ScoreNetwork(std::string fileOutput, bool fUseEdgeScore, bool fAccumulateToInitialNodeScore, bool fResetSeedScoresToInitial, bool fVerbose) 
{
    typeTransfer = IDENTITY;
    //minScore = INFINITY;
    //maxScore = -INFINITY;
    //nIteration = 0;
    repeatCounter = 0;
    iterationCounter = 0;
    flagVerbose = fVerbose;
    flagUseEdgeScore = fUseEdgeScore; // should? be in updateEdgeScore
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore; // should? be in updateNodeScore
    flagResetSeedScoresToInitial = fResetSeedScoresToInitial; 
    outputFile = fileOutput;
    pNetwork = 0;
}

ScoreNetwork::ScoreNetwork(string fileNode, string fileEdge, string fileOutput, bool fUseEdgeScore, bool fAccumulateToInitialNodeScore, bool fResetSeedScoresToInitial, bool fVerbose) 
{
    typeTransfer = IDENTITY;
    //minScore = INFINITY;
    //maxScore = -INFINITY;
    //nIteration = 0;
    repeatCounter = 0;
    iterationCounter = 0;
    flagVerbose = fVerbose;
    flagUseEdgeScore = fUseEdgeScore; // should? be in updateEdgeScore
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore; // should? be in updateNodeScore
    flagResetSeedScoresToInitial = fResetSeedScoresToInitial; 
    outputFile = fileOutput;
    pNetwork = new Graph();
    getNetwork().loadNodes(fileNode); 
    getNetwork().loadEdges(fileEdge); 
    //cout << "Snet file:" << fileOutput << " " << outputFile << endl;
    //cout << "Snet fverbose:" << fVerbose << " " << flagVerbose << endl;
    //cout << "Snet fAccumulate: " << this->flagAccumulateToInitialNodeScore << endl;
    //boost::tie(minScore, maxScore) = getMinAndMaxNodeScores();
}

ScoreNetwork::~ScoreNetwork() {
}

void ScoreNetwork::run(int nRepetition, int nIteration, float tError) 
{ 
    float error = INFINITY;
    initializeScoring();
    //cout << "Snet fAccumulate in run: " << flagAccumulateToInitialNodeScore << endl;
    for(repeatCounter = 1; repeatCounter<=nRepetition and error > tError; ++repeatCounter) {
	if(flagVerbose)
	    cout << "((repetition " << repeatCounter << "))" << std::endl;
	initializeRepetition();
	for(iterationCounter = 1; iterationCounter<=nIteration and error > tError; ++iterationCounter) {
	    if(flagVerbose)
		cout << "(iteration " << iterationCounter << ")" << std::endl;
	    initializeIteration();
	    updateNetwork();
	    error = calculateErrorAndUpdateScores();
	    finalizeIteration();
	} 
	finalizeRepetition();
    }
    finalizeScoring();
    //cout << "Snet run file: " << this->outputFile << endl;
    getNetwork().outputScores(outputFile);
}

/*
// check which is called
void ScoreNetwork::initializeIteration() {
    cout << endl << "initialize of base" << endl;  
}

//virtual void ScoreNetwork::finalizeIteration() { }

// should it be in destructor
void ScoreNetwork::finalizeScoring() {
    scaleNodeScores();
}
*/
/*
void ScoreNetwork::updateNodeScore(Vertex v) {}
void ScoreNetwork::updateEdgeScore(Edge e) {}
*/

void ScoreNetwork::updateNetwork() 
{ 
    VertexIterator it, itEnd;
    //OutEdgeIterator et, etEnd;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
        updateNodeScore(*it); 
	/*
	for(boost::tie(et, etEnd) = getNetwork().getEdgeIteratorOfVertex(*it); et != etEnd; ++et) 
	{
	    // Make sure that we update edges only once since it is a an undirected graph
	    if(getNetwork().getVertexIndex(*it) >= getNetwork().getVertexIndex(getNetwork().getTargetVertex(*et))) {
		updateEdgeScore(*et);
	    }
	}
	*/
    }
    //if(flagVerbose)	
    //        printNetwork("After updateNetwork");
}

float ScoreNetwork::calculateErrorAndUpdateScores() {
    float error = 0.0, sumErrorNode = 0.0; //, sumErrorEdge = 0.0;
    VertexIterator it, itEnd;
    OutEdgeIterator et, etEnd;
    for(boost::tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
        if(isnan(getNetwork().getVertexScore(*it)) or isinf(getNetwork().getVertexScore(*it))) {
	   cerr << "NAN or INF score" << endl;
	}
        error = getVertexScoreUpdated(*it) - getNetwork().getVertexScore(*it); 
        sumErrorNode += sq(error);
	getNetwork().setVertexScore(*it, getVertexScoreUpdated(*it)); 
	// For the time being, not considering updation of edge scores
	/*
	for(boost::tie(et, etEnd) = getNetwork().getEdgeIteratorOfVertex(*it); et != etEnd; ++et) 
	{
	    // Make sure that we update edges only once since it is a an undirected graph
	    if(getNetwork().getVertexIndex(*it) >= getNetwork().getVertexIndex(getNetwork().getTargetVertex(*et))) {
                error = getNetwork().getEdgeScore(*et);
		sumErrorEdge += sq(error);
		//!pEdge->data.score = pEdge->data.scoreUpdated;	            
	    }
	}
	*/
    }
    //cerr << "minScore: " << minScore << " maxScore: " << maxScore << endl;
    //return max(float(sqrt( (sumErrorEdge/nEdge)*(sumErrorNode/nNode) )), float( (sqrt(sumErrorNode/nNode)+ sqrt(sumErrorEdge/nEdge)) /2)); //float(sqrt(sumErrorEdge/nEdge)); //float(sqrt(sumErrorNode/nNode)+sqrt(sumErrorEdge/nEdge))/2;
    return float(sqrt( sumErrorNode/getNetwork().getSize() ));
}

float ScoreNetwork::transferScore(float score, float a, float b) {
    switch(typeTransfer) {
    case IDENTITY:
        return score;
    case POLYNOMIAL:
        return a*score+b; // a=0.5, b=0.1, 0.2 - 0.6
    case LOGARITHMIC:
        return float(-(1/10)*log(score)); // 0.0 - 1.0 "J"
    case EXPONENTIAL:
        return float(exp(0.75*score)/10); // 0.1 - 0.2 //float(1-exp(-0.75*score)); // 0.1 - 0.5
    default:
        //cout << "Warning: unidentified type" << endl;
	throw TypeException(); //typeEx;
        return 0.0;
    }
}


