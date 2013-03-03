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

#include "Netween.hpp"

#include <iostream>
#include <vector>

using namespace boost;
using namespace std;

Netween::Netween() 
{
    flagAccumulateToInitialNodeScore = false;
    flagVerbose = false;
}

Netween::Netween(string fileNode, string fileEdge, string fileOutput, float seedScoreThreshold, bool fAccumulateToInitialNodeScore, bool fVerbose)
{
    flagAccumulateToInitialNodeScore = fAccumulateToInitialNodeScore;
    flagVerbose = fVerbose;
    outputFile = fileOutput;
    network.loadNodes(fileNode); 
    network.loadEdges(fileEdge);//, true); 
    VertexIterator it, itEnd;
    for(tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) 
    {
	if(getNetwork().getVertexScore(*it) > seedScoreThreshold)
	    setSeed.insert(*it);
    }
}


Netween::~Netween() {
}


void Netween::run() { 

    // Initialize node involvement counts
    initializeNodeCounts();

    // Get nodes included in the possible shortest paths of seed pairs 
    getIncludedNodes();

    // Get node involvement counts
    getNodeCounts();

    // Calculate scores
    calculateScores();

    // Output scores
    getNetwork().scaleVertexScores(SCALE_BETWEEN_ZERO_AND_ONE);
    getNetwork().outputScores(outputFile);

    return;
}


void Netween::initializeNodeCounts() {
    VertexIterator it, itEnd;
    for(tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) {
	mapVertexToLocalCount[*it] = 0;
	mapVertexToGlobalCount[*it] = 0;
    }
}


void Netween::getIncludedNodes() {
    unordered_set<Vertex>::iterator seedIt, seedItEnd;
    unordered_set<Vertex>::iterator setIt, setItEnd;
    Vertex v, v_prev;
    PredecessorList * pMapPredecessors;

    // For each seed find shortest paths: Ps,t (s € seeds, t € nodes)
    for(seedIt = setSeed.begin(), seedItEnd = setSeed.end(); seedIt != seedItEnd; ++seedIt) {
	pMapPredecessors = new PredecessorList();
	if(flagVerbose)
	    cout << "Checking seed " << getNetwork().getVertexName(*seedIt) << endl;
	getNetwork().getAllShortestPaths(*seedIt, *pMapPredecessors); 
	vertexToMapPredecessors[*seedIt] = pMapPredecessors;
	// Record nodes (I) involved in Ps,s' (s, s' € seeds) where P denotes shortest path
	for(setIt = setSeed.begin(), setItEnd = setSeed.end(); setIt != setItEnd; ++setIt) {
	    v = *setIt;
	    // Skip the node if it is the seed itself
	    if(v == *seedIt)
		continue;
	    v_prev = (*pMapPredecessors)[v][0];
	    //cout << "v " << getNetwork().getVertexName(v) << ", v_prev " << getNetwork().getVertexName(v_prev) << endl;
	    // Add all nodes in the (all possible) shortest path(s) Ps,s'
	    if(flagVerbose)
		cout << "- Including " << getNetwork().getVertexName(v) << endl;
	    setIncluded.insert(v);
	    while(v_prev != v) { 
		// Need not to continue checking if initial seed node is reached
		if(v_prev == *seedIt) 
		    break;
		for(unsigned int i=0; i < (*pMapPredecessors)[v].size(); ++i) {
		    if(flagVerbose)
			cout << "- Including " << getNetwork().getVertexName((*pMapPredecessors)[v][i]) << endl;
		    setIncluded.insert((*pMapPredecessors)[v][i]);
		}
		v_prev = v;
		v = (*pMapPredecessors)[v][0];
	    }
	}
    }

    return;
}
	

void Netween::getNodeCounts() {
    unordered_set<Vertex>::iterator seedIt, seedItEnd;
    unordered_set<Vertex>::iterator setIt, setItEnd;
    VertexIterator it, itEnd;
    Vertex vTarget;
    PredecessorList * pMapPredecessors;

    // For each i € I, check how many times i is involved in all possible Ps,t (s € seeds, t € nodes) 
    // Count i only once for all possible distinct Ps,t for a given pair s,t
    for(setIt = setIncluded.begin(), setItEnd = setIncluded.end(); setIt != setItEnd; ++setIt) {
	if(flagVerbose)
	    cout << "Checking included node " << getNetwork().getVertexName(*setIt) << endl;
	for(seedIt = setSeed.begin(), seedItEnd = setSeed.end(); seedIt != seedItEnd; ++seedIt) {
	    pMapPredecessors = vertexToMapPredecessors[*seedIt];
	    if(flagVerbose)
		cout << "- Checking seed " << getNetwork().getVertexName(*seedIt) << endl;
	    for(tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) {
		vTarget = *it;
		if(vTarget == *seedIt || vTarget == *setIt)
		    continue;
		if(flagVerbose)
		    cout << "- Checking shortest path to " << getNetwork().getVertexName(vTarget) << endl;
		// If included is seed, count the seed node on the path of all targets in inclusion analysis
		if(*seedIt == *setIt) {
		    //if(flagVerbose)
		    //	cout << "-- Counting " << getNetwork().getVertexName(*setIt) << endl;
		    //if(setSeed.find(vTarget) != seedItEnd) { 
		    //	mapVertexToLocalCount[*setIt] += 2; // this will not be counted again because of the constraint of vTarget != *setIt
		    //}
		    //mapVertexToGlobalCount[*setIt] += 1;
		}
		else {
		    if(isVertexIncludedInsideThePathOfGivenVertex(*setIt, vTarget, *pMapPredecessors)) {
			if(flagVerbose)
			    cout << "-- Counting " << getNetwork().getVertexName(*setIt) << endl;
			// For every pair s,s' i is counted twice (this number is diveded by two later in score calculation)
			if(setSeed.find(vTarget) != seedItEnd) {
			    mapVertexToLocalCount[*setIt] += 1;
			}
			mapVertexToGlobalCount[*setIt] += 1;
		    }
		}
	    }
	}
    }

    return;
}


bool Netween::isVertexIncludedInsideThePathOfGivenVertex(Vertex vToBeChecked, Vertex vTarget, PredecessorList & mapPredecessors) {
    Vertex v=vTarget, v_prev;
    v_prev = mapPredecessors[v][0];
    // The target vertex is not checked intentionaly, trying to see if it is on the path to the target
    while(v_prev != v) { 
	for(unsigned int i=0; i < mapPredecessors[v].size(); ++i) {
	    if(vToBeChecked == mapPredecessors[v][i])
		return true;
	}
	v_prev = v;
	v = mapPredecessors[v][0];
    }
    return false;
}

    /*
	// Below code fails in the following setting (was called after each seed's shortest path calculation with newly included nodes):
	// Let x a node included in a shortest path of sw,sv but not in su,sv and su,sw
	// then global counts of x when x is on su, t may be ignored because at the time of checking su it was not in setIncluded
	for(tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) {
	    v = *it;
	    // Skip the check if it is the seed node itself or the node is not in the included set
	    if(v == *seedIt || setIncluded.find(v) != setIncluded.end())
		continue;
	    if(flagVerbose)
		cout << "- Checking shortest path to " << getNetwork().getVertexName(v) << endl;
	    v_prev = mapPredecessors[v][0];
	    while(v_prev != v) { 
		for(unsigned int i=0; i < mapPredecessors[v].size(); ++i) {
		    if(setIncluded.find(v) != setIncluded.end()) {
			if(flagVerbose)
			    cout << "-- Counting " << getNetwork().getVertexName(v) << endl;
			if(setSeed.find(*it) != seedItEnd) {
			    mapVertexToLocalCount[v] += 1;
			}
			mapVertexToGlobalCount[v] += 1;
		    }
		}
		v_prev = v;
		v = mapPredecessors[v][0];
	    }
	}
	setIncluded.clear();
	//mapPredecessors.clear();
    */


void Netween::calculateScores() {
    VertexIterator it, itEnd;
    float score = 0;

    for(tie(it, itEnd) = getNetwork().getVertexIterator(); it != itEnd; ++it) {
	// Divide local scores by 2 since a node should be counted only once for each pair s,s'
	mapVertexToLocalCount[*it] /= 2;
	// Subtract duplicated local score contributions from global scores
       	mapVertexToGlobalCount[*it] -= mapVertexToLocalCount[*it]; 
	if(flagVerbose)
	    cout << getNetwork().getVertexName(*it)  << ": " << mapVertexToLocalCount[*it] << " " << mapVertexToGlobalCount[*it] << endl;
	if(mapVertexToLocalCount[*it] == 0) 
	    score = 0;
	else
	    score = float(mapVertexToLocalCount[*it]) / mapVertexToGlobalCount[*it]; 
	if(flagAccumulateToInitialNodeScore) {
	    score += getNetwork().getVertexScore(*it);
	}
	getNetwork().setVertexScore(*it, score);
    }

    return;
}



