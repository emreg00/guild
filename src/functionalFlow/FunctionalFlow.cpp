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

#include "FunctionalFlow.h"
#include <fstream>
#include <cstdlib> 
#include <ctime>
#include <vector>

FunctionalFlow::FunctionalFlow() {
}

FunctionalFlow::FunctionalFlow(string fileProtein, string fileInteraction, string fOutput, float tSeedScore, bool fRemoveSeedEffect, bool fDebug) {
    flagRemoveSeedEffect = fRemoveSeedEffect;
    flagDebug = fDebug;
    seedScoreThreshold = tSeedScore;
    loadProteins(fileProtein); 
    loadInteractions(fileInteraction); 
    fileOutput = fOutput;
    //if(flagDebug)
    //	printNetwork("Initial: ");
}

FunctionalFlow::~FunctionalFlow() {
}

void FunctionalFlow::loadProteins(string const &fileName) {
    fstream file;
    file.open(fileName.c_str(), ios::in);
    string name, dummy;
    float reserve = 0.0;
    if(!file) cerr << "Warning: file can not be opened" << endl; 
    while(file >> name >> reserve) {
        if(reserve > seedScoreThreshold) {
	    if(flagRemoveSeedEffect) {
		mapSource[addVertex("*"+name, Data(INFINITY, 0))] = NULL;
		addVertex(name, Data(0, 0));
	    } else {
		mapSource[addVertex(name, Data(INFINITY, 0))] = NULL;
	    }
    	} else {
    		addVertex(name, Data(reserve, 0));
    	}
    }
    file.close();
}

void FunctionalFlow::loadInteractions(string const &fileName) {
    fstream file;
    file.open(fileName.c_str(), ios::in);
    string name1, name2;
    float capacity = 0.0;
    if(!file) cerr << "Warning: file can not be opened" << endl;
    std::tr1::unordered_map<int, void*>::iterator it1, it2;
    //hash_map<int, void*>::iterator it1, it2;
    while(file >> name1 >> capacity >> name2) {
	if(flagRemoveSeedEffect) {
	    if(mapNameToId.find(("*"+name1).c_str()) == mapNameToId.end())
		it1 = mapSource.end();
	    else
		it1 = mapSource.find(mapNameToId.find(("*"+name1).c_str())->second);
	    if(mapNameToId.find(("*"+name2).c_str()) == mapNameToId.end())
		it2 = mapSource.end();
	    else
		it2 = mapSource.find(mapNameToId.find(("*"+name2).c_str())->second);
	    if(it1 == mapSource.end() && it1 == mapSource.end()) {
		addEdge(name1, name2, Data(0, capacity));
	    }
	    if(it1 != mapSource.end()) {
		addEdge("*"+name1, name2, Data(0, capacity));
		addEdge(name1, name2, Data(0, capacity));
	    }
	    if(it2 != mapSource.end()) {
		addEdge(name1, "*"+name2, Data(0, capacity));
		addEdge(name1, name2, Data(0, capacity));
	    }
	} else {
	    addEdge(name1, name2, Data(0, capacity));
	}
    }
    //cout << nEdge << endl;
    file.close();
}

void FunctionalFlow::run(int nIteration) {
    //filterNetwork(1);
    //if(flagDebug)
    //	printNodes("Filtered: ");
    //if(flagDebug)
    //	printNetwork("Before iteration");
    for(int i=0; i<nIteration; i++) {
	if(flagDebug)
	    cerr << endl << "----Iteration: " << i << endl;
	updateNetwork();
	//if(flagDebug)
	//    printNetwork("After updateNetwok");
    }	
}

// filter network w.r.t. given parameters 
void FunctionalFlow::filterNetwork(int degreeNodeMin, int degreeNodeMax, float capacityMin, bool flagRemoveSelfEdges) {	
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = pMapIdToNode;
	MapIntToEdge *pMapEdge;
	Vertex *pVertex;
	Edge *pEdge;
	int vDegree; //not uInt !! because of conversion during comparison
	float vReserve, eCapacity;
	vector<int> listNodeToRemove;
	vector< pair<int, int> > listEdgeToRemove;
	for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
		pVertex = itNode->second;		
		pMapEdge = pVertex->pMapIdToEdge;  
		vReserve = pVertex->data.reserve;
		vDegree = pVertex->degree; 		
		if(pVertex->degree != int(pMapEdge->size())) cerr << "Warning: degree inconsistency" << endl;
		if(pVertex->id != itNode->first) cerr << "Warning: vertex id inconsistency" << endl;
		if((vDegree < degreeNodeMin) or (vDegree > degreeNodeMax) ) {
			listNodeToRemove.push_back(itNode->first);
		} else {
			for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
				pEdge = itEdge->second;
				eCapacity = pEdge->data.capacity;
				if(pEdge->idSource != itNode->first or pEdge->idTarget != itEdge->first) cerr << "Warning: edge vertex ids inconsistency" << endl;
				if(flagRemoveSelfEdges and (pEdge->idSource == pEdge->idTarget)) {
					listEdgeToRemove.push_back( pair<int, int>(itNode->first, pEdge->idTarget) );
				} else if(eCapacity < capacityMin) {
					listEdgeToRemove.push_back( pair<int, int>(itNode->first, pEdge->idTarget) );
				}
			}
		}
	}
	for(unsigned int i = 0; i < listNodeToRemove.size(); i++ ) {
        //cout << listNodeToRemove[i] << "\t";
		removeVertex(listNodeToRemove[i]);
	}
	for(unsigned int i = 0; i < listEdgeToRemove.size(); i++ ) {
        //cout << listEdgeToRemove[i].first << "," << listEdgeToRemove[i].second << "\t";
		removeEdge(listEdgeToRemove[i].first, listEdgeToRemove[i].second);
	}
}

void FunctionalFlow::updateNetwork() {
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToEdge::iterator itEdge, itEndEdge;
    MapIntToVertex *pMapNode = pMapIdToNode;
    Vertex *pVertex, *pNeighbor;
    MapIntToEdge *pMapEdge;
    Edge *pEdge;
    //hash_map<int,void*> mapIdProcessed;
    float sumCapacity = 0.0, fluxCurrent = 0.0;
	
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
	//mapIdProcessed[itNode->first] == NULL;
	pVertex = itNode->second;
	if(flagDebug)
	    cout << "Checking node: " << mapIdToName[pVertex->id] << " with reserve " << pVertex->data.reserve << endl;
	pMapEdge = pVertex->pMapIdToEdge;
	sumCapacity = 0.0;
	for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
	    pEdge = itEdge->second;
	    //if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
		pNeighbor = (*pMapNode)[pEdge->idTarget];
		//cout << "- Capacity of edge to neighbor " << mapIdToName[pNeighbor->id] << ": " << pEdge->data.capacity << endl;
		if(pVertex->data.reserve > pNeighbor->data.reserve) {
		    sumCapacity += pEdge->data.capacity;
		}
	    //}
	}
	for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
	    pEdge = itEdge->second;
	    //if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
		pNeighbor = (*pMapNode)[pEdge->idTarget];
		if(flagDebug)
		    cout << "- Calculating score for neighbor: " << mapIdToName[pNeighbor->id] << endl;
		if(pVertex->data.reserve > pNeighbor->data.reserve) {
		    fluxCurrent = 0.0;
		    if(pVertex->data.reserve < sumCapacity) {
			fluxCurrent = pEdge->data.capacity * (pVertex->data.reserve / sumCapacity);
			pNeighbor->data.reserveUpdated += fluxCurrent;
			pVertex->data.reserveUpdated -= fluxCurrent;
		    }
		    else {
			fluxCurrent = pEdge->data.capacity;
			pNeighbor->data.reserveUpdated += fluxCurrent;
			pVertex->data.reserveUpdated -= fluxCurrent;
		    }
		    if(flagDebug)
			cout << "-- Score propagated " << fluxCurrent << ", Neighbor's new reserve " << pNeighbor->data.reserveUpdated << endl;
		}
	    //}
	}
    }
	
    //if(flagDebug)
    //	printNetwork("After update");
    
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
	pVertex = itNode->second;
	if(mapSource.find(itNode->first) != mapSource.end()) {
	    pVertex->data.reserveUpdated = INFINITY;
	}
	pVertex->data.reserve = pVertex->data.reserveUpdated;
    }
}

void FunctionalFlow::outputScores() {
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToVertex *pMapNode = pMapIdToNode;
    Vertex *pVertex;
    std::tr1::unordered_map<int, void*> mapIdProcessed;
    //hash_map<int,void*> mapIdProcessed;
    ofstream fOut;
    fOut.open(fileOutput.c_str());

    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	mapIdProcessed[itNode->first] = NULL;
    	pVertex = itNode->second;
	if(mapIdToName[pVertex->id][0] != '*')
	    fOut << mapIdToName[pVertex->id] << "\t" << (*pMapNode)[pVertex->id]->data.reserveUpdated << endl;	            
    }
    fOut.close();
}

void FunctionalFlow::printNodesAndScores() {
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToVertex *pMapNode = pMapIdToNode;
    Vertex *pVertex;
    std::tr1::unordered_map<int, void*> mapIdProcessed;
    //hash_map<int,void*> mapIdProcessed;

    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	mapIdProcessed[itNode->first] = NULL;
    	pVertex = itNode->second;
	//cout << mapIdToName[pVertex->id] << "\t" << (*pMapNode)[pVertex->id]->data.reserveUpdated << endl;	            
	if(mapIdToName[pVertex->id][0] != '*')
	    cout << mapIdToName[pVertex->id] << "\t" << (*pMapNode)[pVertex->id]->data.reserveUpdated << endl;	            
    }
}

void FunctionalFlow::printNetwork(string explanation) {
    MapIntToVertex::iterator itNode, itEndNode;
    MapIntToEdge::iterator itEdge, itEndEdge;
    MapIntToVertex *pMapNode = pMapIdToNode;
    Vertex *pVertex;
    MapIntToEdge *pMapEdge;
    Edge *pEdge;
    std::tr1::unordered_map<int, void*> mapIdProcessed;
    //hash_map<int,void*> mapIdProcessed;
    
    if(explanation != "") cout << "----" << explanation << ": " << endl;
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	pVertex = itNode->second;
	//cout << mapIdToName[pVertex->id] << endl;
    	pMapEdge = pVertex->pMapIdToEdge;
	for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
	    //cout << mapIdToName[pEdge->idSource] << " " << mapIdToName[pEdge->idTarget] << endl;
            if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
	            cout << mapIdToName[pEdge->idSource] << " " << (*pMapNode)[pEdge->idSource]->data.reserveUpdated 
	 	        	 //<< "(" << (*pMapNode)[pEdge->idSource]->data.reserve << ")\t" 
                                 << "\t"
	 	        	 << mapIdToName[pEdge->idTarget] << " " << (*pMapNode)[pEdge->idTarget]->data.reserveUpdated 
		        	 //<< "(" << (*pMapNode)[pEdge->idTarget]->data.reserve << ")\t" 
                                 << "\t"
                     << "\t" << pEdge->data.capacity << endl;	            
            }
        }
	mapIdProcessed[itNode->first] = NULL;
    }
}

void FunctionalFlow::printNodes(string explanationPrefix) {
    //cout << explanationPrefix << *this << endl;
    MapIntToVertex *pMap = this->pMapIdToNode;
    MapIntToVertex::iterator it, itEnd;
    MapIntToString map2 = this->mapIdToName;
    itEnd = pMap->end(); 
    for(it = pMap->begin(); it != itEnd; ++it)
	cout << "[" << map2[it->first] << "]" << *(it->second) << " - ";
    cout << endl;
}

void FunctionalFlow::printEdges() {
	MapIntToVertex::iterator itNode, itEndNode;
	MapIntToEdge::iterator itEdge, itEndEdge;
	MapIntToVertex *pMapNode = pMapIdToNode;
    Vertex *pVertex;
	MapIntToEdge *pMapEdge;
    Edge *pEdge;
    std::tr1::unordered_map<int, void*> mapIdProcessed;
    //hash_map<int,void*> mapIdProcessed;
    
    for(itNode = pMapNode->begin(), itEndNode = pMapNode->end(); itNode != itEndNode; ++itNode) {
    	mapIdProcessed[itNode->first] = NULL;
    	pVertex = itNode->second;
    	pMapEdge = pVertex->pMapIdToEdge;
	    for(itEdge = pMapEdge->begin(), itEndEdge = pMapEdge->end(); itEdge != itEndEdge; ++itEdge) {
            pEdge = itEdge->second;
            if(mapIdProcessed.find(pEdge->idTarget) == mapIdProcessed.end()) {
	            cout << mapIdToName[pEdge->idSource] << "\t" << mapIdToName[pEdge->idTarget] 
	                 << "\t" << pEdge->data.capacity << endl; 
            }
        }
    }
}

