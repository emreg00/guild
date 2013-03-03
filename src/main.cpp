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

//#include "ScoreNetwork.hpp"
#include "Netshort.hpp"
#include "Netscore.hpp"
#include "Netrank.hpp"
#include "Netzcore.hpp"
#include "Netrandom.hpp"
//#include "Netzscore.hpp"
#include "Netween.hpp"
#include "Netlink.hpp"

#include "functionalFlow/FunctionalFlow.h"

#include <iostream>
#include <sstream>
#include <ctime>
#include <unistd.h>

using namespace std;

void runScoreNetwork();
void runNetshort(string, string, string);
void runNetscore(string, string, string, unsigned int, unsigned int);
void runNetrank(string, string, string, unsigned int);
void runNetzcore(string, string, string, string, unsigned int, unsigned int, unsigned int);
void runNetrandom(string, string, string );
void runNetzscore(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled, unsigned int nRepetition, unsigned int nIteration); 
void runNetz1score(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled, unsigned int nRepetition, unsigned int nIteration); 
void runNetween(string fileNode, string fileEdge, string fileOutput, float threshold);
void runNetlink(string fileNode, string fileEdge, string fileOutput, float threshold);
void functionalFlow(string fileProtein, string fileInteraction, string fileOutput, int nIteration, float seedScoreThreshold);

void printHelp(char *name) {
    cout << endl << name  
	 << " [ Copyleft (GPLv3) - 2011 - Emre Guney (Universitat Pompeu Fabra) ]\n" << endl
	 << " Arguments: " << endl
	 //<< "-s <prioritization_method>{NetScore:s|NetZcore:z|NetShort:d|fFlow:f|NetRank:r|NetRandom:x|NetZScore:h|NetZ1Score:1|NetWeen:w|NetLink:l}\n"
	 << "\t -s <prioritization_method>{NetScore:s|NetZcore:z|NetShort:d|fFlow:f|NetRank:r}\n"
	 << "\t -n <node_file>\n" 
	 << "\t -e <edge_file>\n" 
	 << "\t -o <output_file>\n"
	 << "\t -i <number_of_iterations>\n"
	 << "\t -r <number_of_repetitions>\n"
	 << "\t -t <seed_score_threshold>\n"
	 << "\t -x <number_of_sampled_graphs>\n" 
	 << "\t -d <sampled_graph_prefix>\n"
	 << "\t -h\n"
	 << endl;
}

int main(int argc, char **argv)
{
    static const char *optString = "n:e:i:s:o:r:x:d:t:h?";
    int opt = 0;
    string fileNode, fileEdge, fileOutput, prefixSampling;
    enum ScoringMethod { NETSCORE = 's', NETSHORT = 'd', NETZCORE = 'z', FFLOW = 'f', NETRANK = 'r', NETRANDOM = 'x', NETZSCORE = 'h', NETZ1SCORE = '1', NETWEEN = 'w', NETLINK = 'l' };
    unsigned int nIteration = 1, nRepetition = 1, nSampled = 0;  
    char scoring = 's';
    float threshold = 0.01;

    if(argc < 9) {
	printHelp(argv[0]);
	return 1;
    }

    opt = getopt( argc, argv, optString );
    while( opt != -1 ) {
	istringstream iss;
        switch( opt ) {
            case 'n':
            	fileNode = string(optarg); 
                break;
            case 'e':
            	fileEdge = string(optarg);
                break;
            case 'i':
                iss.str(optarg);
            	iss >> nIteration;
                break;
            case 'r':
                iss.str(optarg);
            	iss >> nRepetition;
                break;
            case 's':
                iss.str(optarg);
                iss >> scoring;
                break;
            case 'o':
            	fileOutput = string(optarg);
                break;
	    case 'x':
                iss.str(optarg);
            	iss >> nSampled;
                break;
            case 'd':
            	prefixSampling = string(optarg); 
                break;
            case 't':
                iss.str(optarg);
            	iss >> threshold;
                break;
	    case 'h':
	    case '?':
		printHelp(argv[0]);
                return 0;
            default:
                break;
        }
        opt = getopt( argc, argv, optString );
    }
        
    cout << "Arguments: scoring type " << scoring << ", nRepetition " << nRepetition << ", nIteration " << nIteration << ", nodeFile " << fileNode << ", edgeFile " << fileEdge << ", outputFile " << fileOutput << endl;
    clock_t t1 = clock();
    if(scoring == NETSHORT)
	runNetshort(fileNode, fileEdge, fileOutput);
    else if (scoring == NETSCORE)
	runNetscore(fileNode, fileEdge, fileOutput, nRepetition, nIteration);
    else if (scoring == NETRANK)
	runNetrank(fileNode, fileEdge, fileOutput, nIteration);
    else if (scoring == NETZCORE)
	runNetzcore(fileNode, fileEdge, fileOutput, prefixSampling, nSampled, nRepetition, nIteration); 
    else if (scoring == NETRANDOM)
	runNetrandom(fileNode, fileEdge, fileOutput); 
    else if (scoring == NETZSCORE)
	runNetzscore(fileNode, fileEdge, fileOutput, prefixSampling, nSampled, nRepetition, nIteration); 
    else if (scoring == NETZ1SCORE)
	runNetz1score(fileNode, fileEdge, fileOutput, prefixSampling, nSampled, nRepetition, nIteration); 
    else if (scoring == NETWEEN)
	runNetween(fileNode, fileEdge, fileOutput, threshold); 
    else if (scoring == NETLINK)
	runNetlink(fileNode, fileEdge, fileOutput, threshold); 
    else if (scoring == FFLOW)
	functionalFlow(fileNode, fileEdge, fileOutput, nIteration, threshold);
    else
	//runScoreNetwork();
	cerr << "Unrecognized scoring type!" << endl;
    clock_t t2 = clock();
    cout << "Time: " << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;
    return 0;
}

void runNetlink(string fileNode, string fileEdge, string fileOutput, float threshold) 
{
    // flagAccumulateToInitialNodeScore, flagVerbose
    Netlink sN(fileNode, fileEdge, fileOutput, threshold, true, false);
    sN.run(1, 1);
}

void runNetween(string fileNode, string fileEdge, string fileOutput, float threshold) 
{
    // seedScoreThreshold, flagAccumulateToInitialNodeScore, flagVerbose
    Netween sN(fileNode, fileEdge, fileOutput, threshold, true, false);
    sN.run();
}

void runNetz1score(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled, unsigned int nRepetition, unsigned int nIteration) 
{
    // flagUseEdgeScore, flagAccumulateToInitialNodeScore, flagResetSeedScoresToInitial, flagVerbose
    Netzcore sNz(fileNode, fileEdge, fileOutput, prefixSampledGraphs, nSampled, true, true, false, false);
    Netscore sNs(fileOutput, true, true, false, false);
    sNs.setPNetwork(sNz.getPNetwork());

    // Normalize score in the begining once
    sNz.run(1, 1);
    sNs.run(nRepetition, nIteration);

    sNz.deletePNetwork();
}

void runNetzscore(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled, unsigned int nRepetition, unsigned int nIteration) 
{
    // flagUseEdgeScore, flagAccumulateToInitialNodeScore, flagResetSeedScoresToInitial, flagVerbose
    //Netzscore sN(fileNode, fileEdge, fileOutput, prefixSampledGraphs, nSampled, true, false, false, true);
    Netzcore sNz(fileNode, fileEdge, fileOutput, prefixSampledGraphs, nSampled, true, false, false, false);
    Netscore sNs(fileOutput, true, true, false, false);
    sNs.setPNetwork(sNz.getPNetwork());

    // Normalize score at each repetition
    for(unsigned int r=0; r < nRepetition; ++r)
    {
	sNz.run(1, 1);
	sNs.run(1, nIteration);
    }
    sNz.deletePNetwork();
}

void runNetzcore(string fileNode, string fileEdge, string fileOutput, string prefixSampledGraphs, unsigned int nSampled, unsigned int nRepetition, unsigned int nIteration) 
{
    // flagUseEdgeScore, flagAccumulateToInitialNodeScore, flagResetSeedScoresToInitial, flagVerbose
    Netzcore sN(fileNode, fileEdge, fileOutput, prefixSampledGraphs, nSampled, true, true, false, false);
    sN.run(nRepetition, nIteration);
}

void runNetscore(string fileNode, string fileEdge, string fileOutput, unsigned int nRepetition, unsigned int nIteration) 
{
    // flagUseEdgeScore, flagAccumulateToInitialNodeScore, flagResetSeedScoresToInitial, flagVerbose
    Netscore sN(fileNode, fileEdge, fileOutput, true, true, false, false);
    sN.run(nRepetition, nIteration);
}

void runNetshort(string fileNode, string fileEdge, string fileOutput) 
{
    //clock_t t1 = clock();
    // flagAccumulateToInitialNodeScore
    Netshort sN(fileNode, fileEdge, fileOutput, true);
    //clock_t t2 = clock();
    //cout << "Time to load graph: " << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;
    sN.run();
    //cout << sN.testRun() << endl;
}

void runNetrank(string fileNode, string fileEdge, string fileOutput, unsigned int nIteration) 
{
    // flagAccumulateToInitialNodeScore
    Netrank sN(fileNode, fileEdge, fileOutput, true);
    sN.run(20, 0.85); // number of pagerank iterations, damping factor
}

void runNetrandom(string fileNode, string fileEdge, string fileOutput) 
{
    Netrandom sN(fileNode, fileEdge, fileOutput, true);
    sN.run();
}

void functionalFlow(string fileProtein, string fileInteraction, string fileOutput, int nIteration, float seedScoreThreshold) {
	// flagRemoveSeedEffect, flagDebug
	FunctionalFlow fFlow = FunctionalFlow(fileProtein, fileInteraction, fileOutput, seedScoreThreshold, false, false);
	fFlow.run(nIteration);
	fFlow.outputScores();
}

/*
void runScoreNetwork() 
{
    //ScoreNetwork sN("../../data/toy_data/test_proteins_small.txt", "../../data/toy_data/test_interactions_small.txt", "../../data/out_test.txt");
    ScoreNetwork sN("../../data/toy_data/test_proteins_small.txt", "../../data/toy_data/test_interactions_small.txt", "../../data/out_test.txt", true, true, false, false);
    if(sN.getNetwork().getSize() ==  5)
	cout << "Correcto" << endl;
    else {
	cout << sN.getNetwork().getSize() << endl;
    }

    map<Vertex, float>::iterator it, itEnd;
    map<Vertex, float> mapB = sN.getNetwork().calculateBetweenness();

    for(it=mapB.begin(), itEnd=mapB.end(); it != itEnd; ++it)
	cout << sN.getNetwork().getVertexIndex(it->first) << ": " << it->second << endl;

    map<Vertex, float> mapC = sN.getNetwork().calculateClusteringCoefficient();

    for(it=mapC.begin(), itEnd=mapC.end(); it != itEnd; ++it)
	cout << sN.getNetwork().getVertexIndex(it->first) << ": " << it->second << endl;

    map<Vertex, float> mapS = sN.getNetwork().calculateShortestPath(sN.getNetwork().getVertex("v1"));

    for(it=mapS.begin(), itEnd=mapS.end(); it != itEnd; ++it)
	cout << sN.getNetwork().getVertexIndex(it->first) << ": " << it->second << endl;
}
*/

