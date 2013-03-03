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

#ifndef SCORENETWORK_HPP_29045684568920568920456590268
#define SCORENETWORK_HPP_29045684568920568920456590268

#include <cmath> // for INFINITY

#include "graph/Graph.hpp"

#include "Exceptions.hpp"

//#define CONTROL 0
#define sq(x) ((x)*(x))

//typedef boost::unordered_map<Vertex, float> VertexToFloat;
//typedef std::map<Vertex, float> VertexToFloat;
typedef boost::unordered_map<Vertex, float> VertexToFloat;

class ScoreNetwork {
public:
    ScoreNetwork(); 
    ScoreNetwork(std::string fileOutput, bool fUseEdgeScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false);
    ScoreNetwork(std::string fileNode, std::string fileEdge, std::string fileOutput, bool fUseEdgeScore = true, bool fAccumulateToInitialNodeScore = true, bool fResetSeedScoresToInitial = false, bool fVerbose = false);
    virtual ~ScoreNetwork();
    //Graph & getNetwork() { return network; };
    Graph & getNetwork() { return *pNetwork; };
    Graph * getPNetwork() { return pNetwork; };
    void setPNetwork(Graph * pG) { pNetwork = pG; };
    void deletePNetwork() { delete pNetwork; };
    void run(int nRepetition, int nIteration, float tError = PRECISION); 
    void scaleNodeScoresWrapper(ScaleType type) { scaleNodeScores(type); }
    static const float getPrecision() { return PRECISION; }
    /*
    //void outputGraph(Graph const & g, std::string fileName);
    void outputGraph(Graph *g, std::string fileName);
    void printNetwork(std::string explanation = "");
    void printGraph(Graph const & g, std::string explanation = "");
    */
    
protected:
    //int nIteration;
    int repeatCounter;
    int iterationCounter;
    //float minScore;
    //float maxScore;
    bool flagUseEdgeScore;
    bool flagAccumulateToInitialNodeScore;
    bool flagResetSeedScoresToInitial;
    bool flagVerbose;
    std::string outputFile;

    //void filterNetwork(int degreeNodeMin = -INT_INFINITY, int degreeNodeMax = INT_INFINITY, float scoreNodeMin = -INFINITY, float scoreNodeMax = INFINITY, float scoreEdgeMin = -INFINITY, float scoreEdgeMax = INFINITY, bool flagRemoveSelfEdges = true);
    //void scaleScores(float scoreAllowedMaxNode = INFINITY, float scoreAllowedEdgeMax = INFINITY, bool flagUseZScoring = false);
    //void resetScores();

    float calculateErrorAndUpdateScores();
    void updateNetwork();
    void scaleNodeScores(ScaleType typeScale) { getNetwork().scaleVertexScores(typeScale); };
    float transferScore(float score, float a = 0.5, float b = 0.1); 

    float getVertexScoreUpdated(Vertex const v) const { return mapScoreUpdated.at(v); };
    void setVertexScoreUpdated(Vertex const v, float vData) { mapScoreUpdated[v] = vData; };

    virtual void updateNodeScore(Vertex v) {}; // = 0; // not making it abstract (interface) for testing purposes
    virtual void updateEdgeScore(Edge e) {}; 
    virtual void initializeScoring() {}; //= 0;
    virtual void finalizeScoring() {};
    virtual void initializeRepetition() {}; // = 0; 
    virtual void finalizeRepetition() {};
    virtual void initializeIteration() {};
    virtual void finalizeIteration() {}; 


private:
    // CONSTANTS (MACROS)
    //static const float MIN_SCORE = 0.0;
    //static const float MAX_SCORE = 1.0;
    static const float PRECISION = 0.00001f;
    //static const enum TransferType { IDENTITY, POLYNOMIAL, LOGARITHMIC, EXPONENTIAL };
    enum TransferType { IDENTITY, POLYNOMIAL, LOGARITHMIC, EXPONENTIAL } typeTransfer;
    
    // MEMBERS
    //Graph network;
    Graph *pNetwork;
    //unordered_set<string> setSource;
    VertexToFloat mapScoreUpdated;
};

#endif // SCORENETWORK_HPP_

