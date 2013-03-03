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

#ifndef GRAPH_HPP_13458158201345810349850213485
#define GRAPH_HPP_13458158201345810349850213485


#include <boost/bimap.hpp>
#include <boost/config.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>

#include "../Exceptions.hpp"

#include <map>
#include <vector>

#include <cmath> // for INFINITY

#include <boost/graph/visitors.hpp>
#include <boost/graph/graph_traits.hpp>

typedef enum ScaleType { SCALE_BY_MAX_SCORE, SCALE_BETWEEN_ZERO_AND_ONE, SCALE_BETWEEN_INITIAL_MIN_AND_MAX_SCORE } ScaleType;

template <class Pair>
bool pairSortBySecondPredicate(const Pair& lhs, const Pair& rhs);

// Define features as built-in types
struct vertex_score_t { 
    typedef boost::vertex_property_tag kind;
};

struct edge_score_t { 
    typedef boost::edge_property_tag kind;
};
//typedef property< vertex_score_t, float> VertexAttributes;

// Moved to Netscore
//typedef std::pair<float, int> Message;
//typedef boost::tuple<float, int, int> Message;
//typedef boost::unordered_map<unsigned int,  Message> UIntToMessage;

// Define features as bundled types - somewhat suggested
struct VertexAttributes {
    ////unsigned int index; // using built-in index attribute
    ////std::string id; // not necessary to store it here since there is the mapIdToIndex in graph
    //float scoreInitial; // Should the error need not be calculated, this would be unnecessary - one would use score as initial scores and scoreUpdated for updated scores
    float score;
    //float scoreUpdated;
    ////UIntToMessage mapMessage;
    //UIntToMessage* pMapMessage; // naming consideration: p_mapMessage;
    //VertexAttributes(unsigned int i=0, std::string id="", float s=0):index(i), id(id), score(s) {};
    VertexAttributes(float s=0):score(s) {};
    ////VertexAttributes(std::string id="", float s=0):id(id), score(s) {};
    ////std::ostream & operator<<(std::ostream & output) { output << id << " (" << score << ")\n"; return output; };
};

struct EdgeAttributes {
    //float weight;
    float score;
    //EdgeAttributes(float s=0, float w=1):score(s),weight(w) {};
    EdgeAttributes(float s=0):score(s) {};
    //std::ostream & operator<<(std::ostream & output) { output << "edge " << " (" << weight << ")\n"; return output; };
};

// Define graph
typedef boost::property<boost::vertex_index_t, unsigned int, VertexAttributes > VertexProperties;
//typedef boost::property<boost::edge_index_t, unsigned int, EdgeAttributes> EdgeProperties;
//typedef boost::property<boost::edge_weight_t, float, EdgeAttributes> EdgeProperties;
typedef EdgeAttributes EdgeProperties; // to make this work: typedef boost::property<EdgeAttributes> EdgeProperties;
//typedef adjacency_list < vecS, hash_setS, undirectedS, VertexAttributes, EdgeAttributes> UndirectedGraph; //! strange but works
typedef boost::adjacency_list < boost::vecS, boost::hash_setS, boost::undirectedS, VertexProperties, EdgeProperties> UndirectedGraph;

typedef boost::graph_traits < UndirectedGraph >::vertex_descriptor Vertex;
typedef boost::graph_traits < UndirectedGraph >::edge_descriptor Edge;
typedef boost::graph_traits< UndirectedGraph >::vertex_iterator VertexIterator;
typedef boost::graph_traits< UndirectedGraph >::adjacency_iterator AdjVertexIterator;
typedef boost::graph_traits < UndirectedGraph >::out_edge_iterator OutEdgeIterator;
typedef boost::graph_traits < UndirectedGraph >::edge_iterator EdgeIterator;

typedef boost::property_map<UndirectedGraph, boost::vertex_index_t>::type IndexMap;

typedef boost::bimap<std::string, unsigned int> BiStrToUInt;
typedef boost::unordered_map<unsigned int, Vertex> UIntToVertex;
typedef BiStrToUInt::left_value_type StrToUInt;
typedef BiStrToUInt::right_value_type UIntToStr;
typedef BiStrToUInt::left_const_iterator IdUIntIterator;
typedef BiStrToUInt::right_const_iterator UIntIdIterator;

typedef std::map<Vertex, std::vector<Vertex> > PredecessorList;

class Graph {
    //friend ostream& operator<<(ostream& output, const Graph& g);
public:
    // MACROS
    //const int MAX_V = 10;

    // Consturctor & Copy cons & Desctructor
    Graph(int n_node=0);
    Graph(Graph const &g);
    virtual ~Graph();
    // Add node & edge
    Vertex addVertex(std::string const &id, float vData);
    Edge addEdge(std::string const &idSource, std::string const &idTarget, float eData);
    void loadNodes(std::string const &fileName) throw(GenericError);
    void loadEdges(std::string const &fileName, bool flagInverseWeights = false) throw(GenericError);    
    //void loadEdgeScores(std::string const &fileName) throw(GenericError);    
    // Getters (general)
    unsigned int getSize() const { return boost::num_vertices(container); }
    unsigned int getNumberOfEdges() const { return boost::num_edges(container); }
    // Vertex Getters & Setters // suggestion: inlines better in cpp 
    Vertex getVertex(std::string const &id) const { return (mapIndexToVertex.find(mapIdToIndex.left.at(id)))->second; }; // was using instead of at: (.left.find(id))->second
    Vertex getVertex(unsigned int index) const { return mapIndexToVertex.find(index)->second; }; 
    //float getVertexScoreInitial(Vertex const v) const { return container[v].scoreInitial; };
    //void setVertexScoreInitial(Vertex const v, float vData) { container[v].scoreInitial = vData; };
    float getVertexScore(Vertex const v) const { return container[v].score; };
    void setVertexScore(Vertex const v, float vData) { container[v].score = vData; };
    //float getVertexScoreUpdated(Vertex const v) const { return container[v].scoreUpdated; };
    //void setVertexScoreUpdated(Vertex const v, float vData) { container[v].scoreUpdated = vData; };
    ////UIntToMessage & getVertexMessageMap(Vertex const v) { return container[v].mapMessage; };
    //UIntToMessage * getVertexMessageMap(Vertex const v) { return container[v].pMapMessage; };
    //void createVertexMessageMap(Vertex const v) { container[v].pMapMessage = new UIntToMessage(); };
    ////float getVertexScore(std::string const &id) const { getVertexScore(getVertex(id)); }; // does not work
    float getVertexScore(std::string const &id) const { return container[(mapIndexToVertex.find(mapIdToIndex.left.at(id)))->second].score; };
    unsigned int getVertexIndex(std::string const &id) const { return mapIdToIndex.left.at(id); };
    unsigned int getVertexIndex(Vertex const v) const { return mapIndex[v]; };
    std::string getVertexName(Vertex const v) const { return mapIdToIndex.right.at(mapIndex[v]); };
    std::string getVertexName(unsigned int index) const { return mapIdToIndex.right.at(index); };
    // Edge Getters & Setters
    Edge getEdge(Vertex vSource, Vertex vTarget) { return edge(vSource, vTarget, container).first; };
    //Edge getEdge(std::string const & idSource, std::string & idTarget) { edge(getVertex(idSource), getVertex(idTarget), container); };
    float getEdgeScore(Vertex vSource, Vertex vTarget) const { return container[edge(vSource, vTarget, container).first].score; };
    float getEdgeScore(Edge e) const { return container[e].score; };
    void setEdgeScore(Edge e, float score) { container[e].score = score; };
    //float getEdgeWeight(Edge e) const { return container[e].weight; };
    Vertex getSourceVertex(Edge e) const { return source(e, container); };
    Vertex getTargetVertex(Edge e) const { return target(e, container); };
    // Iteration
    std::pair<VertexIterator, VertexIterator> getVertexIterator() const { return vertices(this->container); };     
    std::pair<AdjVertexIterator, AdjVertexIterator> getAdjacentVertexIteratorOfVertex(Vertex const v) const { return adjacent_vertices(v, this->container); };
    std::pair<EdgeIterator, EdgeIterator> getEdgeIterator() const { return edges(this->container); };
    std::pair<OutEdgeIterator, OutEdgeIterator> getEdgeIteratorOfVertex(Vertex const v) const { return out_edges(v, this->container); };
    // Various graph algorithms
    std::map<Vertex, float> calculateBetweenness();
    std::map<Vertex, float> calculateClusteringCoefficient();
    std::map<Vertex, float> calculateShortestPath(Vertex v);
    void calculateShortestPath(Vertex v, std::map<Vertex, float> & vertexToFloat, std::map<Vertex, Vertex> & vertexToVertex);
    template <class Visitor> void bfsSearch(Vertex v, Visitor bfsVisitor); 
    PredecessorList getAllShortestPaths(Vertex v);
    void getAllShortestPaths(Vertex v, PredecessorList & vertexToVertices);
    void getAllShortestPaths(Vertex v, std::map<Vertex, float> & vertexToFloat, PredecessorList & vertexToVertices);
    //void calculatePageRank(std::map<Vertex, float> & vertexToFloat, unsigned int nIteration = 20, float dFactor = 0.85);
    //void calculatePageRankWithPriors(std::map<Vertex, float> & vertexToFloat, unsigned int nIteration = 20, float dFactor = 0.85);
    void calculatePageRankWithWeightedPriors(std::map<Vertex, float> & vertexToFloat, unsigned int nIteration = 20, float dFactor = 0.85);
    // Processing graph scores
    void scaleVertexScores(ScaleType typeScale);
    std::pair<float, float> getMinAndMaxNodeScores();
    // Graph outputting & printing 
    void outputScores(std::string outputFile) const throw(GenericError);
    void printAdjacencyList() const;
    void print(bool includeScores) const;
    void print() const;
    void printVertices() const;

private:
    UndirectedGraph container; // Boost undirected graph
    //property_map<UndirectedGraph, vertex_score_t>::type mapScore = get(vertex_score_t(), this->container); 
    BiStrToUInt mapIdToIndex; // Vertex name to index bi-mapping
    UIntToVertex mapIndexToVertex; // Vertex index to descriptor mapping 
    IndexMap mapIndex; // built-in Vertex descriptor to index property map

    float minScore;
    float maxScore;
};


template <class PredecessorListMap, class DistanceMap, class WeightMap, class Tag>
struct all_predecessors_recorder 
    : public boost::base_visitor<all_predecessors_recorder < PredecessorListMap, DistanceMap, WeightMap, Tag> >
{
    typedef Tag event_filter;
    
    all_predecessors_recorder(PredecessorListMap pa, DistanceMap da, WeightMap wa) : m_predecessor(pa), m_distance(da), m_weight(wa) { }
    template <class Edge, class Graph>
	void operator()(Edge e, const Graph& g) {
	    typename boost::graph_traits<Graph>::vertex_descriptor
			    u = boost::source(e, g), v = boost::target(e, g);
	    //put(m_distance, v, get(m_distance, u) + 1);
	    if(boost::get(m_distance, u) + get(m_weight, e) <= boost::get(m_distance, v))
		//put(m_predecessor, target(e, g), source(e, g));
		(m_predecessor[target(e,g)]).push_back(source(e,g));
	}
private:
    PredecessorListMap m_predecessor;
    DistanceMap m_distance;
    WeightMap m_weight;
};

template <class PredecessorListMap, class DistanceMap, class WeightMap, class Tag>
all_predecessors_recorder<PredecessorListMap, DistanceMap, WeightMap, Tag>
record_all_predecessors(PredecessorListMap pa, DistanceMap da, WeightMap wa, Tag) {
    return all_predecessors_recorder<PredecessorListMap, DistanceMap, WeightMap, Tag> (pa, da, wa);
}


/*
inline
float Graph::getVertexScore(std::string id) 
{
    Vertex v;
    v = (mapIndexToVertex.find(mapIdToIndex.left.at(id)))->second;
    return container[v].score;
}

inline
float Graph::getVertexScore(Vertex v) 
{
    return container[v].score;
}

inline
unsigned int Graph::getVertexIndex(string id) {
    return mapIdToIndex.left.at(id);
}

inline
std::pair<VertexIterator, VertexIterator> Graph::getVertexIterator() {
    return vertices(this->container);
}

inline
std::pair<EdgeIterator, EdgeIterator> Graph::getEdgeIteratorOfVertex(Vertex v) {
    return out_edges(v, this->container);
}
*/


void test_Graph();
void test_boost_create_and_iterate();

#endif //GRAPH_HPP_

