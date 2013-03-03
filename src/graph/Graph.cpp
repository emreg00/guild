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

#include "Graph.hpp"

#include <ctime>
#include <fstream>

#include "clustering_coefficient.hpp"
//#include <boost/graph/page_rank.hpp>
//#include "page_rank_with_initial_ranks.hpp"
//#include "page_rank_with_priors.hpp"
#include "page_rank_with_weighted_priors.hpp"

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/betweenness_centrality.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>

using namespace boost; 
using namespace std; 

Graph::Graph(int n_node) 
{
    //cout << "Graph" << endl;
    //this->container(n_node); 
    this->container = UndirectedGraph(n_node);
    mapIndex = get(vertex_index, this->container); // get Vertex to Index mapping
    minScore = INFINITY;
    maxScore = -INFINITY;
}

Graph::Graph(Graph const &g) 
{
}

Graph::~Graph() 
{
}

Edge Graph::addEdge(string const &idSource, string const &idTarget, float eData) 
{
    // get u & v
    Vertex u,v;
    Edge e;
    IdUIntIterator it;
    bool found = false;
    //IdVertexIterator it = mapIntexToVertex.left.find(idSource), itEnd = mapIntexToVertex.left.end();
    //if(it != itEnd) u = it->second;
    //else u = addVertex(idSource, )

    //try {
	it = mapIdToIndex.left.find(idSource);
	u = mapIndexToVertex.at(it->second);
	it = mapIdToIndex.left.find(idTarget);
	v =  mapIndexToVertex.at(it->second);
    //} catch (...) {
    //	cerr << "addEdge: missing node" << endl;
    //}

    //tie(e, found) = add_edge(u, v, EdgeProperties(eData, EdgeAttributes(eData)), this->container); // If graph defined with EdgeAttributes instead of EdgeProperties
    //tie(e, found) = add_edge(u, v, eData, this->container); // If graph defined with EdgeProperties - with edge_weight_t
    //tie(e, found) = add_edge(u, v, this->container); // If graph defined with EdgeProperties - edge_index_t instead of edge_weight_t
    tie(e, found) = add_edge(u, v, this->container); // If graph defined with EdgeProperties as EdgeAttributes
    container[e].score = eData;
    //cout << "Added edge: " << getVertexName(v) << ", " << getVertexName(u) << " " << eData << endl;

    //mapEdgeScore mScore = get(&EdgeAttributes::score, this->container);
    //typedef property_map<UndirectedGraph, string VertexAttribtues::id>::type VertexName;
    //VertexName mName = get(VertexAttributes::id, this->container);
    //assert(this->container[u].id == mName(u)); 
    // convert bundled properties to property map (i.e. for algorithms)
    //typedef property_map<UndirectedGraph, double EdgeAttribtues::*>::type EdgeScore;
    //EdgeScore mScore = get(&EdgeAttributes::score, this->container);
    //assert(this->container[e].score == mScore(e)); 
    return e;
}

Vertex Graph::addVertex(string const &id, float vData) 
{
    Vertex v;
    unsigned int i = num_vertices(this->container);
    mapIdToIndex.right.insert(BiStrToUInt::right_value_type(i,id)); //mapIdToIndex.insert(BiStrToUInt::value_type(id, i)); // also using bundled id vertex attribute
    ////v = add_vertex(VertexProperties(i, VertexAttributes(id, vData)), this->container); // If graph defined with VertexAttributes instead of VertexProperties
    v = add_vertex(i, this->container); //v=add_vertex(this->container);
    mapIndexToVertex.insert(make_pair(i, v)); // UIntToVertex::value_type(i, v) 
    ////container[v].id = id;
    //container[v].scoreInitial = vData;
    container[v].score = vData;
    ////cout << get(vertex_index, this->container, v) << endl;
    ////cout << "Added vertex: " << i << " " << id << " " << vData << endl;
    return v;
}

void Graph::loadNodes(string const &fileName) throw(GenericError) 
{
    fstream file;
    file.open(fileName.c_str(), ios::in);
    string name;
    float value = 0.0;
    if(!file) {
	cerr << "Warning: Node file can not be opened " + fileName << endl; 
	throw (new GenericError("File can not be opened: " + fileName)); 
    }
    while(file >> name >> value) {
	addVertex(name, value);
    }
    file.close();
    pair<float,float> tempPair = getMinAndMaxNodeScores();
    minScore = tempPair.first;
    maxScore = tempPair.second;
}

void Graph::loadEdges(string const &fileName, bool flagInverseWeights) throw(GenericError) 
{
    fstream file;
    file.open(fileName.c_str(), ios::in);
    string name1, name2;
    float value = 0.0;
    if(!file) {
	cerr << "Error: Edge file can not be opened " + fileName << endl;
	throw (new GenericError( "File can not be opened: " + fileName)); 
    }
    while(file >> name1 >> value >> name2) {
	if(flagInverseWeights)
	{
	    value = 1/value; // for converting edge relevance to edge weigth e.g. in shortest path
	}
        addEdge(name1, name2, value); 
    }
    file.close();
}

/*
void ScoreNetwork::loadEdgeScores(string const &fileName) throw(GenericError) 
{
    fstream file;
    file.open(fileName.c_str(), ios::in);
    string name1, name2;
    float value = 0.0;
    if(!file) {
	cerr << "Error: file can not be opened " + fileName << endl;
	throw (new GenericError( "File can not be opened: " + fileName)); 
    }
    while(file >> name1 >> name2 >> value) {
        network.setEdgeScore(network.getEdge(network.getVertex(name1), network.getVertex(name2)), value); 
    }
    file.close();
}
*/


map<Vertex, float> Graph::calculateBetweenness() 
{
    map<Vertex, float> vertexToBetweenness;
    associative_property_map< map<Vertex, float> > mapBetweenness(vertexToBetweenness);
    brandes_betweenness_centrality(container, mapBetweenness);
    return vertexToBetweenness;
}

map<Vertex, float> Graph::calculateClusteringCoefficient()
{
    map<Vertex, float> vertexToCoefficient;
    associative_property_map< map<Vertex, float> > mapCoefficient(vertexToCoefficient);
    //clustering_coefficient(container, mapCoefficient, get(vertex_index, container));
    clustering_coefficient(container, mapCoefficient);
    return vertexToCoefficient;
}

void Graph::calculateShortestPath(Vertex v, map<Vertex, float> & vertexToFloat, map<Vertex, Vertex> & vertexToVertex) 
{
    associative_property_map< map<Vertex, Vertex> > mapPredecessor(vertexToVertex);
    associative_property_map< map<Vertex, float> > mapDistance(vertexToFloat);
    dijkstra_shortest_paths(container, v, 
				    predecessor_map(mapPredecessor)
				    .distance_map(mapDistance)
				    .weight_map(get(&EdgeProperties::score, container)));
    return;
}

map<Vertex, float> Graph::calculateShortestPath(Vertex v) 
{
    map<Vertex, Vertex> vertexToVertex;
    associative_property_map< map<Vertex, Vertex> > mapPredecessor(vertexToVertex);
    map<Vertex, float> vertexToFloat;
    associative_property_map< map<Vertex, float> > mapDistance(vertexToFloat);
    //vector< Vertex > vertexToFloat(num_vertices(container));
    //iterator_property_map< vector<Vertex> > mapDistance(vertexToFloat, get(vertex_index, container));
    //make_iterator_map(vertexToFloat.begin(), get(vertex_index, container), vertexToFloat[0]);
    //distance_map(&vertexToFloat[0])

    //vector<double> distances(num_vertices(map));
    //dijkstra_shortest_paths(map, from,
    //	          weight_map(get(&Highway::miles, map))
    // 	          .distance_map(make_iterator_property_map(distances.begin(),
    //						       get(vertex_index, map))));

	
    dijkstra_shortest_paths(container, v, 
				    predecessor_map(mapPredecessor)
				    .distance_map(mapDistance)
				    //.weight_map(get(&EdgeProperties::weight, container)));
				    .weight_map(get(&EdgeProperties::score, container)));
    return vertexToFloat;
}


template <class Visitor>
void Graph::bfsSearch(Vertex v, Visitor bfsVisitor) 
{
    breadth_first_search(container, v, visitor(bfsVisitor));
    return;
}


PredecessorList Graph::getAllShortestPaths(Vertex v) 
{
    PredecessorList vertexToVertices;
    associative_property_map< PredecessorList > mapPredecessorList(vertexToVertices);
    map<Vertex, float> vertexToFloat;
    associative_property_map< map<Vertex, float> > mapDistance(vertexToFloat);

    getAllShortestPaths(v, vertexToFloat, vertexToVertices);

    /*
    map<Vertex, Vertex> vertexToVertex;
    associative_property_map< map<Vertex, Vertex> > mapPredecessor(vertexToVertex);

    dijkstra_shortest_paths(container, v, 
				    predecessor_map(mapPredecessor)
				    .distance_map(mapDistance)
				    .weight_map(get(&EdgeProperties::score, container))
				    .visitor(make_dijkstra_visitor(
						record_all_predecessors(mapPredecessorList, 
								mapDistance, 
								get(&EdgeProperties::score, container),
								on_edge_not_relaxed()))));
    map<Vertex, Vertex>::iterator it, itEnd;
    for(it=vertexToVertex.begin(), itEnd=vertexToVertex.end(); it != itEnd; ++it)
	vertexToVertices[it->first].push_back(it->second);

    */

    return vertexToVertices;
}


void Graph::getAllShortestPaths(Vertex v, PredecessorList & vertexToVertices) 
{
    map<Vertex, float> vertexToFloat;
    associative_property_map< map<Vertex, float> > mapDistance(vertexToFloat);

    getAllShortestPaths(v, vertexToFloat, vertexToVertices);

    return;
}


void Graph::getAllShortestPaths(Vertex v, map<Vertex, float> & vertexToFloat, PredecessorList & vertexToVertices) 
{
    associative_property_map< PredecessorList > mapPredecessorList(vertexToVertices);
    associative_property_map< map<Vertex, float> > mapDistance(vertexToFloat);
    map<Vertex, Vertex> vertexToVertex;
    associative_property_map< map<Vertex, Vertex> > mapPredecessor(vertexToVertex);

    dijkstra_shortest_paths(container, v, 
				    predecessor_map(mapPredecessor)
				    .distance_map(mapDistance)
				    .weight_map(get(&EdgeProperties::score, container))
				    .visitor(make_dijkstra_visitor(
						record_all_predecessors(mapPredecessorList, 
								mapDistance, 
								get(&EdgeProperties::score, container),
								on_edge_not_relaxed()))));

    map<Vertex, Vertex>::iterator it, itEnd;
    for(it=vertexToVertex.begin(), itEnd=vertexToVertex.end(); it != itEnd; ++it)
	vertexToVertices[it->first].push_back(it->second);

    return;
}

/*
// PageRank
void Graph::calculatePageRank(map<Vertex, float> & vertexToFloat, unsigned int nIteration, float dFactor) 
{
    associative_property_map< map<Vertex, float> > mapRank(vertexToFloat);
    //page_rank(container, rank(mapRank));
    page_rank(container, mapRank, graph::n_iterations(nIteration), dFactor);
    return;
}
*/

/*
// PageRank using initial ranks
void Graph::calculatePageRank(map<Vertex, float> & vertexToFloat, unsigned int nIteration, float dFactor) 
{
    map<Vertex, float> vertexToNormalizedScore;
    float sum = 0.0;
    VertexIterator it, itEnd;
    for (tie(it, itEnd) = vertices(this->container); it != itEnd; ++it) {
	sum += getVertexScore(*it);
	//sum += 1; // no initial score consideration
    }
    for (tie(it, itEnd) = vertices(this->container); it != itEnd; ++it) {
	vertexToNormalizedScore[*it] = getVertexScore(*it) / sum;
	//vertexToNormalizedScore[*it] = 1 / sum;
	//cout << getVertexName(*it) << "(" << *it << "): " << vertexToNormalizedScore[*it] << endl;
    }
    associative_property_map< map<Vertex, float> > mapRankInitial(vertexToNormalizedScore);
    associative_property_map< map<Vertex, float> > mapRank(vertexToFloat);
    page_rank(container, mapRankInitial, mapRank, graph::n_iterations(nIteration), dFactor);
    return;
}

// PageRank with priors in addition to using initial ranks (using pre-assignment of initial ranks has very minor effect)
void Graph::calculatePageRankWithPriors(map<Vertex, float> & vertexToFloat, unsigned int nIteration, float dFactor) 
{
    map<Vertex, float> vertexToNormalizedScore;
    float sum = 0.0;
    VertexIterator it, itEnd;
    for (tie(it, itEnd) = vertices(this->container); it != itEnd; ++it) {
	sum += getVertexScore(*it);
    }
    for (tie(it, itEnd) = vertices(this->container); it != itEnd; ++it) {
	vertexToNormalizedScore[*it] = getVertexScore(*it) / sum;
    }
    associative_property_map< map<Vertex, float> > mapRankInitial(vertexToNormalizedScore);
    associative_property_map< map<Vertex, float> > mapRank(vertexToFloat);
    page_rank_with_priors(container, mapRankInitial, mapRank, graph::n_iterations(nIteration), dFactor);
    return;
}
*/

// PageRank with priors in addition to using initial ranks (using pre-assignment of initial ranks has very minor effect)
void Graph::calculatePageRankWithWeightedPriors(map<Vertex, float> & vertexToFloat, unsigned int nIteration, float dFactor) 
{
    map<Vertex, float> vertexToNormalizedScore;
    float sum = 0.0;
    VertexIterator it, itEnd;
    for (tie(it, itEnd) = vertices(this->container); it != itEnd; ++it) {
	sum += getVertexScore(*it);
    }
    for (tie(it, itEnd) = vertices(this->container); it != itEnd; ++it) {
	vertexToNormalizedScore[*it] = getVertexScore(*it) / sum;
    }
    EdgeIterator et, etEnd;
    float max_weight = 0;
    for (tie(et, etEnd) = edges(this->container); et != etEnd; ++et) {
	if(getEdgeScore(*et) > max_weight)
	    max_weight = getEdgeScore(*et);
    }
    for (tie(et, etEnd) = edges(this->container); et != etEnd; ++et) {
	setEdgeScore(*et, getEdgeScore(*et) / max_weight);
    }
    associative_property_map< map<Vertex, float> > mapRankInitial(vertexToNormalizedScore);
    associative_property_map< map<Vertex, float> > mapRank(vertexToFloat);
    page_rank_with_weighted_priors(container, mapRankInitial, get(&EdgeProperties::score, container), mapRank, graph::n_iterations(nIteration), dFactor);
    return;
}

/*
// Stores ranks as scores
void Graph::calculatePageRank() 
{
    map<Vertex, float> vertexToNormalizedScore;
    float sum = 0.0;
    VertexIterator it, itEnd;
    for (tie(it, itEnd) = vertices(this->container); it != itEnd; ++it) {
	sum += getVertexScore(*it);
    }
    for (tie(it, itEnd) = vertices(this->container); it != itEnd; ++it) {
	vertexToNormalizedScore[*it] = getVertexScore(*it) / sum;
    }
    associative_property_map< map<Vertex, float> > mapRankInitial(vertexToNormalizedScore);
    page_rank(container, mapRankInitial, get(&VertexProperties::score, container));
    return;
}
*/

pair<float, float> Graph::getMinAndMaxNodeScores()
{
    VertexIterator it, itEnd;
    float value = 0;
    pair<float, float> result(INFINITY, -INFINITY);
    for(tie(it, itEnd) = getVertexIterator(); it != itEnd; ++it) 
    {
	value = getVertexScore(*it);
	result.second = (value > result.second)?value:result.second;
	result.first = (value < result.first)?value:result.first;
    }
    return result;
}

void Graph::scaleVertexScores(ScaleType typeScale) {
    VertexIterator it, itEnd;
    float value = 0;
    pair<float, float> result;
    switch(typeScale)
    {
	case SCALE_BY_MAX_SCORE:
	    result = getMinAndMaxNodeScores();
	    break;
	case SCALE_BETWEEN_ZERO_AND_ONE:
	    result = getMinAndMaxNodeScores();
	    break;
	case SCALE_BETWEEN_INITIAL_MIN_AND_MAX_SCORE:
	    scaleVertexScores(SCALE_BETWEEN_ZERO_AND_ONE);
	    result = make_pair(minScore, maxScore);
	    break;
	default:
	    cerr << "Unrecognized scaling type" << endl;
	    return;
    }
    for(tie(it, itEnd) = getVertexIterator(); it != itEnd; ++it) 
    {
	value = getVertexScore(*it);
	switch(typeScale)
	{
	    case SCALE_BY_MAX_SCORE:
		value /= result.second;
		break;
	    case SCALE_BETWEEN_ZERO_AND_ONE:
		// If all scores are equal s.t. min==max assign 0.5
		if((result.second - result.first) == 0) 
		    value = 0.5;
		else
		    value = (value - result.first) / (result.second - result.first);
		break;
	    case SCALE_BETWEEN_INITIAL_MIN_AND_MAX_SCORE :
		value = result.first + value * (result.second - result.first);
		break;
	    default:
		cerr << "Unrecognized scaling type" << endl;
		return;
	}
	setVertexScore(*it, value);
    }
    return;
}

void Graph::outputScores(string fileName) const throw(GenericError)
{
    Vertex u;
    VertexIterator vt, vtEnd;
    list< pair<string, float> > scores;
    list< pair<string, float> >::iterator it, itEnd;

    fstream outFile;
    outFile.open(fileName.c_str(), ios::out);
    if(!outFile)
    {
	cerr << "Error in file opening: " << fileName << endl;
	throw (new GenericError("File can not be opened: " + fileName)); 
    }
    for(tie(vt, vtEnd) = getVertexIterator(); vt != vtEnd; ++vt) 
    {
        u = *vt;
        //outFile << getVertexName(u) << "\t" << getVertexScore(u) << endl;
        scores.push_back(make_pair(getVertexName(u), getVertexScore(u)));
    }
    scores.sort(pairSortBySecondPredicate< pair<string, float> >);
    for(it=scores.begin(), itEnd=scores.end(); it != itEnd; ++it)
    {
	outFile << it->first << "\t" << it->second << endl;
    }
    outFile.close();
}

void Graph::printAdjacencyList() const 
{
    //IndexMap mapIndex = get(vertex_index, this->container); // get Vertex to Index mapping
    print_graph(this->container, mapIndex);
}

void Graph::print(bool includeScores) const
{
    Vertex u;
    VertexIterator it, itEnd;
    AdjVertexIterator vt, vtEnd;

    for(tie(it, itEnd) = getVertexIterator(); it != itEnd; ++it) 
    {
        u = *it;
	cout << getVertexName(u);
	if(includeScores) 
	  cout << "(" << getVertexScore(u) << "):";
	for(tie(vt, vtEnd) = getAdjacentVertexIteratorOfVertex(u); vt != vtEnd; ++vt) 
	{
	    cout << "\t" << getVertexName(*vt);
	    if(includeScores)
		cout << "(" << getEdgeScore(u, *vt) << ")";
	}
	cout << endl;
    }
}

void Graph::print() const
{
    Vertex u;
    VertexIterator it, itEnd;
    OutEdgeIterator out, out_end;
    //typedef property_traits<IndexMap>::value_type IndexType;
    //IndexType idx = get(index, this->container, u);
    //IndexMap mapIndex = get(vertex_index, this->container); // get Vertex to Index mapping

    for (tie(it, itEnd) = vertices(this->container); it != itEnd; ++it) {
        u = *it;
        cout << mapIndex[u] << ": " << mapIdToIndex.right.at(mapIndex[u]) << "(" << this->container[u].score << ")\t"; // this->container[u].id <<
        for (tie(out, out_end) = out_edges(u, this->container); out != out_end; ++out) {
	  //cout << mapIndex[source(*out, container)] << "," << mapIndex[target(*out, container)] << "[" << this->container[*out].weight << "]" << "(" << this->container[*out].score << ")\t";
	  cout << mapIndex[source(*out, container)] << "," << mapIndex[target(*out, container)] << "[" << this->container[*out].score << "]" << "(" << this->container[*out].score << ")\t";
        }
	cout << endl;
    }
}

void Graph::printVertices() const
{
    //property_map<UndirectedGraph, vertex_index_t>::type id = get(vertex_index, container);
    cout << "vertices(g) = ";
    graph_traits<UndirectedGraph>::vertex_iterator vi;
    for (vi = vertices(container).first; vi != vertices(container).second; ++vi)
	cout << mapIndex[*vi] <<  " ";
    cout << endl;
}

void test_Graph() {
    //Graph g();
    //Graph *g2;
    //g2 = new Graph();
    //Graph g(2);
    Graph g;
    g.addVertex("u", 0.2);
    g.addVertex("v", 0.5);
    g.addVertex("p", 0.7);
    g.addVertex("q", 0.8);
    g.addEdge("u", "v", 1);
    g.addEdge("u", "q", 0.4);
    g.addEdge("v", "p", 0.6);
    g.printAdjacencyList();
    g.print();
}

void test_boost_create_and_iterate()
{
    Vertex u,v;
    VertexIterator it, itEnd;
    OutEdgeIterator out, out_end;

    const unsigned int V = 10;

    UndirectedGraph undigraph(0);
    bool firstTime = true;

    clock_t t1 = clock();

    //float score = 0.2, weight = 0.5;
    string s = "a";
    
    property_map<UndirectedGraph, vertex_index_t>::type v_index = get(vertex_index, undigraph); 

    for(unsigned int i=0; i<V; i++) 
    {
        //u = add_vertex(VertexProperties(num_vertices(undigraph), VertexAttributes(s, score*i)), undigraph);
        //u=add_vertex(undigraph);
        //u = add_vertex(property<vertex_index_t, int>(num_vertices(undigraph)), undigraph);
        //u = add_vertex(VertexProperties(i, s, score), undigraph);
        //u = add_vertex(VertexAttributes(i, s, score), undigraph);
        u = add_vertex(i, undigraph);
        //u = add_vertex(i, undigraph);
        firstTime=true;
        for(unsigned int j=i; j<V; j++) 
        {
	  if(firstTime) {
	      v=add_vertex(undigraph);
	      //undigraph[v].id = "a";
	      undigraph[v].score = 0.2*j;
	      //v = add_vertex(VertexProperties(num_vertices(undigraph), VertexAttributes(s, score)), undigraph);
	      //v = add_vertex(VertexProperties(j, s, score), undigraph);
	      //v = add_vertex(property<vertex_index_t, int>(j), undigraph);
	      //v = add_vertex(j, undigraph);
	      firstTime=false;
	  }
	  //add_edge(u, v, EdgeProperties(weight), undigraph); 
	  //add_edge(u, v, EdgeProperties(weight, EdgeAttributes(weight*(i+j))), undigraph); 
	  //add_edge(u, v, weight, undigraph); 
	  add_edge(u, v, undigraph); 
	  //add_edge(u, v, undigraph); 
        }
    }

    int c = 0;

    clock_t t2 = clock();
    cout << (t2-t1) << " (" << (t2-t1)/(double)CLOCKS_PER_SEC << "s)" << endl;


    for (tie(it, itEnd) = vertices(undigraph); it != itEnd; ++it, ++c) {
        u = *it;
        //cout << static_cast<unsigned int>(v_index[u]) << endl;
        //cout << get(vertex_index, undigraph, u) << endl;
	//
        //cout << endl << undigraph[u].id;
        for (tie(out, out_end) = out_edges(u, undigraph); out != out_end; ++out) {
	  //cout << "\t" << undigraph[*out].weight << "\t";
	  cout << "\t" << undigraph[*out].score << "\t";
        }
    }

    clock_t t3 = clock();
    cout << (t3-t2) << " (" << (t3-t2)/(double)CLOCKS_PER_SEC << "s)" << endl;
}

template <class Pair>
bool pairSortBySecondPredicate(const Pair& lhs, const Pair& rhs)
{
      return lhs.second < rhs.second;
}

