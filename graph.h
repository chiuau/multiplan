#ifndef SIM_GRAPH_H
#define SIM_GRAPH_H

#include <iostream>
#include <stdexcept>
#include <list>
#include <vector>
#include <unordered_map>
#include <queue>
#include <limits>
#include <cassert>


template<typename V, typename E>
class AdjListGraph {

  struct VertexObj;
  using VertexObjList = std::list<VertexObj>;
  using VertexObjItr = typename VertexObjList::iterator;

  struct EdgeObj;
  using EdgeObjList = std::list<EdgeObj>;
  using EdgeObjItr = typename EdgeObjList::iterator;

public:

  class Vertex;
  using VertexList = std::list<Vertex>;
  using VertexItr = typename VertexList::iterator;

  class Edge;
  using EdgeList = std::list<Edge>;
  using EdgeItr = typename EdgeList::iterator;

  using EVList = std::list<std::pair<Edge,Vertex>>;
  using EVItr = typename EVList::iterator;

  static const Vertex NULL_VERTEX;

private:

  /*
   * VertexObj stores data of an vertex.
   */
  struct VertexObj {
    V elt;                                   // the element stored at this vertex
    int vertex_id;                           // the vertex id.
    VertexObjItr vertex_obj_collection_itr;  // the position of this vertex object in vertex_obj_collection
    VertexItr vertex_collection_itr;         // the position of this vertex in vertex_collection
    EVList outgoing_ev_list;                 // the list of outgoing [ edge , vertex ]
    EVList incoming_ev_list;                 // the list of incoming [ edge , vertex ]

    explicit VertexObj(V elt) : elt{elt} {}
  };

  /*
   * EdgeObj stores the data of an edge.
   */
  struct EdgeObj {
    E elt;                                // the element stored at this edge
    Vertex orig_vertex;                   // the vertex at the origin
    Vertex dest_vertex;                   // the vertex at the destination
    EdgeObjItr edge_obj_collection_itr;   // the position of this edge object in edge_obj_collection
    EdgeItr edge_collection_itr;          // the position of this Edge in edge_collection
    EVItr orig_outgoing_ev_itr;           // the iterator to orig_vertex.outgoing_ev_list
    EVItr dest_incoming_ev_itr;           // the iterator to dest_vertex.incoming_ev_list

    explicit EdgeObj(E elt) : elt{elt} {}
  };

public:

  /*
   * Vertex - a wrapper of a reference to an VertexObj
   */
  class Vertex {

    VertexObj* vertex_obj_ptr;

  public:

    Vertex() : vertex_obj_ptr{nullptr} {}

    explicit Vertex(VertexObj& vertex_obj) : vertex_obj_ptr{&vertex_obj} {}
    
    /*
     * Return the ID of this vertex.
     */
    int id() const { return vertex_obj_ptr->vertex_id; }

    /*
     * Whether the vertex is null.
     */
    bool isNull() const { return vertex_obj_ptr == nullptr; }

    /*
     * Return the element stored at this vertex.
     */
    V& get() const { return vertex_obj_ptr->elt; }

    /*
     * Return a list of outgoing [ edge , vertex ].
     */
    const EVList& outgoingEVList() const { return vertex_obj_ptr->outgoing_ev_list; }

    /*
     * Return a list of incoming [ edge , vertex ].
     */
    const EVList& incomingEVList() const { return vertex_obj_ptr->incoming_ev_list; }

    /*
     * Check whether this vertex is the same as the given vertex
     *
     * v - the given vertex
     * Return true if this vertex is the same as the given vertex
     */
    bool operator==(const Vertex& other) const {
      return vertex_obj_ptr == other.vertex_obj_ptr;
    }

    friend AdjListGraph<V,E>;
  };


  /*
   * Edge - a wrapper of a reference to an EdgeObj
   */
  class Edge {

    EdgeObj* edge_obj_ptr;

  public:

    Edge() : edge_obj_ptr{nullptr} {}

    explicit Edge(EdgeObj& edge_obj) : edge_obj_ptr{&edge_obj} {}

    /*
     * Whether the vertex is null.
     */
    bool isNull() const { return edge_obj_ptr == nullptr; }

    /*
     * Return the element stored at this edge.
     */
    E& get() const { return edge_obj_ptr->elt; }

    /*
     * Return the vertex at the origin of this edge.
     */
    Vertex origVertex() const { return edge_obj_ptr->orig_vertex; }

    /*
     * Return the vertex at the destination of this edge.
     */
    Vertex destVertex() const { return edge_obj_ptr->dest_vertex; }

    /*
     * Check whether this edge is the same as the given edge.
     *
     * edge - the given edge
     * Return true if this edge is the same as the given edge.
     */
    bool operator==(const Edge& other) const {
      return edge_obj_ptr == other.edge_obj_ptr;
    }

    friend AdjListGraph<V,E>;
  };

public:

  AdjListGraph() = default;

  AdjListGraph(const AdjListGraph& other) :
      next_vertex_id{other.next_vertex_id}
  {
    for(auto& v_obj_other : other.vertex_obj_collection) {
      vertex_obj_collection.emplace_back(v_obj_other.elt);
      auto& v_obj = vertex_obj_collection.back();
      vertex_collection.emplace_back(v_obj);
      vid_to_vertex[v_obj_other.vertex_id] = vertex_collection.back();
      v_obj.vertex_id = v_obj_other.vertex_id;
      v_obj.vertex_obj_collection_itr = vertex_obj_collection.end();
      v_obj.vertex_obj_collection_itr--;
      v_obj.vertex_collection_itr = vertex_collection.end();
      v_obj.vertex_collection_itr--;
      // skip v_obj.outgoing_ev_list and v_obj.incoming_ev_list
    }

    for(auto& e_obj_other : other.edge_obj_collection) {
      edge_obj_collection.emplace_back(e_obj_other.elt);
      auto& e_obj = edge_obj_collection.back();
      edge_collection.emplace_back(e_obj);
      auto& e = edge_collection.back();
      e_obj.orig_vertex = vid_to_vertex[e_obj_other.orig_vertex.id()];
      e_obj.dest_vertex = vid_to_vertex[e_obj_other.dest_vertex.id()];
      e_obj.edge_obj_collection_itr = edge_obj_collection.end();
      e_obj.edge_obj_collection_itr--;
      e_obj.edge_collection_itr = edge_collection.end();
      e_obj.edge_collection_itr--;
      // establish the edges_for_decoration
      e_obj.orig_vertex.vertex_obj_ptr->outgoing_ev_list.emplace_back(e, e_obj.dest_vertex);  // ***
      e_obj.orig_outgoing_ev_itr = e_obj.orig_vertex.vertex_obj_ptr->outgoing_ev_list.end();
      e_obj.orig_outgoing_ev_itr--;
      e_obj.dest_vertex.vertex_obj_ptr->incoming_ev_list.emplace_back(e, e_obj.orig_vertex); // ***
      e_obj.dest_incoming_ev_itr = e_obj.dest_vertex.vertex_obj_ptr->incoming_ev_list.end();
      e_obj.dest_incoming_ev_itr--;
    }
    // the orderings of v_obj.outgoing_ev_list and v_obj.incoming_ev_list are not preserved.
  }

  /*
   * Return the largest vertex id.
   */
  int getLargestVertexId() const {
    return next_vertex_id-1;
  }

  /*
   * Return the number of vertex_objs in this graph.
   */
  int getVertexNum() const {
    return vertex_obj_collection.size();
  }

  /*
   * Return the number of edge_objs in this graph.
   */
  int getEdgeNum() const {
    return edge_obj_collection.size();
  }

  /*
   * Return the list of vertices in this graph.
   */
  const VertexList& vertices() const {
    return vertex_collection;
  }

  /*
   * Return the list of edges_for_decoration in this graph.
   */
  const EdgeList& edges() const {
    return edge_collection;
  }

  /*
   * Check whether the vertex of a given ID exists in this graph.
   */
  bool hasVertexId(int id) const {
    return vid_to_vertex.contains(id);
  }

  /*
   * Return the vertex of a given ID in this graph.
   */
  Vertex getVertexById(int id) const {
    return vid_to_vertex.at(id);
  }

  /*
   * Add a new vertex to this graph.
   *
   * x - the element to be stored in the new vertex.
   * Return the newly created vertex.
   */
  Vertex addVertex(V x) {
    vertex_obj_collection.emplace_back(x);
    VertexObj& vertex_obj = vertex_obj_collection.back();

    vertex_obj.vertex_id = next_vertex_id;
    vertex_obj.vertex_obj_collection_itr = vertex_obj_collection.end();
    vertex_obj.vertex_obj_collection_itr--;
    vertex_collection.emplace_back(vertex_obj);
    Vertex& vertex = vertex_collection.back();
    vertex_obj.vertex_collection_itr = vertex_collection.end();
    vertex_obj.vertex_collection_itr--;

    vid_to_vertex[next_vertex_id] = vertex;
    next_vertex_id++;

    return vertex;
  }

  /*
   * Remove a vertex from this graph. All edge_objs that contain
   * v as one of their vertex_objs are also removed.
   *
   * v - a vertex
   */
  void deleteVertex(Vertex v) {
    VertexObj& v_obj = *(v.vertex_obj_ptr);
    // delete the incident edge_objs and vertex_objs
    for(auto [ e2, v2 ] : v_obj.outgoing_ev_list) {
      EdgeObj& e2_obj = *(e2.edge_obj_ptr);
      VertexObj& v2_obj = *(v2.vertex_obj_ptr);
      v2_obj.incoming_ev_list.erase(e2_obj.dest_incoming_ev_itr);
      edge_collection.erase(e2_obj.edge_collection_itr);
      edge_obj_collection.erase(e2_obj.edge_obj_collection_itr);
    }
    for(auto [ e2, v2 ] : v_obj.incoming_ev_list) {
      EdgeObj& e2_obj = *(e2.edge_obj_ptr);
      VertexObj& v2_obj = *(v2.vertex_obj_ptr);
      v2_obj.outgoing_ev_list.erase(e2_obj.orig_outgoing_ev_itr);
      edge_collection.erase(e2_obj.edge_collection_itr);
      edge_obj_collection.erase(e2_obj.edge_obj_collection_itr);
    }
    // actually delete the vertex
    vid_to_vertex.erase(v_obj.vertex_id);
    vertex_collection.erase(v_obj.vertex_collection_itr);
    vertex_obj_collection.erase(v_obj.vertex_obj_collection_itr);  // must be the last statement here, as it deletes *v
  }


  /*
   * Check whether an edge between v1 and v2 exists.
   *
   * v1 - the vertex at the origin
   * v2 - the vertex at the destination
   */
  bool hasEdge(Vertex v1, Vertex v2) const {
    for(auto [ e3, v3 ] : v1.vertex_obj_ptr->outgoing_ev_list) {
      if (v3 == v2) return true;
    }
    return false;
  }


  /*
   * Return an edge between v1 and v2 exists if exist
   *
   * v1 - the vertex at the origin
   * v2 - the vertex at the destination
   */
  Edge getEdge(Vertex v1, Vertex v2) const {
    for(auto [ e3, v3 ] : v1.vertex_obj_ptr->outgoing_ev_list) {
      if (v3 == v2) return e3;
    }
    return Edge();
  }

  /*
   * Add a new edge to this graph. Throw an exception
   * if an edge has already existed between v1 and v2.
   *
   * v1 - the vertex at the origin
   * v2 - the vertex at the destination
   * x  - the element to be stored in the new edge.
   * Return the newly created edge.
   */
  Edge addEdge(Vertex v1, Vertex v2, E x) {
    // check for duplication
    for(auto [ e3, v3 ] : v1.vertex_obj_ptr->outgoing_ev_list) {
      if (v3 == v2) throw std::runtime_error("Error: AdjListGraph<V,E>::addActiveEdge(): edge_obj_ptr has already existed.");
    }
    for(auto [ e3, v3 ] : v2.vertex_obj_ptr->incoming_ev_list) {
      if (v3 == v1) throw std::runtime_error("Error: AdjListGraph<V,E>::addActiveEdge(): edge_obj_ptr has already existed.");
    }
    // create the edge_obj_ptr
    edge_obj_collection.emplace_back(x);
    EdgeObj& edge_obj = edge_obj_collection.back();
    edge_obj.orig_vertex = v1;
    edge_obj.dest_vertex = v2;
    edge_obj.edge_obj_collection_itr = edge_obj_collection.end();
    edge_obj.edge_obj_collection_itr--;
    edge_collection.emplace_back(edge_obj);
    Edge& edge = edge_collection.back();
    edge_obj.edge_collection_itr = edge_collection.end();
    edge_obj.edge_collection_itr--;
    // update the [ edge_obj_ptr, vertex ] references in v1 and v2
    VertexObj& v1_obj = *(v1.vertex_obj_ptr);
    VertexObj& v2_obj = *(v2.vertex_obj_ptr);
    v1_obj.outgoing_ev_list.emplace_back(edge, v2);
    edge_obj.orig_outgoing_ev_itr = v1_obj.outgoing_ev_list.end();
    edge_obj.orig_outgoing_ev_itr--;
    v2_obj.incoming_ev_list.emplace_back(edge, v1);
    edge_obj.dest_incoming_ev_itr = v2_obj.incoming_ev_list.end();
    edge_obj.dest_incoming_ev_itr--;

    return edge;
  }

  /*
   * Remove an edge from this graph.
   *
   * edge - an edge
   */
  void deleteEdge(Edge edge) {
    EdgeObj& edge_obj = *(edge.edge_obj_ptr);
    Vertex v1 = edge.origVertex();
    Vertex v2 = edge.destVertex();
    v1.vertex_obj_ptr->outgoing_ev_list.erase(edge_obj.orig_outgoing_ev_itr);
    v2.vertex_obj_ptr->incoming_ev_list.erase(edge_obj.dest_incoming_ev_itr);
    // actually delete the edge
    edge_collection.erase(edge_obj.edge_collection_itr);
    edge_obj_collection.erase(edge_obj.edge_obj_collection_itr);   // must be the last statement here, as it deletes *edge
  }

  // ------ debug ------

  void debug() {
    std::cout << "------ VertexObj#=" << vertex_obj_collection.size() << " EdgeObj#=" << edge_obj_collection.size() << " ------" << std::endl;
    assert(vid_to_vertex.size() == vertex_obj_collection.size() && "vid_to_vertex and vertex_obj_collection have a different size");
    for(auto& v : vertex_collection) {
      std::cout << "Vertex " << v.id() << "[" << v.get() << "]:";
      assert(getVertexById(v.id()) == v && "getVertexById(v.id()) does not match with v");
      for(auto [ edge, v3 ] : v.outgoingEVList()) {
        std::cout << " (" << edge.origVertex().id() << "," << edge.destVertex().id() << ")[" << edge.get() << "]";
        assert(v3 == edge.destVertex() && "v3 != edge.destVertex()");
      }
      for(auto [ edge, v3 ] : v.incomingEVList()) {
        std::cout << " (" << edge.origVertex().id() << "," << edge.destVertex().id() << ")[" << edge.get() << "]";
        assert(v3 == edge.origVertex() && "v3 != edge.origVertex()");
      }
      std::cout << std::endl;
    }
    std::cout << "Edges:";
    for(auto& e : edge_collection) {
      std::cout << " (" << e.origVertex().id() << "," << e.destVertex().id() << ")" << "[" << e.get() << "]";
      EdgeObj& e_obj = *(e.edge_obj_ptr);
      assert(e == e_obj.orig_outgoing_ev_itr->first && "e != e_obj.orig_outgoing_ev_itr->first");
      assert(e == e_obj.dest_incoming_ev_itr->first && "e != e_obj.dest_incoming_ev_itr->first");
      assert(e.destVertex() == e_obj.orig_outgoing_ev_itr->second && "e.destVertex() == e_obj.orig_outgoing_ev_itr->second");
      assert(e.origVertex() == e_obj.dest_incoming_ev_itr->second && "e.origVertex() == e_obj.dest_incoming_ev_itr->second");
    }
    std::cout << std::endl;
  }

private:

  // ---------------------------------------------------------------------------------
  // Collections
  // ---------------------------------------------------------------------------------

  int next_vertex_id = 0;
  VertexObjList vertex_obj_collection;
  VertexList vertex_collection;
  EdgeObjList edge_obj_collection;
  EdgeList edge_collection;
  std::unordered_map<int,Vertex> vid_to_vertex;

};


/*
 * The default "nullptr" to a vertex
 */

template<typename V, typename E>
const typename AdjListGraph<V,E>::Vertex AdjListGraph<V,E>::NULL_VERTEX{};


// ---------------------------------------------------------------------------------
//    Dijkstra's Algorithm
// ---------------------------------------------------------------------------------

template<typename V, typename E>
class DijkstraAlgm {
  using Vertex = typename AdjListGraph<V,E>::Vertex;
  using VertexList = typename AdjListGraph<V,E>::VertexList;

  const AdjListGraph<V,E>& graph;
  int n;
  const Vertex start_vertex;
  const Vertex goal_vertex;
  int source_id;
  int goal_id;

  std::vector<E> dist;
  std::vector<int> prev;

public:

  DijkstraAlgm(const AdjListGraph<V,E>& graph, const Vertex start_vertex, const Vertex goal_vertex = AdjListGraph<V,E>::NULL_VERTEX) :
      graph{graph}, n{graph.getVertexNum()}, start_vertex{start_vertex}, goal_vertex{goal_vertex}
  {
    source_id = start_vertex.id();
    goal_id = (goal_vertex == AdjListGraph<V,E>::NULL_VERTEX) ? (-1) : goal_vertex.id();

    dist.resize(n, std::numeric_limits<E>::max());
    prev.resize(n, n);    // n means not visited (i.e., no previous vertex)

    calc();
  }

  const Vertex getStartVertex() const { return start_vertex; }

  bool hasGoalVertex() const { return goal_vertex != AdjListGraph<V,E>::NULL_VERTEX; }

  const Vertex getGoalVertex() const { return goal_vertex; }

  VertexList getPath(const Vertex new_goal_vertex = AdjListGraph<V,E>::NULL_VERTEX) {
    int new_goal_id = (new_goal_vertex == AdjListGraph<V,E>::NULL_VERTEX) ? goal_id : (new_goal_vertex.id());
    assert(new_goal_id >= 0 && "Error in DijkstraSolution::getPath(): new_goal_id < 0");

    VertexList path;
    int vid = new_goal_id;
    while (true) {
      path.push_front(graph.getVertexById(vid));
      vid = prev[vid];
      if (vid >= n) {  // the algorithm hasn't visited vid
        return VertexList{};  // just return an empty list to signify no solution
      } else if (vid == source_id) {
        assert(prev[vid] < 0 && "Error in DijkstraSolution::getPath(): prev[source_id] >= 0");  // a valid path must start with source_id
        return path;
      }
    }
  };

  E getCost(const Vertex new_goal_vertex = AdjListGraph<V,E>::NULL_VERTEX) {
    int new_goal_id = (new_goal_vertex == AdjListGraph<V,E>::NULL_VERTEX) ? goal_id : (new_goal_vertex.id());
    assert(new_goal_id >= 0 && "Error in DijkstraSolution::getCost(): new_goal_id < 0");
    return dist[new_goal_id];
  }

private:

  void calc() {
    // initialize the matrices;
    std::vector<bool> visited(n, false);
    dist[source_id] = E{};  // zero value
    prev[source_id] = -1;   // -1 means source_id is the root

    // set up priority queue
    auto cmp = [](std::pair<int, E> &p1, std::pair<int, E> &p2) { return p1.second > p2.second; };
    std::priority_queue<std::pair<int, E>, std::vector<std::pair<int, E>>, decltype(cmp)> queue{cmp};
    queue.push( {source_id, dist[source_id] } );

    // main loop
    while (!queue.empty()) {
      auto [ vid1, shortest_dist ] = queue.top();
      queue.pop();
      if (!visited[vid1]) {
        if (vid1 == goal_id) break;   // if (goal_id == -1) this check is always skipped.
        visited[vid1] = true;
        Vertex vertex1 = graph.getVertexById(vid1);
        for(auto [ edge, vertex2 ] : vertex1.outgoingEVList()) {
          int vid2 = vertex2.id();
          if (!visited[vid2]) {
            E new_cost = shortest_dist + edge.get();
            if (new_cost < dist[vid2]) {
              dist[vid2] = new_cost;
              prev[vid2] = vid1;
              queue.push({vid2, new_cost});
            }
          }
        }
      }  // else vid1 has been visited. Do nothing
    }
  }

};


// ---------------------------------------------------------------------------------
//    Floyd-Warshall Algorithm
// ---------------------------------------------------------------------------------

template<typename V, typename E>
class FloydWarshallAlgm {
  using Vertex = typename AdjListGraph<V, E>::Vertex;
  using VertexList = typename AdjListGraph<V, E>::VertexList;

  AdjListGraph<V, E> &graph;
  int n;

  std::vector<std::vector<E>> dist;
  std::vector<std::vector<int>> next;
  std::unordered_map<int,int> vertex_id_to_i;
  std::vector<int> i_to_vertex_id;


public:

  FloydWarshallAlgm(AdjListGraph<V, E> &graph) :
      graph{graph}, n{graph.getVertexNum()}
  {
    // Create the matrices
    dist.resize(n);
    next.resize(n);
    i_to_vertex_id.resize(n);
    for(int i=0; i<n; i++) {
      dist[i].resize(n, std::numeric_limits<E>::max());
      next[i].resize(n, -1);
    }
    // Find vertex_id to i maps
    int i=0;
    for(auto vertex : graph.vertices()) {
      int vid = vertex.id();
      vertex_id_to_i[vid] = i;
      i_to_vertex_id[i] = vid;
      i++;
    }
    // Initialize the matrices
    for(auto edge : graph.edges()) {
      int i1 = vertex_id_to_i[edge.origVertex().id()];
      int i2 = vertex_id_to_i[edge.destVertex().id()];
      dist[i1][i2] = edge.get();
      next[i1][i2] = i2;
    }
    for(int i=0; i<n; i++) {
      dist[i][i] = E{};
    }
    // Main loop
    for(int k=0; k<n; k++) {
      for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
          if (dist[i][k]!=std::numeric_limits<E>::max() && dist[k][j]!=std::numeric_limits<E>::max()) {
            if (dist[i][k] + dist[k][j] < dist[i][j]) {
              dist[i][j] = dist[i][k] + dist[k][j];
              next[i][j] = next[i][k];
            }
          }
        }
      }
    }
  }

  VertexList getPath(const Vertex orig_vertex, const Vertex dest_vertex) {
    VertexList path;

    int i1 = vertex_id_to_i[orig_vertex.id()];
    int i2 = vertex_id_to_i[dest_vertex.id()];
    if (next[i1][i2] < 0) return path;    // return an empty path

    while(i1 != i2) {
      i1 = next[i1][i2];
      path.push_back(graph.getVertexById(i_to_vertex_id[i1]));
    }
    return path;
  };

  E getCost(const Vertex orig_vertex, const Vertex dest_vertex) {
    int i1 = vertex_id_to_i[orig_vertex.id()];
    int i2 = vertex_id_to_i[dest_vertex.id()];
    return dist[i1][i2];
  }

  void print() {
    std::cout << "------ dist ------" << std::endl;
    for(int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        if (dist[i][j]==std::numeric_limits<E>::max()) {
          std::cout << "XX ";
        } else {
          std::cout << dist[i][j] << " ";
        }
      }
      std::cout << std::endl;
    }
    std::cout << "------ next ------" << std::endl;
    for(int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        std::cout << next[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }

};


// ---------------------------------------------------------------------------------
//    Testing
// ---------------------------------------------------------------------------------


void test_graph_1();

void test_graph_2();



#endif //SIM_GRAPH_H


