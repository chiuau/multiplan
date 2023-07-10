#include <iostream>
#include <iomanip>
#include <random>
#include <set>
#include <unordered_map>
#include <limits>
#include <chrono>
#include <numeric>
#include <optional>

#include "shared.h"
#include "graph.h"
#include "formation_graph.h"

#define USE_CACHE true


/**************************************************************************************
 *   These are the definitions of some helper functions for managing formation graphs
 ************************************************************************************** */


std::ostream& operator<<(std::ostream& out, const FormationGraph& graph) {
  for(int vid=0; vid<graph.getVertexNum(); vid++) {
    auto vertex = graph.getVertexById(vid);
    out << "Vertex " << vertex.id() << ":";
    for (unsigned int a1_id = 0; a1_id < vertex.get().size(); a1_id++) {
      auto a1 = vertex.get()[a1_id];
      out << " " << std::setw(2) << a1;
    }
    out << std::endl;
  }

  for(int vid1=0; vid1<graph.getVertexNum(); vid1++) {
    for(int vid2=0; vid2<graph.getVertexNum(); vid2++) {
      auto vertex1 = graph.getVertexById(vid1);
      auto vertex2 = graph.getVertexById(vid2);
      if (graph.hasEdge(vertex1, vertex2)) {
        out << "Edge: Vertex " << vertex1.id() << " -> " << "Vertex " << vertex2.id() << std::endl;
        auto edge = graph.getEdge(vertex1, vertex2);
        for (unsigned int a1_id = 0; a1_id < vertex1.get().size(); a1_id++) {
          for (unsigned int a2_id = 0; a2_id < vertex2.get().size(); a2_id++) {
            out << " " << std::setw(2) << edge.get()[a1_id][a2_id];
          }
          out << std::endl;
        }
      }
    }
  }
  return out;
}


std::ostream& operator<<(std::ostream& out, const FormationGraphResult& graph_result) {
  std::cout << "[";
  for(auto a : graph_result.assignment) {
    std::cout << " " << a;
  }
  std::cout << " ] (" << graph_result.time << ")";
  return out;
}



std::ostream& operator<<(std::ostream& out, const FormationGraphDynProgResult& graph_result) {
  std::cout << "[";
  for(auto& as : graph_result.assignment) {
    std::cout << "[";
    for(auto a : as) {
      std::cout << " " << a;
    }
    std::cout << " ]";
  }
  std::cout << " ] (" << graph_result.time << ")";
  return out;
}


/**************************************************************************************
 *  The definitions of the member functions in FormationGraphDynProgResult
 ************************************************************************************** */


FormationGraphResult FormationGraphDynProgResult::toFormationGraphResultByFirst(const FormationGraph& graph) const {
  FormationGraphResult result(time, graph.getVertexNum());
  for(int vid=0; vid<graph.getVertexNum(); vid++) {
    assert(!assignment[vid].empty());
    result.assign(vid, assignment[vid][0]);
  }
  return result;
}


FormationGraph FormationGraphDynProgResult::makeReducedGraph(const FormationGraph& old_graph) const {
  FormationGraph new_graph;
  std::vector<std::vector<int>> new2old_assignment(old_graph.getVertexNum());
  std::vector<int> new_vertex_action_num(old_graph.getVertexNum());

  for(int vid=0; vid < old_graph.getVertexNum(); vid++) {
    assert(!assignment[vid].empty());
    auto& old_aids = old_graph.getVertexById(vid).get();
    FormationActions actions;
    for(auto aid : assignment[vid]) {
      actions.push_back(old_aids[aid]);
      new2old_assignment[vid].push_back(aid);
    }
    new_vertex_action_num[vid] = assignment[vid].size();
    auto new_vertex = new_graph.addVertex(actions);
    assert(new_vertex.id() == vid);
  }

  for(auto v1 : old_graph.vertices()) {
    for(auto v2 : old_graph.vertices()) {
      if (old_graph.hasEdge(v1, v2)) {
        auto old_deltas = old_graph.getEdge(v1, v2).get();

        FormationDeltas new_deltas;
        for(int i=0; i<new_vertex_action_num[v1.id()]; i++) {
          new_deltas.emplace_back(new_vertex_action_num[v2.id()], -1);
        }

        for(int new_aid_1=0; new_aid_1 < new_vertex_action_num[v1.id()]; new_aid_1++) {
          for(int new_aid_2=0; new_aid_2 < new_vertex_action_num[v2.id()]; new_aid_2++) {
            new_deltas[new_aid_1][new_aid_2] = old_deltas[new2old_assignment[v1.id()][new_aid_1]][new2old_assignment[v2.id()][new_aid_2]];
          }
        }

        new_graph.addEdge(new_graph.getVertexById(v1.id()), new_graph.getVertexById(v2.id()), new_deltas);
      }
    }
  }

  return new_graph;
}


/**************************************************************************************
 *   These are the definitions of some helper functions for calculating the times
 ************************************************************************************** */


void calc_formation_graph_exact_time_dfs(FormationSchedule& plan,
                                         const FormationGraph& graph,
                                         const std::vector<int>& assignment,
                                         const FormationGraph::Vertex& vertex,
                                         std::vector<FormationTime>& delta_list)
{
  FormationTime sum = std::accumulate(delta_list.begin(), delta_list.end(), 0) + vertex.get()[assignment[vertex.id()]];

  if (sum > plan[vertex.id()]) {
    plan[vertex.id()] = sum;
  }

  for (auto[edge, next_vertex]: vertex.outgoingEVList()) {
    auto delta = edge.get()[assignment[vertex.id()]][assignment[next_vertex.id()]];
    delta_list.push_back(delta);
    calc_formation_graph_exact_time_dfs(plan, graph, assignment, next_vertex, delta_list);
    delta_list.pop_back();
  }
}


FormationSchedule calc_formation_graph_end_time_schedule(const FormationGraph& graph, const FormationGraphResult& result, int init_vertex_id) {
  FormationSchedule plan(graph.getVertexNum(), 0);
  std::vector<FormationTime> delta_list;

  calc_formation_graph_exact_time_dfs(plan, graph, result.getAssignment(), graph.getVertexById(init_vertex_id), delta_list);

  return plan;
}


FormationTime calc_formation_graph_end_time_schedule_max_time(const FormationSchedule& end_time_schedule) {
  FormationTime max_time = 0;
  for(auto t : end_time_schedule) {
    if (t > max_time) max_time = t;
  }
  return max_time;
}


FormationSchedule calc_formation_graph_schedule(const FormationGraph& graph, const FormationGraphResult& result, const FormationSchedule& end_time_schedule) {
  FormationSchedule schedule = end_time_schedule;

  for(signed int i=0; i < end_time_schedule.size(); i++) {
    schedule[i] -= graph.getVertexById(i).get()[result.getActionId(i)];
  }

  return schedule;
}


void advanceFormationGraphSchedule(FormationSchedule& schedule, FormationTime time) {
  for(unsigned int i=0; i < schedule.size(); i++) {
    schedule[i] += time;
  }
}


/**************************************************************************************
 *  The brute force algorithm
 *
 *  This is used to check the correctness of the graph search algorithm and the dynamic programming algorithm.
 ************************************************************************************** */


void formation_graph_brute_force_eval_dfs(FormationTime& max_value, const FormationGraph& graph, const std::unordered_map<int,int>& action_combination, const FormationGraph::Vertex& vertex, FormationTime sum, bool isDebug= false) {
  int a1_id = action_combination.at(vertex.id());
  int cost = sum + vertex.get()[a1_id];
  if (cost > max_value) { max_value = cost; }

  for (auto[edge, vertex2]: vertex.outgoingEVList()) {
    int a2_id = action_combination.at(vertex2.id());
    formation_graph_brute_force_eval_dfs(max_value, graph, action_combination, vertex2, sum + edge.get()[a1_id][a2_id],
                                         isDebug);
  }
}


FormationTime formation_graph_brute_force_eval(const FormationGraph& graph, const std::unordered_map<int,int>& action_combination, const FormationGraph::Vertex init_vertex, bool isDebug = false) {
  FormationTime max_value = 0;
  formation_graph_brute_force_eval_dfs(max_value, graph, action_combination, init_vertex, 0, isDebug);
  return max_value;
}


void formation_graph_brute_force_dfs(FormationTime& min_value, std::vector<int>& min_action_combination, std::unordered_map<int,int>& action_combination, const FormationGraph& graph,
                                     FormationGraph::VertexList::const_iterator& vertex_itr, const FormationGraph::Vertex& init_vertex)
{
  if (action_combination.size() == static_cast<unsigned int>(graph.getVertexNum())) {
    auto v = formation_graph_brute_force_eval(graph, action_combination, init_vertex);
    if (v < min_value) {
      min_value = v;
      for(int vid=0; vid<graph.getVertexNum(); vid++) {
        min_action_combination[vid] = action_combination[vid];
      }
    }
  } else {
    int id = vertex_itr->id();
    int i=0;
    for(auto v : vertex_itr->get()) {
      action_combination[id] = i;
      vertex_itr++;
      formation_graph_brute_force_dfs(min_value, min_action_combination, action_combination, graph, vertex_itr, init_vertex);
      vertex_itr--;
      action_combination.erase(id);
      i++;
    }
  }
}


FormationGraphResult formation_graph_brute_force(const FormationGraph& graph, int init_vertex_id) {
  std::unordered_map<int,int> action_combination;
  std::vector<int> min_action_combination(graph.getVertexNum());
  FormationTime min_value = std::numeric_limits<FormationTime>::max();
  auto vertex_itr = graph.vertices().begin();
  formation_graph_brute_force_dfs(min_value, min_action_combination, action_combination, graph, vertex_itr, graph.getVertexById(init_vertex_id));
  return FormationGraphResult(min_value, min_action_combination);
}


/**************************************************************************************
 *   Formation Graph Dynamic Programming
 *
 *   This implementation is slightly different from the Algorithm 1 in the paper.
 *   In the paper, Algorithm 1 computes every entry in the P table, which is edge_stored_values
 *   in the function below. However, the following implementation avoid computing every
 *   entry in edge_stored_values by doing a depth-first search, which only compute the entries
 *   that are relevant to the value of the entry of the root node in edge_stored_values.
 *
 ************************************************************************************** */

/**************************************************************************************
 *   formation_graph_dynamic_programming_dfs()
 *
 *   The main loop of the dfs of the dynamic programming method. This computes Equation 5 in the paper directly.
 *
 *     graph               - the formation graph
 *     edge_stored_values  - the table for dynamic programming
 *     vertex1             - the current vertex
 *     last_vertex         - the parent of vertex1
 *     last_a_id           - the action id chosen for last_vertex
 *     init_vertex_id      - the vertex id of the root vertex
 ************************************************************************************** */


FormationGraphDynProgResult formation_graph_dynamic_programming_dfs(const FormationGraph& graph,
                                                                    std::vector<std::vector<std::vector<FormationGraphDynProgResult>>>& edge_stored_values,
                                                                    const FormationGraph::Vertex& vertex1,
                                                                    const FormationGraph::Vertex& last_vertex,
                                                                    int last_a_id,
                                                                    int init_vertex_id)
{
  // initialize the min value; the initial min value is infinity (i.e., no value)
  FormationGraphDynProgResult min_value(std::numeric_limits<FormationTime>::max(), graph.getVertexNum());

  for(unsigned int a1_id = 0; a1_id < vertex1.get().size(); a1_id++) {
    auto a1 = vertex1.get()[a1_id];

    // initialize the max value; the initial max value is a1
    FormationGraphDynProgResult max_value(a1, graph.getVertexNum());

    // for each children of vertex1,
    for (auto[edge, vertex2]: vertex1.outgoingEVList()) {
      FormationGraphDynProgResult& v_of_max = edge_stored_values[vertex1.id()][vertex2.id()][a1_id];
      if (v_of_max.getTime() < 0) {  // haven't compute and store the value yet
        // compute and store the value in v_of_max as well as edge_stored_values
        v_of_max = formation_graph_dynamic_programming_dfs(graph, edge_stored_values, vertex2, vertex1, a1_id, init_vertex_id);
        assert(v_of_max.getTime() >= 0);
      }
      // updating the max value
      if (v_of_max.getTime() > max_value.getTime()) {
        max_value.setTime(v_of_max.getTime());
      }
      max_value.assign(v_of_max.getAssignment());  // need to rename the chosen action of the current max value
    }

    // add delta to the max value
    assert((vertex1.id()!=init_vertex_id)?(last_a_id>=0):(last_a_id<0));
    int delta = (vertex1.id()==init_vertex_id)?0:(graph.getEdge(last_vertex, vertex1).get()[last_a_id][a1_id]);
    max_value.addTime(delta);

    // updating the min value
    if (max_value.getTime() < min_value.getTime()) {
      min_value = max_value;
      min_value.assign(vertex1.id(), a1_id);
    }
  }

  return min_value;
}


/**************************************************************************************
 *   formation_graph_dynamic_programming()
 *
 *   The main function of the dynamic programming algorithm.
 *
 *     graph          - the formation graph
 *     init_vertex_id - the vertex id of the root vertex
 ************************************************************************************** */

FormationGraphDynProgResult formation_graph_dynamic_programming(const FormationGraph& graph, int init_vertex_id) {

  // The table for dynamic programming.  There is one entry for each action a1 on an edge (v1, v2).
  // The entry edge_stored_values(v1_id, v2_id, a1_id) stores a solution
  std::vector<std::vector<std::vector<FormationGraphDynProgResult>>> edge_stored_values;

  // initialize the table
  for (int v1_id = 0; v1_id < graph.getVertexNum(); v1_id++) {
    auto v1 = graph.getVertexById(v1_id);
    edge_stored_values.emplace_back(graph.getVertexNum());
    for (int v2_id = 0; v2_id < graph.getVertexNum(); v2_id++) {
      auto v2 = graph.getVertexById(v2_id);
      if (graph.hasEdge(v1, v2)) {
        for (unsigned int a1_id = 0; a1_id < v1.get().size(); a1_id++) {
          edge_stored_values[v1_id][v2_id].emplace_back(-1, graph.getVertexNum());   // -1 means no solution yet
        }
      }
    }
  }

  // start running the algorithm
  auto result = formation_graph_dynamic_programming_dfs(graph,
                                                        edge_stored_values,
                                                        graph.getVertexById(init_vertex_id),
                                                        graph.getVertexById(init_vertex_id), // just a dummy
                                                        -1,   // the root has no parent, and hence no action id
                                                        init_vertex_id); // just a dummy

  return result;
}


/**************************************************************************************
 *  The Graph Search Algorithm
 *
 *  This is a completer search algorithm that is roughly three times faster than the brute force method.
 ************************************************************************************** */

/**************************************************************************************
 *  This is the cache being used in the graph search algorithm
 ************************************************************************************** */


class FormationGraphSearchCache {

  const FormationGraph& graph;
  const int init_vertex_id;
  const std::vector<std::vector<int>>& common_prefix_last_vid;

  std::vector<std::set<int>> free_vertex_ids;
  std::vector<std::vector<int>> inverted_free_vertex_ids;

  std::vector<int> assignment;
  std::vector<std::optional<FormationGraphResult>> cache;

public:

  FormationGraphSearchCache(const FormationGraph& graph, int init_vertex_id, const std::vector<std::vector<int>>& common_prefix_last_vid) :
      graph(graph),
      init_vertex_id(init_vertex_id),
      common_prefix_last_vid(common_prefix_last_vid),
      free_vertex_ids(graph.getVertexNum()),
      inverted_free_vertex_ids(graph.getVertexNum()),
      assignment(graph.getVertexNum(), -1),
      cache(graph.getVertexNum(), std::nullopt)
  {
    calcFreeVertexIds();
    calcInvertedFreeVertexIds();
  }

  void assign(int vid, int aid) {
    assert(assignment[vid] < 0);
    assignment[vid] = aid;
  }

  void unassign(int vid) {
    assert(assignment[vid] >= 0);
    assignment[vid] = -1;
    for(auto affected_vid : inverted_free_vertex_ids[vid]) {
      cache[affected_vid] = std::nullopt;
    }
  }

  bool isCachable(int vid, const std::vector<int>& current_assignment) {
    for(auto affected_vid : free_vertex_ids[vid]) {
      if (assignment[affected_vid] < 0) return false;
      if (assignment[affected_vid] != current_assignment[affected_vid]) return false;
    }
    return true;
  }

  bool isCached(int vid) const {
    #if USE_CACHE
      return cache[vid].has_value();
    #else
      return false;
    #endif
  }

  FormationGraphResult get(int vid) const {
    return *cache[vid];
  }

  void set(int vid, const FormationGraphResult& result) {
    #if USE_CACHE
      cache[vid] = result;
    #endif
  }


private:

  void calcFreeVertexIds() {
    std::set<int> visited_vertex;
    calcFreeVertexIdsSelf(graph.getVertexById(init_vertex_id), visited_vertex);   // ***

    for(int vid=0; vid<graph.getVertexNum(); vid++) {
      for(auto prefix_last_vid : common_prefix_last_vid[vid]) {
        visited_vertex.clear();
        std::vector<bool> mark(graph.getVertexNum(), false);
        calcFreeVertexIdsMark(graph.getVertexById(vid), visited_vertex, prefix_last_vid, mark);
        visited_vertex.clear();
        calcFreeVertexIdsBackwardSearch(graph.getVertexById(prefix_last_vid), visited_vertex, prefix_last_vid, mark);
      }
    }
  }

  void calcFreeVertexIdsSelf(const FormationGraph::Vertex vertex, std::set<int>& visited_vertex) {
    if (!visited_vertex.contains(vertex.id())) {
      free_vertex_ids[vertex.id()].insert(vertex.id());
      visited_vertex.insert(vertex.id());
      for (auto[edge, next_vertex]: vertex.outgoingEVList()) {
        calcFreeVertexIdsSelf(next_vertex, visited_vertex);
      }
    }
  }

  void calcFreeVertexIdsMark(const FormationGraph::Vertex vertex, std::set<int>& visited_vertex, int prefix_last_vid, std::vector<bool>& mark) {
    mark[vertex.id()] = true;
    visited_vertex.insert(vertex.id());

    if (vertex.id() != prefix_last_vid) {
      for (auto[edge, next_vertex]: vertex.outgoingEVList()) {
        if (!visited_vertex.contains(next_vertex.id())) {
          calcFreeVertexIdsMark(next_vertex, visited_vertex, prefix_last_vid, mark);
        }
      }
    }
  }

  void calcFreeVertexIdsBackwardSearch(const FormationGraph::Vertex vertex, std::set<int>& visited_vertex, int prefix_last_vid, std::vector<bool>& mark) {
    assert(mark[vertex.id()]);
    free_vertex_ids[vertex.id()].insert(prefix_last_vid);
    visited_vertex.insert(vertex.id());
    for (auto[edge, previous_vertex]: vertex.incomingEVList()) {
      if (!visited_vertex.contains(previous_vertex.id()) && mark[previous_vertex.id()]) {
        calcFreeVertexIdsBackwardSearch(previous_vertex, visited_vertex, prefix_last_vid, mark);
      }
    }
  }

  void calcInvertedFreeVertexIds() {
    for(int vid1=0 ; vid1 < graph.getVertexNum(); vid1++) {
      for(auto vid2 : free_vertex_ids[vid1]) {
        inverted_free_vertex_ids[vid2].push_back(vid1);
      }
    }
  }

};


void find_common_prefix_in_formation_graph_dfs(std::list<int>& reversed_prefix, const FormationGraph& graph, const FormationGraph::Vertex& vertex, std::set<int>& visited, int init_vertex_id) {
  auto find_ptr = std::find(reversed_prefix.begin(), reversed_prefix.end(), vertex.id());

  if (find_ptr == reversed_prefix.end()) {  // vertex not on the reserved prefix
    if (!visited.contains(vertex.id())) {
      if (reversed_prefix.empty() || reversed_prefix.back() != init_vertex_id) {  // the first prefix has not been established yet.
        reversed_prefix.push_back(vertex.id());
      }
      visited.insert(vertex.id());

      for (auto[edge, previous_vertex]: vertex.incomingEVList()) {
        find_common_prefix_in_formation_graph_dfs(reversed_prefix, graph, previous_vertex, visited, init_vertex_id);
      }
    } // else do nothing
  } else {  // vertex is on the reserved prefix
    // delete the prefix of the reversed_prefix, until *find_ptr
    reversed_prefix.erase(reversed_prefix.begin(), find_ptr);
    // just backtrack
  }
}


std::vector<int> find_common_prefix_in_formation_graph(const FormationGraph& graph, int init_vertex_id) {
  std::vector<int> common_prefix_last_vertex(graph.getVertexNum());

  for(int vid=0; vid<graph.getVertexNum(); vid++) {
    std::set<int> visited;
    std::list<int> reversed_prefix;
    find_common_prefix_in_formation_graph_dfs(reversed_prefix, graph, graph.getVertexById(vid), visited, init_vertex_id);
    common_prefix_last_vertex[vid] = reversed_prefix.front();
  }

  return common_prefix_last_vertex;
}


bool advance_vector_counter(std::vector<int>& counter, const std::vector<int>& counter_size) {
  unsigned int i = 0;
  while(i < counter.size()) {
    counter[i]++;
    if (counter[i] < counter_size[i]) return true;
    counter[i]=0;
    i++;
  }
  return false;
}


FormationGraphResult formation_graph_search_dfs(const FormationGraph& graph,
                                                const std::vector<std::vector<int>>& common_prefix_last_vid,
                                                const FormationGraph::Vertex& current_vertex,
                                                const FormationGraph::Vertex& last_vertex,
                                                std::vector<int>& assignment,
                                                FormationGraphSearchCache& cache,
                                                int init_vertex_id,
                                                double time_limit,
                                                const std::chrono::time_point<std::chrono::high_resolution_clock>& start_time);


FormationGraphResult formation_graph_search_dfs_max(const FormationGraph& graph,
                                                    const std::vector<std::vector<int>>& common_prefix_last_vid,
                                                    const FormationGraph::Vertex& current_vertex,
                                                    const FormationGraph::Vertex& last_vertex,
                                                    std::vector<int>& assignment,
                                                    FormationGraphSearchCache& cache,
                                                    int init_vertex_id,
                                                    double time_limit,
                                                    const std::chrono::time_point<std::chrono::high_resolution_clock>& start_time)
{
  if (duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count() >= time_limit) {
    return FormationGraphResult();  // invalid value
  }

  // set the max value to a1
  int current_vid = current_vertex.id();
  int a1_id = assignment[current_vid];
  assert(a1_id >= 0);
  auto a1 = current_vertex.get()[a1_id];
  FormationGraphResult max_value(a1, graph.getVertexNum());

  if (common_prefix_last_vid[current_vid].empty()) {  // no min_max

    // assert(cache.isCachable(current_vertex.id(), assignment));

    if (cache.isCached(current_vertex.id())) {
      max_value = cache.get(current_vertex.id());
    } else {
      // compute the max value directly
      for (auto[edge, next_vertex]: current_vertex.outgoingEVList()) {
        auto v_of_max = formation_graph_search_dfs(graph, common_prefix_last_vid, next_vertex, current_vertex,
                                                   assignment, cache, init_vertex_id, time_limit, start_time);
        if (!v_of_max.isValid()) { return v_of_max; }

        if (v_of_max.getTime() > max_value.getTime()) {
          max_value.setTime(v_of_max.getTime());
        }
        // __line__();
        max_value.assign(v_of_max.getAssignment());
        // __line__();
      }
      cache.set(current_vertex.id(), max_value);
    }

  } else {  // have min_max

    FormationGraphResult min_value_2(std::numeric_limits<FormationTime>::max(), graph.getVertexNum());

    // prepare the counter
    std::vector<int> aid(common_prefix_last_vid[current_vid].size(), 0);
    std::vector<int> aid_size;
    for(auto vid : common_prefix_last_vid[current_vid]) {
      aid_size.push_back(graph.getVertexById(vid).get().size());
    }

    // start the min_max loop
    bool isContinue = true;
    while(isContinue) {
      // start evaluation

      // assign
      int j=0;
      for(auto vid: common_prefix_last_vid[current_vid]) {
        assert(assignment[vid] < 0);
        assignment[vid] = aid[j];
        cache.assign(vid, aid[j]);
        j++;
      }

      // assert(cache.isCachable(current_vertex.id(), assignment));  // must put this here after the assignment
      FormationGraphResult max_value_2(0, graph.getVertexNum());

      if (cache.isCached(current_vertex.id())) {
        max_value_2 = cache.get(current_vertex.id());
      } else {
        // compute the max value 2
        for (auto[edge, next_vertex]: current_vertex.outgoingEVList()) {
          auto v_of_max_2 = formation_graph_search_dfs(graph, common_prefix_last_vid, next_vertex, current_vertex,
                                                       assignment, cache, init_vertex_id, time_limit, start_time);
          if (!v_of_max_2.isValid()) { return v_of_max_2; }

          if (v_of_max_2.getTime() > max_value_2.getTime()) {
            max_value_2.setTime(v_of_max_2.getTime());
          }
          max_value_2.assign(v_of_max_2.getAssignment());
        }
        cache.set(current_vertex.id(), max_value_2);
      }

      // improve the min value 2
      if (max_value_2.getTime() < min_value_2.getTime()) {
        min_value_2 = max_value_2;
        // remember the current assignment for the min_max
        int j=0;
        for(auto vid: common_prefix_last_vid[current_vid]) {
          min_value_2.assign(vid, aid[j]);
          j++;
        }
      }

      // unassign
      for(auto vid: common_prefix_last_vid[current_vid]) {
        assignment[vid] = -1;
        cache.unassign(vid);
      }

      // next combination
      isContinue = advance_vector_counter(aid, aid_size);
    }

    // improve max_value
    if (min_value_2.getTime() > max_value.getTime()) {
      max_value.setTime(min_value_2.getTime());
    }
    max_value.assign(min_value_2.getAssignment());
  }

  return max_value;
}


FormationGraphResult formation_graph_search_dfs(const FormationGraph& graph,
                                                const std::vector<std::vector<int>>& common_prefix_last_vid,
                                                const FormationGraph::Vertex& current_vertex,
                                                const FormationGraph::Vertex& last_vertex,
                                                std::vector<int>& assignment,
                                                FormationGraphSearchCache& cache,
                                                int init_vertex_id,
                                                double time_limit,
                                                const std::chrono::time_point<std::chrono::high_resolution_clock>& start_time)
{
  if (duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count() >= time_limit) {
    return FormationGraphResult();  // invalid value
  }

  if (current_vertex.id() != init_vertex_id && assignment[current_vertex.id()] >= 0) {  // no min loop

    int last_a_id = assignment[last_vertex.id()];
    assert(last_a_id >= 0);
    int a1_id = assignment[current_vertex.id()];
    assert(a1_id >= 0);
    FormationTime delta = graph.getEdge(last_vertex, current_vertex).get()[last_a_id][a1_id];

    // assert(cache.isCachable(current_vertex.id(), assignment));

    FormationGraphResult max_value;
    if (cache.isCached(current_vertex.id())) {
      max_value = cache.get(current_vertex.id());
    } else {
      max_value = formation_graph_search_dfs_max(graph,
                                                 common_prefix_last_vid,
                                                 current_vertex,
                                                 last_vertex,
                                                 assignment,
                                                 cache,
                                                 init_vertex_id,
                                                 time_limit,
                                                 start_time);
      if (!max_value.isValid()) { return max_value; }

      cache.set(current_vertex.id(), max_value);
    }

    max_value.addTime(delta);
    return max_value;

  } else {  // have min loop
    FormationGraphResult min_value(std::numeric_limits<FormationTime>::max(), graph.getVertexNum());

    for(unsigned int a1_id = 0; a1_id < current_vertex.get().size(); a1_id++) {
      FormationTime delta = (current_vertex.id() == init_vertex_id)?
                            0:
                            (graph.getEdge(last_vertex, current_vertex).get()[assignment[last_vertex.id()]][a1_id]);
      assignment[current_vertex.id()] = a1_id;
      cache.assign(current_vertex.id(), a1_id);

      auto max_value = formation_graph_search_dfs_max(graph,
                                                      common_prefix_last_vid,
                                                      current_vertex,
                                                      last_vertex,
                                                      assignment,
                                                      cache,
                                                      init_vertex_id,
                                                      time_limit,
                                                      start_time);
      if (!max_value.isValid()) { return max_value; }

      assignment[current_vertex.id()] = -1;
      cache.unassign(current_vertex.id());

      max_value.addTime(delta);
      if (max_value.getTime() < min_value.getTime()) {
        min_value = max_value;
        min_value.assign(current_vertex.id(), a1_id);
      }
    }

    return min_value;
  }
}


/**************************************************************************************
 *  This is the main function of the graph search algorithm
 ************************************************************************************** */


FormationGraphResult formation_graph_search(const FormationGraph& graph, int init_vertex_id, double time_limit) {
  auto start_time = std::chrono::high_resolution_clock::now();

  auto common_prefix_last_vertex = find_common_prefix_in_formation_graph(graph, init_vertex_id);

  std::vector<std::vector<int>> common_prefix_last_vid(graph.getVertexNum());
  for(int vid=0; vid<graph.getVertexNum(); vid++) {
    if (common_prefix_last_vertex[vid] != vid) {  // vid is of cause a default last vertex of itself.
      common_prefix_last_vid[common_prefix_last_vertex[vid]].push_back(vid);
    }
  }

  FormationGraphSearchCache cache(graph, init_vertex_id, common_prefix_last_vid);
  std::vector<int> assignment(graph.getVertexNum(), -1);

  auto result = formation_graph_search_dfs(graph,
                                           common_prefix_last_vid,
                                           graph.getVertexById(init_vertex_id),
                                           graph.getVertexById(init_vertex_id), // just a dummy
                                           assignment,
                                           cache, // just a dummy
                                           init_vertex_id,
                                           time_limit,
                                           start_time);
  return result;
}



/**************************************************************************************
 *  Testing
 ************************************************************************************** */


FormationActions calcDefaultFormationActions(int n, std::mt19937& rng, std::uniform_int_distribution<std::mt19937::result_type>& rndgen_action_duration) {
  FormationActions actions;
  for(int i=0; i<n; i++) {
    actions.push_back(rndgen_action_duration(rng));
  }
  return actions;
}


FormationDeltas calcDefaultFormationDeltas(int n, std::mt19937& rng, std::uniform_int_distribution<std::mt19937::result_type>& rndgen_delta) {
  FormationDeltas deltas;
  for(int i=0; i<n; i++) {
    deltas.emplace_back();
    for(int j=0; j<n; j++) {
      deltas.back().push_back(rndgen_delta(rng));
    }
  }
  return deltas;
}


FormationGraph makeRandomGridFormationGraph(std::mt19937& rng, int size_x, int size_y, int action_num) {

  auto calcVertexId = [](int x, int y, int size_x) {
    return y * size_x + x;
  };

  std::uniform_int_distribution<std::mt19937::result_type> rndgen_action_duration(5, 10);
  std::uniform_int_distribution<std::mt19937::result_type> rndgen_delta(1, 4);

  FormationGraph graph;

  for(int y=0; y<size_y; y++) {
    for(int x=0; x<size_x; x++) {
      int id = calcVertexId(x, y, size_x);
      auto v = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
      assert(v.id() == id && "Error in makeRandomGridFormationGraph(): incorrect node id");
    }
  }

  for(int y=0; y<size_y-1; y++) {
    for(int x=0; x<size_x; x++) {
      int id1 = calcVertexId(x, y, size_x);
      int id2 = calcVertexId(x, y + 1, size_x);
      auto v1 = graph.getVertexById(id1);
      auto v2 = graph.getVertexById(id2);
      graph.addEdge(v1, v2, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
    }
  }

  for(int x=0; x<size_x-1; x++) {
    for(int y=0; y<size_y; y++) {
      int id1 = calcVertexId(x, y, size_x);
      int id2 = calcVertexId(x + 1, y, size_x);
      auto v1 = graph.getVertexById(id1);
      auto v2 = graph.getVertexById(id2);
      graph.addEdge(v1, v2, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
    }
  }

  return graph;
}


void pure_formation_graph_expr(const FormationGraph& graph, int init_vertex_id) {

  auto start_time_gg = std::chrono::high_resolution_clock::now();
  auto ans_gg = formation_graph_search(graph, init_vertex_id, EXPR_TIME_LIMIT);
  auto end_time_gg = std::chrono::high_resolution_clock::now();

  auto start_time_dp = std::chrono::high_resolution_clock::now();
  auto ans_dp_tmp = formation_graph_dynamic_programming(graph, init_vertex_id);
  auto end_time_dp_tmp = std::chrono::high_resolution_clock::now();
  auto graph_dp = ans_dp_tmp.makeReducedGraph(graph);
  auto ans_dp = formation_graph_search(graph_dp, init_vertex_id, EXPR_TIME_LIMIT);
  auto end_time_dp = std::chrono::high_resolution_clock::now();

  double exetime_gg = duration_cast<std::chrono::microseconds>(end_time_gg - start_time_gg).count() / 1000000.0;
  double exetime_dp_tmp = duration_cast<std::chrono::microseconds>(end_time_dp_tmp - start_time_dp).count() / 1000000.0;
  double exetime_dp = duration_cast<std::chrono::microseconds>(end_time_dp - start_time_dp).count() / 1000000.0;

  std::cout << "No of cars = " << graph.getVertexNum() << std::endl;
  std::cout << "GSearch makespan = " << (ans_gg.isValid()?(ans_gg.getTime() * EXPR_TIME_STEP):0.0) << std::endl;
  std::cout << "DynProg makespan = " << (ans_dp.isValid()?(ans_dp.getTime() * EXPR_TIME_STEP):0.0) << std::endl;
  std::cout << "GSearch execution time = " << exetime_gg << std::endl;
  std::cout << "DynProg execution time = " << exetime_dp << std::endl;
}


void pure_formation_graph_expr_1() {
  auto& rng = Shared::getInstance().getRng();
  auto graph = makeRandomGridFormationGraph(rng, 4, 5, 2);
  int init_vertex_id = 0;
  pure_formation_graph_expr(graph, init_vertex_id);
}


