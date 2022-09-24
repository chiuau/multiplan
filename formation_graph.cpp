//
// Created by Tsz-Chiu Au on 3/1/22.
//

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


// **************************************************************************************
//   Formation Graph
// **************************************************************************************

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


// **************************************************************************************
//   Formation Graph Brute Force
// **************************************************************************************


void formation_graph_brute_force_eval_dfs(FormationTime& max_value, const FormationGraph& graph, const std::unordered_map<int,int>& action_combination, const FormationGraph::Vertex& vertex, FormationTime sum, bool isDebug= false) {
  int a1_id = action_combination.at(vertex.id());
//  if (isDebug) __pp__("vid=", vertex.id(), " aid=", a1_id, " a=", vertex.get()[a1_id]);
//  if (isDebug) __pp__("sum = ", sum);
  int cost = sum + vertex.get()[a1_id];
  if (cost > max_value) { max_value = cost; }
//  if (isDebug) __pp__("  cost = ", cost);

  for (auto[edge, vertex2]: vertex.outgoingEVList()) {
    int a2_id = action_combination.at(vertex2.id());
//    if (isDebug) __pp__("edge.get()[", a1_id, "][", a2_id, "] = ", edge.get()[a1_id][a2_id]);
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
    // __pp__(v, " ", min_action_combination);  // debug
    if (v < min_value) {
      min_value = v;
      for(int vid=0; vid<graph.getVertexNum(); vid++) {
        min_action_combination[vid] = action_combination[vid];
      }
//      __pp__(min_value, " ", min_action_combination);  // debug
//      formation_graph_brute_force_eval(graph, action_combination, true);  // debug
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


// **************************************************************************************
//   Formation Graph Dynamic Programming
// **************************************************************************************


FormationGraphDynProgResult formation_graph_dynamic_programming_dfs(const FormationGraph& graph,
                                                                    std::vector<std::vector<std::vector<FormationGraphDynProgResult>>>& edge_stored_values,
                                                                    const FormationGraph::Vertex& vertex1,
                                                                    const FormationGraph::Vertex& last_vertex,
                                                                    int last_a_id,
                                                                    int init_vertex_id)
{
  FormationGraphDynProgResult min_value(std::numeric_limits<FormationTime>::max(), graph.getVertexNum());

  for(unsigned int a1_id = 0; a1_id < vertex1.get().size(); a1_id++) {
    auto a1 = vertex1.get()[a1_id];

    FormationGraphDynProgResult max_value(a1, graph.getVertexNum());

    for (auto[edge, vertex2]: vertex1.outgoingEVList()) {
      FormationGraphDynProgResult& v_of_max = edge_stored_values[vertex1.id()][vertex2.id()][a1_id];
      if (v_of_max.getTime() < 0) {
        v_of_max = formation_graph_dynamic_programming_dfs(graph, edge_stored_values, vertex2, vertex1, a1_id, init_vertex_id);
        assert(v_of_max.getTime() >= 0);
      }
      if (v_of_max.getTime() > max_value.getTime()) {
        max_value.setTime(v_of_max.getTime());
      }
      max_value.assign(v_of_max.getAssignment());
    }

//    assert((vertex1.id()>0)?(last_a_id>=0):(last_a_id<0));
//    int delta = (vertex1.id()==0)?0:(graph.getEdge(last_vertex, vertex1).get()[last_a_id][a1_id]);
    assert((vertex1.id()!=init_vertex_id)?(last_a_id>=0):(last_a_id<0));
    int delta = (vertex1.id()==init_vertex_id)?0:(graph.getEdge(last_vertex, vertex1).get()[last_a_id][a1_id]);

    max_value.addTime(delta);
    if (max_value.getTime() < min_value.getTime()) {
      min_value = max_value;
      min_value.assign(vertex1.id(), a1_id);
    }
  }

  return min_value;
}


FormationGraphDynProgResult formation_graph_dynamic_programming(const FormationGraph& graph, int init_vertex_id) {

  // create the dynamic table
  std::vector<std::vector<std::vector<FormationGraphDynProgResult>>> edge_stored_values;
  for (int v1_id = 0; v1_id < graph.getVertexNum(); v1_id++) {
    auto v1 = graph.getVertexById(v1_id);
    edge_stored_values.emplace_back(graph.getVertexNum());
    for (int v2_id = 0; v2_id < graph.getVertexNum(); v2_id++) {
      auto v2 = graph.getVertexById(v2_id);
      if (graph.hasEdge(v1, v2)) {
        for (unsigned int a1_id = 0; a1_id < v1.get().size(); a1_id++) {
          edge_stored_values[v1_id][v2_id].emplace_back(-1, graph.getVertexNum());
        }
      }
    }
  }

  auto result = formation_graph_dynamic_programming_dfs(graph,
                                                        edge_stored_values,
                                                        graph.getVertexById(init_vertex_id),
                                                        graph.getVertexById(init_vertex_id), // just a dummy
                                                        -1,
                                                        init_vertex_id); // just a dummy

  return result;
}


// **************************************************************************************
//   Formation Graph Search
// **************************************************************************************

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
//    __vv__(free_vertex_ids);
//    __vv__(inverted_free_vertex_ids);
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

//template<typename V, typename E>
//void find_common_prefix_in_tree_dfs(std::unordered_map<int,std::vector<int>>& common_prefix,
//                                    const AdjListGraph<V,E>& graph,
//                                    const typename AdjListGraph<V,E>::Vertex& vertex,
//                                    std::vector<int>& current_path)
//{
//  auto& last_prefix = common_prefix[vertex.id()];
//  if (last_prefix.empty()) {
//    last_prefix = current_path;
//  } else {
//    int i=0;
//    for(int vid : last_prefix) {
//      if (i >= current_path.size() || current_path[i] != vid) break;
//      i++;
//    }
//    last_prefix.clear();
//    for(int j=0; j<i; j++) {
//      last_prefix.push_back(current_path[j]);
//    }
//  }
//
//  for (auto[edge, vertex2]: vertex.outgoingEVList()) {
//    current_path.push_back(vertex2.id());
//    find_common_prefix_in_tree_dfs(common_prefix, graph, vertex2, current_path);
//    current_path.pop_back();
//  }
//}
//
//// this is a brute force method, suitable for a tree only.
//template<typename V, typename E>
//std::unordered_map<int,std::vector<int>> find_common_prefix_in_tree_by_brute_force(const AdjListGraph<V,E>& graph, const typename AdjListGraph<V,E>::Vertex root) {
//  std::unordered_map<int,std::vector<int>> common_prefix(graph.getVertexNum());
//  std::vector<int> current_path;
//  current_path.push_back(root.id());
//  find_common_prefix_in_tree_dfs(common_prefix, graph, root, current_path);
//  return common_prefix;
//}
//
//
//std::vector<int> find_common_prefix_in_formation_graph_by_brute_force(const FormationGraph& graph) {
//  auto common_prefix = find_common_prefix_in_tree_by_brute_force(graph, graph.getVertexById(0));
//
//  std::vector<int> common_prefix_last_vertex(graph.getVertexNum());
//  for(int vid=0; vid<graph.getVertexNum(); vid++) {
//    common_prefix_last_vertex[vid] = common_prefix[vid].back();
//  }
//
//  return common_prefix_last_vertex;
//}


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



FormationGraphResult formation_graph_search_dfs(const FormationGraph& graph,
                                                const std::vector<std::vector<int>>& common_prefix_last_vid,
                                                const FormationGraph::Vertex& current_vertex,
                                                const FormationGraph::Vertex& last_vertex,
                                                std::vector<int>& assignment,
                                                FormationGraphSearchCache& cache,
                                                int init_vertex_id);


FormationGraphResult formation_graph_search_dfs_max(const FormationGraph& graph,
                                                    const std::vector<std::vector<int>>& common_prefix_last_vid,
                                                    const FormationGraph::Vertex& current_vertex,
                                                    const FormationGraph::Vertex& last_vertex,
                                                    std::vector<int>& assignment,
                                                    FormationGraphSearchCache& cache,
                                                    int init_vertex_id)
{
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
                                                   assignment, cache, init_vertex_id);
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
                                                       assignment, cache, init_vertex_id);
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
                                                int init_vertex_id)
{
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
                                                 init_vertex_id);
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
                                                      init_vertex_id);
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


FormationGraphResult formation_graph_search(const FormationGraph& graph, int init_vertex_id) {
  // auto common_prefix_last_vert_bf = find_common_prefix_in_formation_graph_by_brute_force(graph);
  auto common_prefix_last_vertex = find_common_prefix_in_formation_graph(graph, init_vertex_id);
  // assert(common_prefix_last_vert_bf == common_prefix_last_vertex);

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
                                           init_vertex_id);
  return result;
}



// **************************************************************************************
//   testing
// **************************************************************************************

#define TEST_EQ_DOMAIN_SIZE 4


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


FormationGraph makeFormationGraph0(int action_num, std::mt19937& rng) {

  std::uniform_int_distribution<std::mt19937::result_type> rndgen_action_duration(5, TEST_EQ_DOMAIN_SIZE + 5);
  std::uniform_int_distribution<std::mt19937::result_type> rndgen_delta(0, TEST_EQ_DOMAIN_SIZE);

  FormationGraph graph;

  auto v1 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v2 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v3 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v4 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));

  assert(v1.id() == 0);

  graph.addEdge(v1, v2, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v2, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v3, v4, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));

  return graph;
}


FormationGraph makeFormationGraph1(int action_num, std::mt19937& rng) {

  std::uniform_int_distribution<std::mt19937::result_type> rndgen_action_duration(5, TEST_EQ_DOMAIN_SIZE + 5);
  std::uniform_int_distribution<std::mt19937::result_type> rndgen_delta(0, TEST_EQ_DOMAIN_SIZE);

  FormationGraph graph;

  auto v1 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v2 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v3 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v4 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));

  assert(v1.id() == 0);

  graph.addEdge(v1, v2, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v1, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v2, v4, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v3, v4, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));

  return graph;
}

FormationGraph makeFormationGraph1_1(int action_num, std::mt19937& rng) {

  std::uniform_int_distribution<std::mt19937::result_type> rndgen_action_duration(5, TEST_EQ_DOMAIN_SIZE + 5);
  std::uniform_int_distribution<std::mt19937::result_type> rndgen_delta(0, TEST_EQ_DOMAIN_SIZE);

  FormationGraph graph;

  auto v1 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v2 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v3 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v4 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));

  assert(v1.id() == 0);

  graph.addEdge(v1, v2, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v1, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v2, v4, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  // graph.addEdge(v3, v4, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));

  return graph;
}


FormationGraph makeFormationGraph1_2(int action_num, std::mt19937& rng) {

  std::uniform_int_distribution<std::mt19937::result_type> rndgen_action_duration(5, TEST_EQ_DOMAIN_SIZE + 5);
  std::uniform_int_distribution<std::mt19937::result_type> rndgen_delta(0, TEST_EQ_DOMAIN_SIZE);

  FormationGraph graph;

  auto v1 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v2 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v3 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v4 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v5 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));

  assert(v1.id() == 0);

  graph.addEdge(v1, v2, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v1, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v2, v4, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v3, v5, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));

  return graph;
}

FormationGraph makeFormationGraph2(int action_num, std::mt19937& rng) {

  std::uniform_int_distribution<std::mt19937::result_type> rndgen_action_duration(5, TEST_EQ_DOMAIN_SIZE + 5);
  std::uniform_int_distribution<std::mt19937::result_type> rndgen_delta(0, TEST_EQ_DOMAIN_SIZE);

  FormationGraph graph;

  auto v1 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v2 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v3 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v4 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));

  assert(v1.id() == 0);

  graph.addEdge(v1, v2, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v1, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v2, v4, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v3, v4, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v2, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));

  return graph;
}

FormationGraph makeFormationGraph3(int action_num, std::mt19937& rng) {

  std::uniform_int_distribution<std::mt19937::result_type> rndgen_action_duration(20, TEST_EQ_DOMAIN_SIZE + 20);
  std::uniform_int_distribution<std::mt19937::result_type> rndgen_delta(5, TEST_EQ_DOMAIN_SIZE+20);

  FormationGraph graph;

  auto v1 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v2 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v3 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v4 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v5 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v6 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));

  assert(v1.id() == 0);

  graph.addEdge(v1, v2, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v1, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v2, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v2, v4, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v3, v4, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v4, v5, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v5, v6, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));

  return graph;
}


void test_formation_graph_1() {
  auto& rng = Shared::getInstance().getRng();

  int init_vertex_id = 0;

  for(int trial=0; trial<200; trial++) {
    auto graph = makeFormationGraph3(TEST_EQ_DOMAIN_SIZE, rng);
    // std::cout << graph;

    auto start_time_bf = std::chrono::high_resolution_clock::now();
    auto ans_bf = formation_graph_brute_force(graph, init_vertex_id);
    auto end_time_bf = std::chrono::high_resolution_clock::now();

    auto start_time_gg = std::chrono::high_resolution_clock::now();
    auto ans_gg = formation_graph_search(graph, init_vertex_id);
    auto end_time_gg = std::chrono::high_resolution_clock::now();

    auto start_time_dp = std::chrono::high_resolution_clock::now();
    auto ans_dp = formation_graph_dynamic_programming(graph, init_vertex_id);
    auto dp_new_graph = ans_dp.makeReducedGraph(graph);
    auto ans_dp_new_graph = formation_graph_search(dp_new_graph, init_vertex_id);
    auto end_time_dp = std::chrono::high_resolution_clock::now();

    std::cout << ans_bf.getTime() << " " << ans_gg.getTime() << " " << ans_dp.getTime() << " (" << ans_dp.getTime() - ans_gg.getTime() << ")";

    // auto new_ans_dp_time_1 = calc_formation_graph_end_time_schedule_max_time(calc_formation_graph_end_time_schedule(graph, ans_dp.toFormationGraphResultByFirst(graph)));
    // std::cout << " " << new_ans_dp_time_1 << " (" << new_ans_dp_time_1 - ans_gg.getTime() << ")";
    auto new_ans_dp_time_2 = calc_formation_graph_end_time_schedule_max_time(
        calc_formation_graph_end_time_schedule(dp_new_graph, ans_dp_new_graph, init_vertex_id));
    std::cout << " " << new_ans_dp_time_2 << " (" << new_ans_dp_time_2 - ans_gg.getTime() << ")";

    std::cout << " [" << ans_dp_new_graph.getTime() << "]";

    std::cout << "   ";

    std::cout << " " << duration_cast<std::chrono::microseconds>(end_time_bf - start_time_bf).count();
    std::cout << " " << duration_cast<std::chrono::microseconds>(end_time_gg - start_time_gg).count();
    std::cout << " " << duration_cast<std::chrono::microseconds>(end_time_dp - start_time_dp).count();

    if (ans_bf.getTime() != ans_gg.getTime()) {
      std::cout << "  ans_bf != ans_gg (diff = " << ans_bf.getTime() - ans_gg.getTime() << ")" << std::endl;
      break;
    } else if (ans_gg.getTime() < ans_dp.getTime()) {
      std::cout << "  ans_gg < ans_dp" << std::endl;
      break;
    } else if (!ans_gg.hasAllAssigned()) {
      std::cout << "  !ans_gg.hasAllAssigned()" << std::endl;
      __pp__("bf's assignment: ", ans_bf.getAssignment());
      __pp__("gg's assignment: ", ans_gg.getAssignment());
      break;
    } else if (ans_gg.getAssignment() != ans_bf.getAssignment()) {
      auto plan = calc_formation_graph_end_time_schedule(graph, ans_gg, init_vertex_id);
      auto max_time = calc_formation_graph_end_time_schedule_max_time(plan);
      if (max_time != ans_bf.getTime()) {
        std::cout << "  ans_gg.getAssignment() != ans_bf.getAssignment()" << std::endl;
        __pp__("bf's assignment: ", ans_bf.getAssignment());
        __pp__("gg's assignment: ", ans_gg.getAssignment());
        __pp__("And the max time is different: ", max_time);
        break;
      }
    }
    std::cout << std::endl;
  }

  __pp__("okay");
}


// ---------------------------------------------------------------------------------------

FormationGraph makeFormationGraph10(int action_num, std::mt19937& rng) {

  std::uniform_int_distribution<std::mt19937::result_type> rndgen_action_duration(5, TEST_EQ_DOMAIN_SIZE + 5);
  std::uniform_int_distribution<std::mt19937::result_type> rndgen_delta(0, TEST_EQ_DOMAIN_SIZE);

  FormationGraph graph;

  auto v0 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v1 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v2 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));
  auto v3 = graph.addVertex(calcDefaultFormationActions(action_num, rng, rndgen_action_duration));

  assert(v0.id() == 0);

  graph.addEdge(v1, v0, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v1, v2, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v0, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
  graph.addEdge(v2, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));

//  graph.addEdge(v0, v1, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
//  graph.addEdge(v0, v2, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
//  graph.addEdge(v1, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));
//  graph.addEdge(v2, v3, calcDefaultFormationDeltas(action_num, rng, rndgen_delta));

  return graph;
}



void test_formation_graph_2() {
  auto& rng = Shared::getInstance().getRng();

  auto graph = makeFormationGraph10(TEST_EQ_DOMAIN_SIZE, rng);
  int init_vertex_id = 1;

  std::cout << graph << std::endl;

  __line__("formation_graph_brute_force");
  auto ans_bf = formation_graph_brute_force(graph, init_vertex_id);

  __line__("formation_graph_search");
  auto ans_gg = formation_graph_search(graph, init_vertex_id);

  __line__("formation_graph_dynamic_programming");
  auto ans_dp_tmp = formation_graph_dynamic_programming(graph, init_vertex_id);
  auto graph_dp = ans_dp_tmp.makeReducedGraph(graph);
  auto ans_dp = formation_graph_search(graph_dp, init_vertex_id);

  std::cout << ans_bf.getTime() << " " << ans_gg.getTime() << " " << ans_dp.getTime() << " " << "[" << ans_dp_tmp.getTime() << "]" << " (" << ans_dp.getTime() - ans_gg.getTime() << ")";

  std::cout << std::endl;

  __vv__(ans_bf);
  __vv__(ans_gg);
  __vv__(ans_dp);
  __vv__(ans_dp_tmp);

  std::cout << std::endl;

  if (ans_bf.getTime() != ans_gg.getTime()) {
    std::cout << "Error: ans_bf != ans_gg (diff = " << ans_bf.getTime() - ans_gg.getTime() << ")" << std::endl;
  } else if (ans_gg.getTime() < ans_dp_tmp.getTime()) {
    std::cout << "Error: ans_gg < ans_dp_tmp" << std::endl;
  } else if (!ans_gg.hasAllAssigned()) {
    std::cout << "Error: !ans_gg.hasAllAssigned()" << std::endl;
    __pp__("bf's assignment: ", ans_bf.getAssignment());
    __pp__("gg's assignment: ", ans_gg.getAssignment());
  } else if (ans_gg.getAssignment() != ans_bf.getAssignment()) {
    auto plan = calc_formation_graph_end_time_schedule(graph, ans_gg, init_vertex_id);
    auto max_time = calc_formation_graph_end_time_schedule_max_time(plan);
    if (max_time != ans_bf.getTime()) {
      std::cout << "Error: ans_gg.getAssignment() != ans_bf.getAssignment()" << std::endl;
      __pp__("bf's assignment: ", ans_bf.getAssignment());
      __pp__("gg's assignment: ", ans_gg.getAssignment());
      __pp__("And the max time is different: ", max_time);
    }
  }

  std::cout << std::endl;

}


