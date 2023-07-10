#ifndef SIM_FORMATION_GRAPH_H
#define SIM_FORMATION_GRAPH_H

#include <iostream>
#include <iomanip>
#include <random>
#include <unordered_map>
#include <limits>

#include "shared.h"
#include "graph.h"

#define EXPR_TIME_STEP 0.025
#define EXPR_TIME_LIMIT 60.0


using FormationTime = int;
using FormationActions = std::vector<FormationTime>;
using FormationDeltas = std::vector<std::vector<FormationTime>>;
using FormationGraph = AdjListGraph<FormationActions,FormationDeltas>;
using FormationSchedule = std::vector<FormationTime>;


/**************************************************************************************
 *  FormationGraphResult
 *
 *  This class holds the result of the brute-force search and the graph search
 ************************************************************************************** */

class FormationGraphResult {

  FormationTime time;
  std::vector<int> assignment;

public:

  FormationGraphResult() : time(-1.0) {}

  FormationGraphResult(FormationTime time, int vertex_num) :
    time(time), assignment(vertex_num, -1)
  {
    // do nothing
  }

  FormationGraphResult(FormationTime time, const std::vector<int>& assignment) :
    time(time), assignment(assignment)
  {
    // do nothing
  }


  FormationTime getTime() const { return time; }

  void setTime(FormationTime t) { time = t; }

  void addTime(FormationTime t) { time += t; }

  bool isValid() const { return time >= 0.0; }

  const std::vector<int>& getAssignment() const { return assignment; }

  int getActionId(int plan_id) const { return assignment[plan_id]; }


  bool isAssigned(int vid) const { return assignment[vid] >= 0; }

  bool hasAllAssigned() const {
    for(auto aid : assignment) {
      if (aid < 0) return false;
    }
    return true;
  }

  void assign(int vid, int aid) {
    assert(!isAssigned(vid));        // *** why it failed? ***
    assignment[vid] = aid;
  }

  void assign(const std::vector<int>& more_assignment) {
    assert(assignment.size() == more_assignment.size());
    for(unsigned int vid=0; vid<assignment.size(); vid++) {
      if (more_assignment[vid] >= 0) {
        assert(!isAssigned(vid));       // *** why it failed? ***
        assignment[vid] = more_assignment[vid];
      }
    }
  }

  friend std::ostream& operator<<(std::ostream& out, const FormationGraphResult& graph);

};



/**************************************************************************************
 *  FormationGraphDynProgResult
 *
 *  This class holds the result of the dynamic programming algorithm
 ************************************************************************************** */

class FormationGraphDynProgResult {

  FormationTime time;
  std::vector<std::vector<int>> assignment;

public:

  FormationGraphDynProgResult(FormationTime time, int vertex_num) :
      time(time), assignment(vertex_num)
  {
    // do nothing
  }

  FormationGraphDynProgResult(FormationTime time, const std::vector<std::vector<int>>& assignment) :
      time(time), assignment(assignment)
  {
    // do nothing
  }


  FormationTime getTime() const { return time; }

  void setTime(FormationTime t) { time = t; }

  void addTime(FormationTime t) { time += t; }

  bool isValid() const { return time >= 0.0; }

  const std::vector<std::vector<int>>& getAssignment() const { return assignment; }

  std::vector<int> getActionIds(int plan_id) const { return assignment[plan_id]; }


  bool isAssigned(int vid) const { return !assignment[vid].empty(); }

  bool hasAllAssigned() const {
    for(auto aids : assignment) {
      if (aids.empty()) return false;
    }
    return true;
  }

  void assign(int vid, int aid) {
    if (std::find(assignment[vid].begin(), assignment[vid].end(), aid) == assignment[vid].end()) {
      assignment[vid].push_back(aid);
    }
  }

  void assign(int vid, const std::vector<int>& aids) {
    for(auto aid : aids) {
      assign(vid, aid);
    }
  }

  void assign(const std::vector<std::vector<int>>& more_assignment) {
    assert(assignment.size() == more_assignment.size());
    for(unsigned int vid=0; vid<assignment.size(); vid++) {
      assign(vid, more_assignment[vid]);
    }
  }

  FormationGraphResult toFormationGraphResultByFirst(const FormationGraph& graph) const;

  FormationGraph makeReducedGraph(const FormationGraph& old_graph) const;

  friend std::ostream& operator<<(std::ostream& out, const FormationGraphDynProgResult& graph);

};


/**************************************************************************************
 *  Helper functions
 ************************************************************************************** */

std::ostream& operator<<(std::ostream& out, const FormationGraph& graph);

std::ostream& operator<<(std::ostream& out, const FormationGraphResult& graph_result);

std::ostream& operator<<(std::ostream& out, const FormationGraphDynProgResult& graph_result);


FormationSchedule calc_formation_graph_end_time_schedule(const FormationGraph& graph, const FormationGraphResult& result, int init_vertex_id);

FormationTime calc_formation_graph_end_time_schedule_max_time(const FormationSchedule& end_time_schedule);

FormationSchedule calc_formation_graph_schedule(const FormationGraph& graph, const FormationGraphResult& result, const FormationSchedule& end_time_schedule);

void advanceFormationGraphSchedule(FormationSchedule& schedule, FormationTime time);


/**************************************************************************************
 *  The main functions of the three algorithms
 ************************************************************************************** */

FormationGraphResult formation_graph_brute_force(const FormationGraph& graph, int init_vertex_id);

FormationGraphDynProgResult formation_graph_dynamic_programming(const FormationGraph& graph, int init_vertex_id);

FormationGraphResult formation_graph_search(const FormationGraph& graph, int init_vertex_id, double time_limit = EXPR_TIME_LIMIT);


/**************************************************************************************
 *  Tests
 ************************************************************************************** */

void pure_formation_graph_expr_1();


#endif //SIM_FORMATION_GRAPH_H
