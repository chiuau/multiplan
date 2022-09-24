#include "graph.h"


void test_graph_1() {
  AdjListGraph<int,int> graph;

  // This example is copied from Page 528 on the Algorithm book

  auto v0 = graph.addVertex(0);
  auto v1 = graph.addVertex(1);
  auto v2 = graph.addVertex(2);
  auto v3 = graph.addVertex(3);
  auto v4 = graph.addVertex(4);

  graph.addEdge(v0, v1, 10);
  graph.addEdge(v0, v2, 5);
  graph.addEdge(v1, v2, 2);
  graph.addEdge(v1, v3, 1);
  graph.addEdge(v2, v1, 3);
  auto e = graph.addEdge(v2, v3, 9);
  graph.addEdge(v2, v4, 2);
  graph.addEdge(v3, v4, 4);
  graph.addEdge(v4, v0, 7);
  graph.addEdge(v4, v3, 6);

  graph.debug();

  DijkstraAlgm dijkstra{graph, v0, v3};

  std::cout << "Shortest distance = " << dijkstra.getCost() << std::endl;
  std::cout << "Shortest path = [";
  for(auto vertex : dijkstra.getPath()) {
    std::cout << " " << vertex.id();
  }
  std::cout << " ]" << std::endl;


  graph.deleteVertex(v1);
  graph.debug();
  graph.deleteEdge(e);
  graph.debug();

}


void test_graph_2() {
  AdjListGraph<int, int> graph;

  // This example is copied from Page 556 on the Algorithm book

  auto v1 = graph.addVertex(1);
  auto v2 = graph.addVertex(2);
  auto v3 = graph.addVertex(3);
  auto v4 = graph.addVertex(4);
  auto v5 = graph.addVertex(5);

  graph.addEdge(v1, v2, 3);
  graph.addEdge(v1, v3, 8);
  graph.addEdge(v1, v5, -4);
  graph.addEdge(v2, v4, 1);
  graph.addEdge(v2, v5, 7);
  graph.addEdge(v3, v2, 4);
  graph.addEdge(v4, v1, 2);
  graph.addEdge(v4, v3, -5);
  graph.addEdge(v5, v4, 6);

  graph.debug();

  FloydWarshallAlgm floyd_warshall{graph};
  floyd_warshall.print();

  std::cout << "Shortest distance = " << floyd_warshall.getCost(v2,v5) << std::endl;
  std::cout << "Shortest path = [";
  for(auto vertex : floyd_warshall.getPath(v2,v5)) {
    std::cout << " " << vertex.id();
  }
  std::cout << " ]" << std::endl;

}
