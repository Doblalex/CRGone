// A C++ program to find biconnected components in a given undirected graph
#pragma once

#include "util.hpp"

#define NIL -1

// A class that represents an directed graph
class Graph {
    int V; // No. of vertices
    list<int>* adj; // A dynamic array of adjacency lists
 
    // A Recursive DFS based function used by SCC()
    void SCCUtil(int u, int disc[], int low[],
                 stack<int>* st, bool stackMember[], vector<vector<int>>& sccs);
 
public:
    Graph(int V); // Constructor
    void addEdge(int v,
                 int w); // function to add an edge to graph
    vector<vector<int>> SCC(); // prints strongly connected components
};
