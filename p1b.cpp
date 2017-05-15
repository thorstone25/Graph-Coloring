// Code to read graph instances from a file.  Uses the Boost Graph Library (BGL).

#include <iostream>
#include <limits.h>
#include "d_except.h"
#include <fstream>

#include <boost/graph/adjacency_list.hpp>

#define LargeValue 99999999

using namespace std;
using namespace boost;

int const NONE = -1;  // Used to represent a node that does not exist

struct VertexProperties;
struct EdgeProperties;

typedef adjacency_list<vecS, vecS, bidirectionalS, VertexProperties, EdgeProperties> Graph;

struct VertexProperties
{
   pair<int,int> cell; // maze cell (x,y) value
   Graph::vertex_descriptor pred;
   bool visited;
   bool marked;
   int weight;
};

// Create a struct to hold properties for each edge
struct EdgeProperties
{
   int weight;
   bool visited;
   bool marked;
};

void initializeGraph(Graph &g, ifstream &fin)
// Initialize g using data from fin.  
{
   int n, e;
   int j,k;

   fin >> n >> e;
   Graph::vertex_descriptor v;
   
   // Add nodes.
   for (int i = 0; i < n; i++)
      v = add_vertex(g);
   
   for (int i = 0; i < e; i++)
   {
      fin >> j >> k;
      add_edge(j,k,g);  // Assumes vertex list is type vecS
   }
}

void setNodeWeights(Graph &g, int w)
// Set all node weights to w.
{
   pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
   
   for (Graph::vertex_iterator vItr= vItrRange.first; vItr != vItrRange.second; ++vItr)
   {
      g[*vItr].weight = w;
   }
}

int main()
{
   char x;
   ifstream fin;
   string fileName;
   
   // Read the name of the graph from the keyboard or
   // hard code it here for testing.
   
   fileName = "/Users/wmeleis/2560-code/tree2/tree/graph1.txt";
   
   //   cout << "Enter filename" << endl;
   //   cin >> fileName;
   
   fin.open(fileName.c_str());
   if (!fin)
   {
      cerr << "Cannot open " << fileName << endl;
      exit(1);
   }
   
   try
    {
      cout << "Reading graph" << endl;
      Graph g;
      initializeGraph(g,fin);

      cout << "Num nodes: " << num_vertices(g) << endl;
      cout << "Num edges: " << num_edges(g) << endl;
      cout << endl;
      
      // cout << g;
   }
}
