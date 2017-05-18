// Code to read graph instances from a file.  Uses the Boost Graph Library (BGL).

#include <iostream>
#include <limits.h>
#include "d_except.h"
#include <fstream>
#include <boost/graph/adjacency_list.hpp>

#define LargeValue 99999999

using namespace std;
using namespace boost;

//int const NONE = -1;  // Used to represent a node that does not exist

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

int getGraphConflicts(Graph &g)
{
	/* loop through all the edges, +1 for every edge with equal node values*/
	// conflict placeholder
	int conflicts = 0;
	
	// edge iteration range
	pair<Graph::edge_iterator, Graph::edge_iterator> eItrRange = edges(g);
	
	// vertex 
	
	// loop though all edges
	for (Graph::edge_iterator eItr= eItrRange.first; eItr != eItrRange.second; ++eItr)
	{
		// get the target
		Graph::vertex_descriptor v = target(*eItr, g);
		
		// get the source
		Graph::vertex_descriptor u = source(*eItr, g);
		
		// if the weights are equal, +1 conflict!
		if(g[v].weight == g[u].weight)
		{
			conflicts++;
		}
	}
	return conflicts;
}

void findBestColoring(Graph &g, int m, Graph::vertex_iterator &current, Graph::vertex_iterator &end, Graph &b, clock_t startTime, int t_limit)
{
    // if time limit is exceeded, leave the function
    unsigned long diff = clock()-startTime;
    if( (float)t_limit <  (float)diff / CLOCKS_PER_SEC) {
        cout << "Time limit exceeded." << endl;
        return;
    }
    
    /* if this is the end of the vertices, compute the conflicts, and compare to the best graph so far*/
	// if at the end of the list of vertices
	if(current == end)
	{
		// get the number of conflicts of this graph
		int g_conflicts = getGraphConflicts(g);
		
		// get the number of conflicts of the best graph so far (b)
		int b_conflicts = getGraphConflicts(b);
		
		// if this graph has fewer conflicts, set b to be this graph
		b = (g_conflicts < b_conflicts ? g : b);
	}

	/* if not the end, loop through the different colors for this node, and do all permutations for the rest of the nodes */
	else
	{
		for( int i = 0; i < m; i++)
		{
			// set this node to color i
			g[*current].weight = i;
			
			// set color of the remaining nodes
			++current;
			findBestColoring(g,m,current,end,b, startTime, t_limit);
		}
	}
	
}

/* loop through all permutations of coloring to find the one with the least conflicts*/
void exhaustiveColoring(Graph &g, int m, int t)
{
    // get time reference for start
    clock_t startTime;//, endTime;
    startTime = clock();
    
    // get reference to a new graph
	Graph b = g;
	
	// get iterator to the vertices
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
	
	// iterate through all vertices and check all possible colorings
	findBestColoring(g,m,vItrRange.first,vItrRange.second,b, startTime, t);
	
	// set g to be the best graph
	g = b;
    
    unsigned long diff = clock()-startTime;
    cout << endl << "Exhaustive Algorithm Total Runtime: " << (float) diff / CLOCKS_PER_SEC << "s" << endl;
}

void generateOutput (Graph &g, string filename) {
    ofstream myfile;
    string filenameext = filename + ".output";
    myfile.open (filenameext.c_str());
    
    myfile << "Total conflicts: " << getGraphConflicts(g) << endl;
    for(int counter = 0; counter < num_vertices(g); counter++){
        myfile << counter << " : " << g[counter].weight << endl;
    }
    
    myfile.close();
}


int main()
{
	//char x;
	ifstream fin;
	string filenameext, filename;

	// Read the name of the graph from the keyboard or
	// hard code it here for testing.

	filenameext = "color12-3.input";

	//   cout << "Enter filename" << endl;
	//   cin >> fileName;

    filename = filenameext.substr(0, filenameext.find_last_of("."));
    
	fin.open(filenameext.c_str());
	if (!fin)
	{
		cerr << "Cannot open " << filenameext << endl;
		exit(1);
	}

	try
	{
		int m; // number of colors
		cout << "Reading graph" << endl;
		fin >> m;
		Graph g;
		initializeGraph(g,fin);

		cout << "Num colors: " << m << endl;
		cout << "Num nodes: " << num_vertices(g) << endl;
		cout << "Num edges: " << num_edges(g) << endl;
		cout << endl;
		
		exhaustiveColoring(g,m,600);
        

		// cout << g;
        generateOutput(g, filename);
        
		exit(0);
	}
	catch (indexRangeError &ex) 
	{ 
		cout << ex.what() << endl; exit(1);
	}
	catch (rangeError &ex)
	{
		cout << ex.what() << endl; exit(1);
	}
}
