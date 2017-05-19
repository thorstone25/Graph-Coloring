// Code to read graph instances from a file.  Uses the Boost Graph Library (BGL).

#include <iostream>
#include <limits.h>
#include "d_except.h"
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <time.h>
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
	int color;
};

// Create a struct to hold properties for each edge
struct EdgeProperties
{
	int weight;
	bool visited;
	bool marked;
};

void setNodesUnvisited(Graph &g)
// Set all node colors to w.
{
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);

	for (Graph::vertex_iterator vItr= vItrRange.first; vItr != vItrRange.second; ++vItr)
	{
		g[*vItr].visited = false;
	}
}

void initializeGraph(Graph &g, ifstream &fin)
// Initialize g using data from fin.  
{
	int n, e;
	int j,k;

	fin >> n >> e;
	Graph::vertex_descriptor v;

	// Add nodes.
	for (int i = 0; i < n; i++)
	{
		v = add_vertex(g);
	}

	for (int i = 0; i < e; i++)
	{
		fin >> j >> k;
		add_edge(j,k,g);  // Assumes vertex list is type vecS
	}
	
	setNodesUnvisited(g);
}

void setNodeColors(Graph &g, int c)
// Set all node colors to w.
{
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);

	for (Graph::vertex_iterator vItr= vItrRange.first; vItr != vItrRange.second; ++vItr)
	{
		g[*vItr].color = c;
	}
}

int getGraphConflicts(Graph &g)
{
	/* loop through all the edges, +1 for every edge with equal node values*/
	// conflict placeholder
	int conflicts = 0;
	
	// edge iteration range
	pair<Graph::edge_iterator, Graph::edge_iterator> eItrRange = edges(g);

	// loop though all edges
	for (Graph::edge_iterator eItr= eItrRange.first; eItr != eItrRange.second; ++eItr)
	{
		// get the target
		Graph::vertex_descriptor v = target(*eItr, g);
		
		// get the source
		Graph::vertex_descriptor u = source(*eItr, g);
		
		// if the colors are equal, +1 conflict!
		if(g[v].color == g[u].color)
		{
			conflicts++;
		}
	}
	return conflicts;
}

void exhaustiveColoring(Graph &g, int m, Graph &b, clock_t startTime, int t_limit)
{

	// if time limit is exceeded, leave the function
    unsigned long diff = clock()-startTime;
    if( (float)t_limit <  (float)diff / CLOCKS_PER_SEC) 
	{
        return;
    }
    	
	/* if all vertices have been visited, compute the conflicts, and compare to the best graph so far
	   if a vertex has yet to be visited, loop over all the possible color combinations*/
		
	// get iterator to the vertices
	pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
	
	// get reference for the next unvisited node
	Graph::vertex_iterator vItr;
	
	// find the next unvisited node
	for (vItr= vItrRange.first; vItr != vItrRange.second; ++vItr)
	{
		// if this node has not been visited, exit loop
		if( (g[*vItr].visited ) == false)
		{
			break;
		}
	}
	
	// if at the end
	if(vItr == vItrRange.second)
	{
		// get the number of conflicts of this graph
		int g_conflicts = getGraphConflicts(g);
		
		// get the number of conflicts of the best graph so far (b)
		int b_conflicts = getGraphConflicts(b);
		
		// if this graph has fewer conflicts, set b to be this graph
		if(g_conflicts < b_conflicts)
		{
			b = g;
		}
	}

	/* if not the end, loop through the different colors for this node, and do all permutations for the rest of the nodes */
	else
	{
		// set the node as visited
		g[*vItr].visited = true;
		
		// loop through different colors
		for( int i = 0; i < m; i++)
		{
			// set this node to color i
			g[*vItr].color = i + 1;
			
			// set color of the remaining nodes
			exhaustiveColoring(g, m, b, startTime, t_limit);
		}
		
		// unvisit the node when finished going through all the colors
		g[*vItr].visited =  false;
	}
	
}

void generateOutput (Graph &g, string filename, float runtime) 
{
    ofstream myfile;
    string filenameext = filename + ".output";
    myfile.open (filenameext.c_str());
    
    myfile << "Total conflicts: " << getGraphConflicts(g) << endl;
    for(int counter = 0; counter < num_vertices(g); counter++){
        myfile << counter << " : " << g[counter].color << endl;
    }
    myfile << "Time to completion: " << runtime << " seconds." << endl;
    
    myfile.close();
}


int p1b(string filenameext)
{
	char x;
	ifstream fin;
	string filename;

	// Read the name of the graph from the keyboard or
	// hard code it here for testing.

	filename = filenameext.substr(0, filenameext.find_last_of("."));
	fin.open(filenameext.c_str());
	
	//   cout << "Enter filename" << endl;
	//   cin >> fileName;

	if (!fin)
	{
		cerr << "Cannot open " << filenameext << endl;
		exit(1);
	}

	try
	{
		int m; // number of colors
		 // get time reference for start
	    clock_t startTime;
	    startTime = clock();
	    
		cout << "Reading graph" << endl;
		fin >> m;
		Graph g, b;
		initializeGraph(g,fin);
		b = g;
		setNodeColors(b,1); // set all nodes in b to the same color for maximum conflict

		cout << "Num colors: " << m << endl;
		cout << "Num nodes: " << num_vertices(g) << endl;
		cout << "Num edges: " << num_edges(g) << endl;
		cout << endl;
		
		exhaustiveColoring(g,m,b,startTime,600);
		unsigned long diff = clock()-startTime;
    	float runTime = (float) diff / CLOCKS_PER_SEC;
		cout << endl << "Exhaustive Algorithm Total Runtime: " << runTime << "s" << endl;
		generateOutput(b, filename, runTime);
		// cout << g;
		// exit(0);
		return 0;
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

int main()
{
	p1b("color12-3.input");
	p1b("color12-4.input");
	p1b("color24-4.input");
	p1b("color24-5.input");
	p1b("color48-5.input");
	p1b("color48-6.input");
	p1b("color96-6.input");
	p1b("color96-7.input");
	p1b("color192-6.input");
	p1b("color192-7.input");
	p1b("color192-8.input");
}
