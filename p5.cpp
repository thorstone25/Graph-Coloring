// Project 5 Local Search: Graph Coloring Steepest Decsent


#include <iostream>
#include <limits.h>
#include "d_except.h"
#include <fstream>
#include <boost/graph/adjacency_list.hpp>
#include <time.h>
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
    int color;
    int degree;
    int ID;
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

void setNodeColors(Graph &g, int c)
// Set all node colors to w.
{
    pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
    
    for (Graph::vertex_iterator vItr= vItrRange.first; vItr != vItrRange.second; ++vItr)
    {
        g[*vItr].color = c;
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
        g[v].degree = 0;
        g[v].ID = i;
    }
    pair<Graph::edge_descriptor, bool> newEdge;
    pair<Graph::edge_descriptor, bool> newEdge2;
    for (int i = 0; i < e; i++)
    {
        fin >> j >> k;
        newEdge = add_edge(j,k,g);  // Assumes vertex list is type vecS
        newEdge2 = add_edge(k,j,g);
        g[target(newEdge.first,g)].degree++; // increase the degree of the target node
        g[source(newEdge.first,g)].degree++; // increase the degree of the source node
    }
    
    setNodesUnvisited(g);
    setNodeColors(g,0);
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

// performs Bubble Sort because I'm lazy and I hate reprogramming sorting algorithms
void sortNodes(vector<Graph::vertex_descriptor> &V, Graph &g)
{
    Graph::vertex_descriptor temp;
    for (int i = 0; i < V.size() - 1; i++)
    {
        for (int j = 1; j < V.size() - i; j++)
        {
            if(g[V[j]].degree > g[V[j-1]].degree)
            {
                temp = V[j];
                V[j] = V[j-1];
                V[j-1] = temp;
            }
        }
    }
    
}

void greedyColoring(Graph &g, int m, Graph &b)
{
    // order all of the nodes in the graph by their degree \\
    // vector of vertex descriptors
    vector<Graph::vertex_descriptor> V;
    
    // resize to number of nodes
    V.resize(num_vertices(g));
    
    // input the nodes into the vector
    pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(g);
    int i = 0;
    for (Graph::vertex_iterator vItr= vItrRange.first; vItr != vItrRange.second; ++vItr)
    {
        V[i] = *vItr;
        i++;
    }
    
    // sort the nodes by degree
    sortNodes(V,g);
    
    // get a vector to store conflicts for each color
    vector<int> conflicts;
    
    // color the nodes in order of highest degree to lowest degree by fewest generated conflicts \\
    // for every node in order of degree ...
    for (int i = 0; i < V.size(); i++)
    {
        // reset to 0 conflicts
        conflicts.clear(); // remove all elements from the vector
        conflicts.resize(m,0); // initialize conflicts to m entries each of value "0"
        
        
        // iterate across the adjacent nodes and get number of conflicts for each color
        pair<Graph::adjacency_iterator, Graph::adjacency_iterator>  vItrRange = adjacent_vertices(V[i], g);
        int j = 0;
        for (Graph::adjacency_iterator vItr= vItrRange.first; vItr != vItrRange.second; ++vItr)
        {
            j++;
            if(g[*vItr].color) // if this nodes has a color (0 if unassigned) and it's color
            {
                conflicts[g[*vItr].color - 1]++; // increment the number of conflicts for this color
            }
        }
        
        // iterate over conflicts and choose color with fewest conflicts
        int best_color = 1;
        long fewest_conflicts = num_vertices(g);
        for(int c = 0; c < m; c++)
        {
            if(conflicts[c] < fewest_conflicts)
            {
                best_color = c + 1; // choose this color
                fewest_conflicts = conflicts[c]; // update best color so far
            }
        }
        
        /*
         // DEBUG: PRINT GRAPH COLORS target node, and conflicts \\
         
         cout << endl << "Target Node: " << g[V[i]].ID << " with " << j << " options checked ----- " << endl;
         for(int j = 0; j< V.size(); j++)
         {
         cout<< "Node " << g[V[j]].ID << ": color " << g[V[j]].color << " ; degree " << g[V[j]].degree << endl;
         }
         
         cout << endl;
         
         for(int c = 0; c < m; c++)
         {
         cout << conflicts[c] << " conflicts with color " << c + 1<< endl;
         }
         
         // DEBUG \\
         */
        // set the node to the color with the fewest conflicts
        g[V[i]].color = best_color;
        
        // #noregrats
    }
    
    // return the best graph
    b = g;
    
    return;
}

Graph getRandomNeighbor(Graph &g)
{
    Graph n = g;
    
    // iterate through all vertices and pick a random one
    pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(n);
    int v = rand() % num_vertices(n);
    Graph::vertex_iterator vItr= vItrRange.first;
    while (vItr != vItrRange.second && v != 0)
    {
        vItr++;
        v--;
    }
    
    // iterate across the adjacent vertices and pick a random one
    pair<Graph::adjacency_iterator, Graph::adjacency_iterator>  aItrRange = adjacent_vertices(*vItr, n);
    int num_adjacent_vertices = 0;
    for (Graph::adjacency_iterator aItr = aItrRange.first; aItr != aItrRange.second; ++aItr) {
        num_adjacent_vertices++;
    }
    int a = rand() % num_adjacent_vertices;
    Graph::adjacency_iterator aItr = aItrRange.first;
    while (aItr != aItrRange.second && a != 0)
    {
        ++aItr;
        a--;
    }
    
    
    //cout << "Switching colors" << endl;
    
    //switch the node and adjacent node's colors
    int temp = n[*aItr].color;
    //cout << "old aItr color " << n[*aItr].color << endl;
    n[*aItr].color = n[*vItr].color;
    //cout << "new aItr color " << n[*aItr].color << endl;
    //cout << "old vItr color " << n[*vItr].color << endl;
    n[*vItr].color = temp;
    //cout << "new vItr color " << n[*vItr].color << endl;
    
    //cout << "old graph conflicts" << getGraphConflicts(g) << endl;
    //cout << "new neighbor conflicts" << getGraphConflicts(n) << endl;
    
    return n;
}

double acceptanceProbability(int currentConflicts, int neighborConflicts, double dk) {
    // If the neighbor is better, accept it with probablity 1
    if (neighborConflicts < currentConflicts) {
        return 1.0;
    }
    // If the neighbor is worse, calculate the acceptance probability
    return exp((currentConflicts - neighborConflicts) / dk);
}

void simulatedAnnealingColoring(Graph &g, int m, clock_t startTime, int t_limit)
{
    cout << "Started simulated annealing" << endl;
    
    Graph c = g;
    double dk = 10000;
    double rate = 0.005;
    
    while (dk > 1)
    {
        //cout << "dk is" << dk << endl;

        // if time limit is exceeded, leave the function
        unsigned long diff = clock()-startTime;
        if( (float)t_limit <  (float)diff / CLOCKS_PER_SEC) {
            cout << "Finished simulated annealing" << endl;
            return;
        }
        
        //Find a random neighbor of g (called n)
        Graph n = getRandomNeighbor(g);
    
        //cout << "current best conflicts" << getGraphConflicts(c) << endl;
        //cout << "current neighbor conflicts" << getGraphConflicts(n) << endl;
        
        //accept the neighbor using probability
        double ap = acceptanceProbability(getGraphConflicts(g), getGraphConflicts(n), dk);
        double r = ((double) rand() / (RAND_MAX));
        //cout << "acceptance prob: " << ap << " and random prob is: " << r << endl;

        if (ap > r)
        {
            //cout << "Accepting neighbor" << endl;
            g = n;
            if (getGraphConflicts(g) < getGraphConflicts(c)) {
                cout << "Updating best graph" << endl;
                c = g;
            }
        }
        dk = dk * (1 - rate);
    }
    
}

Graph getBestNeighbor(Graph &g, clock_t startTime, int t_limit)
{
    Graph n = g;
    Graph b = g;
    
    // iterate through all vertices
    pair<Graph::vertex_iterator, Graph::vertex_iterator> vItrRange = vertices(n);
    for (Graph::vertex_iterator vItr= vItrRange.first; vItr != vItrRange.second; ++vItr)
    {

        // iterate across the adjacent vertices
        pair<Graph::adjacency_iterator, Graph::adjacency_iterator>  aItrRange = adjacent_vertices(*vItr, n);
        for (Graph::adjacency_iterator aItr = aItrRange.first; aItr != aItrRange.second; ++aItr)
        {
            
            // if time limit is exceeded, leave the function
            unsigned long diff = clock()-startTime;
            if( (float)t_limit <  (float)diff / CLOCKS_PER_SEC) { return n; }
            
            //cout << "Switching colors" << endl;
            
            //switch the node and adjacent node's colors
            int temp = n[*aItr].color;
            //cout << "old aItr color " << n[*aItr].color << endl;
            n[*aItr].color = n[*vItr].color;
            //cout << "new aItr color " << n[*aItr].color << endl;
            //cout << "old vItr color " << n[*vItr].color << endl;
            n[*vItr].color = temp;
            //cout << "new vItr color " << n[*vItr].color << endl;

            //cout << "current best conflicts" << getGraphConflicts(b) << endl;
            //cout << "current neighbor conflicts" << getGraphConflicts(n) << endl;

            //if this coloring has fewer conflics, set b to be this coloring
            if (getGraphConflicts(n) < getGraphConflicts(b)) {
                cout << "Updating best neighbor" << endl;
                b = n;
            } else {
                //switch back the colors
                int temp = n[*aItr].color;
                n[*aItr].color = n[*vItr].color;
                n[*vItr].color = temp;
            }
        }
    }
    return b;
}

void steepestDescentColoring(Graph &g, int m, clock_t startTime, int t_limit)
{
    cout << "Started steepest descent" << endl;

    bool done = false;
    
    while (!done)
    {
        // if time limit is exceeded, leave the function
        unsigned long diff = clock()-startTime;
        if( (float)t_limit <  (float)diff / CLOCKS_PER_SEC) { return; }
        
        //Find the best neighbor of g (called n)
        Graph n = getBestNeighbor(g, startTime, t_limit);
        
        //if the conflicts in n are less than conflicts in g, set g to n
        //else, we are done and g is the best
        
        if (getGraphConflicts(n) < getGraphConflicts(g))
        {
            cout << "Updating best graph" << endl;
            g = n;
        } else {
            done = true;
        }
    }
    cout << "Finished steepest descent" << endl;

    return;
    
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


int p2b(string filenameext)
{
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
        setNodeColors(g,0);
        b = g;
        
        cout << "Num colors: " << m << endl;
        cout << "Num nodes: " << num_vertices(g) << endl;
        cout << "Num edges: " << num_edges(g) << endl;
        
        greedyColoring(g,m,b);
        unsigned long diff = clock()-startTime;
        float runTime = (float) diff / CLOCKS_PER_SEC;
        cout << "Greedy Algorithm Total Runtime: " << runTime << "s" << endl << endl;
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

int p5(string filenameext) {
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
        //setNodeColors(g,0);
        b = g;
        setNodeColors(b,0); // set all nodes in b to the same color for maximum conflict
        
        cout << "Num colors: " << m << endl;
        cout << "Num nodes: " << num_vertices(g) << endl;
        cout << "Num edges: " << num_edges(g) << endl;
        cout << endl;

        //exhaustiveColoring(g,m,b,startTime,30);
        greedyColoring(g, m, b);
        cout << "Found initial solution" << endl;
        cout << "Initial conflicts: " << getGraphConflicts(b) << endl;
        //steepestDescentColoring(b, m, startTime, 300);
        simulatedAnnealingColoring(b, m, startTime, 300);
        cout << "Conflicts after local search: " << getGraphConflicts(b) << endl;

        unsigned long diff = clock()-startTime;
        float runTime = (float) diff / CLOCKS_PER_SEC;
        //cout << "Steepest Descent Algorithm Total Runtime: " << runTime << "s" << endl << endl;
        cout << "Simulated Annealing Algorithm Total Runtime: " << runTime << "s" << endl << endl;
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
    p5("color12-3.input");
    p5("color12-4.input");
    p5("color24-4.input");
    p5("color24-5.input");
    p5("color48-5.input");
    p5("color48-6.input");
    p5("color96-6.input");
    p5("color96-7.input");
    p5("color192-6.input");
    p5("color192-7.input");
    p5("color192-8.input");
}
