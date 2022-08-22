// Assignment 2 20T3 COMP2521: Social Network Analysis
// This program was written by Bejai Cobbing (z5313120)

// CentralityMeasures.c
// This file contains functions relating to CentralityMeasures

#include <stdio.h>
#include <stdlib.h>

#include "CentralityMeasures.h"
#include "Dijkstra.h"
#include "PQ.h"
#include "Graph.h"

static void trace_possible_paths(Vertex curr_vertex, int *numShortestPaths, 
    int *num_not_v_paths, Vertex v, bool goes_through_v, PredNode **pred);

// Find the closeness centrality for each vertex in the given graph and return
// the results in a NodeValues structure
NodeValues closenessCentrality(Graph g) {
	NodeValues nvs;
	nvs.numNodes = GraphNumVertices(g);
	nvs.values = malloc(sizeof(double)*nvs.numNodes);
	
	// Find the closeness centrality for each vertex u in the graph
	for (Vertex u = 0; u < nvs.numNodes; u++) {
	    ShortestPaths pathsRep = dijkstra(g, u);
	    
	    // Iterate through the dist[] array and find the sum of the length of
	    // the shortest paths between u and every other node
	    int n = 0; // The number of nodes u can reach, as defined in spec
	    int sum_shortest_paths = 0; 
	    for (Vertex v = 0; v < nvs.numNodes; v++) {
	        sum_shortest_paths += pathsRep.dist[v]; 
	        if (pathsRep.dist[v] != 0) {
	            n++;
	        }	    
	    }
	    
	    n++; // We need to do this because u can reach itself but it has a value
	         // of 0 in the dist[] array
	    
	    if (n == 1) { // Vertex is isolated; its closeness value should be zero
	        nvs.values[u] = 0;
	        freeShortestPaths(pathsRep);
	        continue;
	    }
	    
	    // Calculate the closeness centrality for the vertex
	    double closenessCentrality = ((double)(n - 1) / (nvs.numNodes - 1))
	        * ((double)(n-1) / sum_shortest_paths);
	    
	    nvs.values[u] = closenessCentrality;
        freeShortestPaths(pathsRep);
	}
	return nvs;
}

// Find the  betweenness centrality for each vertex in the given graph 
// and returns the results in a NodeValues structure
NodeValues betweennessCentrality(Graph g) {
	NodeValues nvs;
	nvs.numNodes = GraphNumVertices(g);
	nvs.values = malloc(sizeof(double)*nvs.numNodes);
	
	// Find the betweenness centrality for each vertex v in the graph
	for (Vertex v = 0; v < nvs.numNodes; v++) {
	    
	    // Start from every vertex s and find the shortest path to every other 
	    // vertex t
	    for (Vertex s = 0; s < nvs.numNodes; s++) {
	        if (s == v) continue;
	        ShortestPaths pathsRep = dijkstra(g, s);
	        
	        // Loop through all of the possible values of t
	        for (Vertex t = 0; t < nvs.numNodes; t++) {
	            if (t == v || t == s) continue;
	            
	            // Create arrays of length of 1 to be passed into recursive 
	            // function, so that the original values can be changed
	            int numShortestPaths[1] = {0};
	            int num_not_v_paths[1] = {0};
	            bool goes_through_v = false; // True if path goes through v
	            
	            // Trace all possible shortest paths from t to s (going 
	            // backwards) by using the pred[] array recursively
	            trace_possible_paths(t, numShortestPaths, num_not_v_paths, v, 
	                goes_through_v, pathsRep.pred);
	           
	            int num_v_paths = numShortestPaths[0] - num_not_v_paths[0];
	            
	            // Calculate the betweenness centrality for these values of s 
	            // and t and add it to the value already stored in the array
	            double BetweennessCentrality = (double)num_v_paths / 
	                numShortestPaths[0];
	            nvs.values[v] += BetweennessCentrality;
	        }
	        freeShortestPaths(pathsRep);
	    }
	}
	return nvs;
}

// Trace all possible shortest paths from t to s (going backwards) by using the 
// pred[] array recursively
static void trace_possible_paths(Vertex curr_vertex, int *numShortestPaths, 
    int *num_not_v_paths, Vertex v, bool goes_through_v, PredNode **pred) {

    PredNode *curr_pred = pred[curr_vertex]; 
        
    // Base Case
    if (curr_pred == NULL) {
        numShortestPaths[0] ++;
        
        if (!goes_through_v) {
            num_not_v_paths[0] ++; 
        }
    }
    
    // Recursively find the lengths of the shortest paths of all of the 
    // predecessors of curr_vertex 
    while (curr_pred != NULL) {
        
        // If curr_vertex is equal to vertex v then set set the flag to True
        if (curr_vertex == v) {
            goes_through_v = true;
        }
        
        trace_possible_paths(curr_pred->v, numShortestPaths, num_not_v_paths,
            v, goes_through_v, pred);
        curr_pred = curr_pred->next;
    }
    return;
}

// Find the  normalised  betweenness centrality for each vertex in the given 
// graph and return the results in a NodeValues structure
NodeValues betweennessCentralityNormalised(Graph g) {
	NodeValues normalised;
	normalised.numNodes = GraphNumVertices(g);
	normalised.values = malloc(sizeof(double)*normalised.numNodes);
	NodeValues un_normalised = betweennessCentrality(g);
	
	// For every non-normalised value, calculate the normalised value and
	// transfer it into the array
	int num_nodes = normalised.numNodes;
	for (Vertex i = 0; i < num_nodes; i++) {
	    double un_normalised_value = un_normalised.values[i];
	    double normalised_value = un_normalised_value * 
	        ((double) 1 / ((num_nodes - 1) * (num_nodes - 2)));
	    normalised.values[i] = normalised_value;
	}
	freeNodeValues(un_normalised);
	return normalised;
}

// Print the values in the given NodeValues structure
void showNodeValues(NodeValues nvs) {
    for (int i = 0; i < nvs.numNodes; i++) {
        printf("%d: %lf\n", i, nvs.values[i]);
    }
    return;
}

// Free all memory associated with the given NodeValues structure.
void freeNodeValues(NodeValues nvs) {
    free(nvs.values);
    return;
}
