#include "paraGraph.h"
#include "graph.h"

/**
	Given a graph, a deltamu per node, the max deltamu value, and the id
	of the node with the max deltamu, decompose the graph into clusters. 
        Returns for each vertex the cluster id that it belongs to inside decomp.
	NOTE: deltamus are given as integers, floating point differences
	are resolved by node id order

**/
void decompose(graph *g, int *decomp, int* dus, int maxVal, int maxId) {
    VertexSet* frontier = newVertexSet(DENSE, 1, g->num_nodes);
    addVertex(frontier, maxId); // vertex with maxDu grows first
    int iter = 0;
    bool* claimed_by_ball = new bool[g->num_nodes]();
    int numNodes = g->num_nodes;
    #pragma omp parallel for default(none) shared(decomp, numNodes)
    for (int i = 0; i< numNodes; i++) {
        decomp[i] = -1;
    }
    
    claimed_by_ball[maxId] = true;
    decomp[maxId] = maxId;

    VertexSet* newFrontier;
    while (frontier->size > 0) {
        bool* updatedIn = new bool[g->num_nodes]();
        newFrontier = newVertexSet(SPARSE, 1, numNodes);
        #pragma omp parallel for default(none) shared(g, frontier, newFrontier, claimed_by_ball, decomp, updatedIn, numNodes)
        for(int i=0; i < numNodes; i++) {
            if (frontier->denseVertices[i]) {
                const Vertex* start = outgoing_begin(g, i);
                const Vertex* end = outgoing_end(g, i);
                for (const Vertex* v=start; v!=end; v++) {
                    #pragma omp critical
                    if (!claimed_by_ball[*v] || updatedIn[*v]) {
                        if (decomp[*v] == -1 || decomp[i] < decomp[*v]) {
                            claimed_by_ball[*v] = true;
                            updatedIn[*v] = true;
                            decomp[*v] = decomp[i];
                            addVertex(newFrontier, *v);
                        }
                    }
                } 
            }
        }

        iter++;
        
        #pragma omp parallel for default(none) shared(newFrontier, claimed_by_ball, decomp, maxVal, numNodes, dus, iter)
        for(int i=0; i<numNodes; i++) {
            if (claimed_by_ball[i] == false) {
                #pragma omp critical
                if (iter > (maxVal - dus[i])) {
                    claimed_by_ball[i] = true;
                    decomp[i] = i;
                    addVertex(newFrontier, i);
                }
            }
        }
        delete[] updatedIn;
        freeVertexSet(frontier);
        frontier = newFrontier;
   }

    delete[] claimed_by_ball;
}
