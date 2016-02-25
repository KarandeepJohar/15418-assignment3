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
    VertexSet* frontier = newVertexSet(SPARSE, 1, g->num_nodes);
    addVertex(frontier, maxId); // vertex with maxDu grows first
    int iter = 0;
    bool* claimed_by_ball = new bool[g->num_nodes]();
    //bool* visited = new bool[g->num_nodes]();
    //int* claimed_by = new int[g->num_nodes]();

    for (int i = 0; i< g->num_nodes; i++) {
        decomp[i] = -1;
    }

    VertexSet* newFrontier;
    while (frontier->size > 0) {
       //foreach src vertex on frontier {A
        newFrontier = newVertexSet(SPARSE, 1, g->num_nodes);
        for(int i=0; i < frontier->numNodes; i++) {
            if (frontier->denseVertices[i]) {
                //printf("Checking %d\n",i);
                if (!claimed_by_ball[i]) {
                    claimed_by_ball[i] = true;
                    decomp[i] = i;
                }
                const Vertex* start = outgoing_begin(g, i);
                const Vertex* end = outgoing_end(g, i);
          //foreach dest vertex of src {
                for (const Vertex* v=start; v!=end; v++) {
                    //printf("Dest %d of %d\n", *v, i);
                    if (!claimed_by_ball[*v]) {
                        //mark dest as claimed by ball of src vertex
                        if (decomp[*v] == -1){ //|| i < decomp[*v]) {
                            //printf("Claiming %d by %d\n",*v,decomp[i]);
                            claimed_by_ball[*v] = true;
                            decomp[*v] = decomp[i];
                        //add dest to new frontier
                            addVertex(newFrontier, *v);
                        }
                    }
                } 
            }
        }

       iter++;

        // start growing all balls i at the next iter with 
        //     // unvisited center i and with maxDu - dus[i] < iter 
        //foreach vertex i in not_visited {
        for(int i=0; i<g->num_nodes; i++) {
            if (claimed_by_ball[i] == false) {
                if (iter > (maxVal - dus[i])) {
                    //printf("Expanding ball %d\n",i);
                    claimed_by_ball[i] = true;
                    decomp[i] = i;
                    addVertex(newFrontier, i);
                }
            }
        }
        freeVertexSet(frontier);
        frontier = newFrontier;
   }
}
