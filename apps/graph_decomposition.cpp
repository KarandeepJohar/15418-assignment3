#include "paraGraph.h"
#include "vertex_set.h"
#include "graph.h"

/**
	Given a graph, a deltamu per node, the max deltamu value, and the id
	of the node with the max deltamu, decompose the graph into clusters. 
        Returns for each vertex the cluster id that it belongs to inside decomp.
	NOTE: deltamus are given as integers, floating point differences
	are resolved by node id order

**/
class Decompose
{
  public:
    Graph graph;
    int* output;
    bool* thisloop;
    int nodes;

    Decompose(Graph graph, int* output, bool* thisloop, int nodes) :
      graph(graph), output(output), thisloop(thisloop), nodes(nodes) {};

    bool update(Vertex src, Vertex dst) {
        bool flag = false;
        //#pragma omp critical
        if (output[dst]==-1 || (thisloop[dst] && output[src]<output[dst])) {
            __sync_bool_compare_and_swap(&thisloop[dst], false, true);
            __sync_bool_compare_and_swap(&output[dst], -1, output[src]);
            __sync_bool_compare_and_swap(&output[dst], output[dst], output[src]);
            flag=true;
        }

        if (flag) 
          return true;
      return false;
    }

    bool updateNoWorries(Vertex src, Vertex dst) {
         bool flag = false;
        if (output[dst]==-1 || (thisloop[dst] && output[src]<output[dst])) {
            thisloop[dst] = true;
            output[dst] = output[src];
            flag=true;
        }
        if (flag)
          return true;
      return false;
    }

    bool cond(Vertex dst) {
        return true;
    }
};

class vMapFunction
{
  public:
    Graph graph;
    int* output;
    int* DUs;
    int iterat;
    int maxV;

    vMapFunction(Graph graph, int* output, int* DUs, int iterat, int maxV) :
      graph(graph), output(output), DUs(DUs), iterat(iterat), maxV(maxV) {} ;

    bool operator()(Vertex dst) {
      bool flag = false;
    #pragma omp critical
      if (output[dst] == -1 && (iterat > (maxV - DUs[dst]))) {
        output[dst] = dst;
          flag = true;
      }
      if (flag)
          return true;
      return false;
    }
};


void decompose(graph *g, int *decomp, int* dus, int maxVal, int maxId) {
    VertexSet* frontier = newVertexSet(DENSE, 1, g->num_nodes);
    addVertex(frontier, maxId); // vertex with maxDu grows first
    int iter = 0;
    int numNodes = g->num_nodes;
    #pragma omp parallel for default(none) shared(decomp, numNodes)
    for (int i = 0; i< numNodes; i++) {
        decomp[i] = -1;
    }

    bool* newDense = new bool[num_nodes(g)]();
    #pragma omp parallel for default(none) shared(numNodes, newDense)
    for (int i=0; i< numNodes; i++) {
      newDense[i] = true;
    }

    VertexSet* allfrontier = newVertexSet(DENSE, g->num_nodes, g->num_nodes, newDense, 0);
    decomp[maxId] = maxId;
        
    bool* updatedIn = new bool[g->num_nodes]();


    VertexSet* newFrontier;
    VertexSet* newFrontier2;

    while (frontier->size > 0) {
        //newFrontier = newVertexSet(DENSE, 1, numNodes);
        //newFrontierV = newVertexSet(DENSE, 1, numNodes);
        //#pragma omp parallel for default(none) shared(g, frontier, newFrontier,decomp, updatedIn, numNodes)
        //for(int i=0; i < numNodes; i++) {
        //    if (frontier->denseVertices[i]) {
        //        const Vertex* start = outgoing_begin(g, i);
        //        const Vertex* end = outgoing_end(g, i);
        //        for (const Vertex* v=start; v!=end; v++) {
        //            #pragma omp critical
        //            if (decomp[*v]==-1 || updatedIn[*v] ) {
        //                if (decomp[*v] == -1 || decomp[i] < decomp[*v]) {
        //                    updatedIn[*v] = true;
        //                    decomp[*v] = decomp[i];
        //                    addVertex(newFrontier, *v);
        //                }
        //            }
        //        } 
        //    }
        //}
        
        Decompose f(g, decomp, updatedIn, numNodes);
        newFrontier = edgeMap<Decompose>(g, frontier, f);
        iter++;
        
        //#pragma omp parallel for default(none) shared(newFrontier,decomp, maxVal, numNodes, dus, iter)
        //for(int i=0; i<numNodes; i++) {
        //    if (decomp[i] == -1) {
        //        #pragma omp critical
        //        if (iter > (maxVal - dus[i])) {
        //            decomp[i] = i;
        //            addVertex(newFrontier, i);
        //        }
        //    }
        //}
        vMapFunction vmf(g, decomp, dus, iter, maxVal);
        newFrontier2 = vertexMap<vMapFunction>(allfrontier, vmf);
        freeVertexSet(frontier);
        //Merge two frontiers:
        frontier = vertexUnion(newFrontier, newFrontier2);

        #pragma omp parallel for
        for (int i=0; i<numNodes; i++){
          updatedIn[i]=false;
        }
   }

        delete[] updatedIn;
        freeVertexSet(allfrontier);
}
