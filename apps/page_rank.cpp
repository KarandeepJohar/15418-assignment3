#include "apps/page_rank.h"

#include <stdlib.h>
#include <cmath>
#include <omp.h>

#include "paraGraph.h"

#include "graph.h"
#include "graph_internal.h"

#include <utility>

template <class T>
struct State
{
  State(Graph g, T damping_, T convergence_)
    : graph(g), damping(damping_), convergence(convergence_)
  {
    int numNodes = num_nodes(g);
    pcurr = (T*)(malloc(sizeof(T) * numNodes));
    pnext = (T*)(malloc(sizeof(T) * numNodes));
    diff = (T*)(malloc(sizeof(T) * numNodes));
    psize = (T*)(malloc(sizeof(T) * numNodes));
    #pragma omp parallel for
    for (int i = 0; i < numNodes; i++) {
      pcurr[i] = 1.0 / numNodes;
      pnext[i] = 0.0;
      psize[i] = outgoing_size(graph, i);
    }
  }

  ~State()
  {
    free(diff);
    free(pnext);
    free(pcurr);
    free(psize);
  }

  bool update(Vertex s, Vertex d)
  {
    float add = pcurr[s] / psize[s];
    #pragma omp atomic
    pnext[d] += add;
    return true;
  }
  bool updateNoWorries(Vertex s, Vertex d) {
    pnext[d] += pcurr[s] / psize[s];
    return true;
  }
  bool cond(Vertex v)
  {
    return true;
  }

  T getError()
  {
    int numNodes = num_nodes(graph);

    T error = 0.0;
    #pragma omp parallel for reduction(+:error) schedule(static)
    for (int i = 0; i < numNodes; i++)
      error += diff[i];

    return error;
  }

  Graph graph;
  T* pcurr;
  T* pnext;
  T* diff;
  T* psize;
  T damping;
  T convergence;
};

template <typename T>
struct Local
{
  Local(Graph graph_, T* pcurr_, T* pnext_, T* diff_, T damping_, T addConst_)
    : graph(graph_), pcurr(pcurr_), pnext(pnext_), diff(diff_), damping(damping_), addConst(addConst_)
  {}

  bool operator()(Vertex i)
  {
    pnext[i] = (damping * pnext[i]) + addConst;
    diff[i] = fabs(pnext[i] - pcurr[i]);
    pcurr[i] = 0.f;
    return true;
  }

  Graph graph;
  T* pcurr;
  T* pnext;
  T* diff;

  T damping;
  T convergence;
  T addConst;
};


void pageRank(Graph g, float* solution, float damping, float convergence)
{
  int numNodes = num_nodes(g);
  State<float> s(g, damping, convergence);
  VertexSet* frontier = newVertexSet(SPARSE, numNodes, numNodes);
  // #pragma omp parallel for
  for (int i = 0; i < numNodes; i++) {
    addVertex(frontier, i);
  }
  float addConst = (1 - damping) / numNodes;
  float error = INFINITY;
  while (error > convergence) {
    Local<float> local(g, s.pcurr, s.pnext, s.diff, damping, addConst);

    VertexSet* frontier2 = edgeMap<State<float> >(g, frontier, s);
    // VertexSet* frontier3 = vertexMap<Local<float> >(frontier2, local);
    #pragma omp parallel for
    for (int i = 0; i < numNodes; ++i)
    {
      local.pnext[i] = (damping * local.pnext[i]) + addConst;
      local.diff[i] = fabs(local.pnext[i] - local.pcurr[i]);
      local.pcurr[i] = 0.f;

    }

    freeVertexSet(frontier);
    // freeVertexSet(frontier2);
    frontier = frontier2;

    error = s.getError();
    std::swap(s.pcurr, s.pnext);
  }

  freeVertexSet(frontier);
  memcpy(solution, s.pcurr, sizeof(float) * numNodes);
}
