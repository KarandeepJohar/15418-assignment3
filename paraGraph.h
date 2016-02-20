#ifndef __PARAGRAPH_H__
#define __PARAGRAPH_H__

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "vertex_set.h"
#include "graph.h"
#include <set>
#include "mic.h"

/*
 * edgeMap --
 * 
 * Students will implement this function.
 * 
 * The input argument f is a class with the following methods defined:
 *   bool update(Vertex src, Vertex dst)
 *   bool cond(Vertex v)
 *
 * See apps/bfs.cpp for an example of such a class definition.
 * 
 * When the argument removeDuplicates is false, the implementation of
 * edgeMap need not remove duplicate vertices from the VertexSet it
 * creates when iterating over edges.  This is a performance
 * optimization when the application knows (and can tell ParaGraph)
 * that f.update() guarantees that duplicate vertices cannot appear in
 * the output vertex set.
 * 
 * Further notes: the implementation of edgeMap is templated on the
 * type of this object, which allows for higher performance code
 * generation as these methods will be inlined.
 */
template <class F>
static VertexSet *edgeMap(Graph g, VertexSet *u, F &f,
    bool removeDuplicates=true)
{
  VertexSet *trueResult = newVertexSet(SPARSE, u->numNodes, u->numNodes);
  //if (removeDuplicates) {
  //    # pragma omp parallel for
  //    for (int i=0; i<u->size; i++) {
  //    Vertex vertex=u->vertices[i];
  //    const Vertex* start = outgoing_begin(g, vertex);
  //    const Vertex* end = outgoing_end(g, vertex);
  //    for(const Vertex* v=start; v!=end; v++) {                
  //      if (f.cond(*v) && f.update(vertex, *v)) { 
  //          #pragma omp critical
  //          addVertex(trueResult, *v);
  //      }
  //    }
  //  }
  //} else {
  #pragma omp parallel
  {
    std::vector<int> vec_private;
    # pragma omp parallel for
    for (int i=0; i<u->size; i++) {
  	  Vertex vertex=u->vertices[i];
      const Vertex* start = outgoing_begin(g, vertex);
      const Vertex* end = outgoing_end(g, vertex);
      for(const Vertex* v=start; v!=end; v++) {
        if (f.cond(*v) && f.update(vertex, *v))
        { 	
          vec_private.push_back(*v);
    	  }
      }
    }
    #pragma omp critical
    trueResult->vertices.insert(trueResult->vertices.end(), vec_private.begin(), vec_private.end());  
  }
  if(removeDuplicates) {
    std::set<int> s(trueResult->vertices.begin(), trueResult->vertices.end());
    trueResult->vertices.assign( s.begin(), s.end());
  }
  trueResult->size = trueResult->vertices.size();
  //}
  return trueResult;
}



/*
 * vertexMap -- 
 * 
 * Students will implement this function.
 *
 * The input argument f is a class with the following methods defined:
 *   bool operator()(Vertex v)
 *
 * See apps/kBFS.cpp for an example implementation.
 * 
 * Note that you'll call the function on a vertex as follows:
 *    Vertex v;
 *    bool result = f(v)
 *
 * If returnSet is false, then the implementation of vertexMap should
 * return NULL (it need not build and create a vertex set)
 */
template <class F>
static VertexSet *vertexMap(VertexSet *u, F &f, bool returnSet=true)
{
  if (returnSet) {
    //std::cout << "Size of result" << u->size;
    VertexSet *trueResult = newVertexSet(SPARSE, u->size, u->numNodes);
    //#pragma omp declare reduction (merge : std::vector<int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
    //std::vector<int> vec;
    #pragma omp parallel
    {
       std::vector<int> vec_private;
       #pragma omp for nowait //fill vec_private in parallel
       for (int i=0; i< u->size; i++) {
          Vertex vertex = u->vertices[i];
          if(f(vertex)) {
            vec_private.push_back(vertex);
          }  
       }
       #pragma omp critical
       trueResult->vertices.insert(trueResult->vertices.end(), vec_private.begin(), vec_private.end());  
    }
    //trueResult->vertices.insert(trueResult->vertices.end(), vec.begin(), vec.end());
    trueResult->size = trueResult->vertices.size();
    return trueResult;
  }
  else{
  	//# pragma omp parallel for
  	for (int i=0; i< u->size; i++) 
  	  f(u->vertices[i]);
  }
  return NULL;
}

#endif /* __PARAGRAPH_H__ */
