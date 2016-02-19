#ifndef __PARAGRAPH_H__
#define __PARAGRAPH_H__

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "vertex_set.h"
#include "graph.h"

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
  for (int i=0; i<u->size; i++) {
    const Vertex* start = outgoing_begin(g, u->vertices[i]);
    const Vertex* end = outgoing_end(g, u->vertices[i]);
    for(const Vertex* v=start; v!=end; v++) {

      if (f.cond(*v) && f.update(u->vertices[i], *v))
      {	
        addVertex(trueResult, *v);
    	}
    }
  }
  // printf("size %d\n",trueResult->size );

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
    VertexSet *trueResult = newVertexSet(SPARSE, u->capacity, u->numNodes);
    bool result;
    for (int i=0; i< u->size; i++) {
      result = f(u->vertices[i]);
      if(result) {
        addVertex(trueResult, u->vertices[i]);
      }
    }
    // printf("size %d\n",trueResult->size );
    return trueResult;
  }
  else{
  	for (int i=0; i< u->size; i++) 
  	  f(u->vertices[i]);
  }
  return NULL;
}

#endif /* __PARAGRAPH_H__ */
