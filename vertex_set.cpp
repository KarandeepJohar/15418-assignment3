#include "vertex_set.h"

#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <stdio.h>
#include "mic.h"
#include <algorithm>
/**
 * Creates an empty VertexSet with the given type and capacity.
 * numNodes is the total number of nodes in the graph.
 * 
 * Student may interpret type however tey wish.  It may be helpful to
 * use different representations of VertexSets under different
 * conditions, and they different conditions can be indicated by 'type'
 */
VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes)
{
  VertexSet* vs = new VertexSet;
  vs->size = 0;
  vs->capacity = capacity;
  // vs->vertices.reserve(capacity);
  vs->numNodes = numNodes;
  vs->type = type;
  // printf("capacity%d\n",capacity);
  return vs;
}

void freeVertexSet(VertexSet *set)
{

  // delete set->vertices;
  delete set;
}

void addVertex(VertexSet *set, Vertex v)
{ std::vector<Vertex>::iterator it;

  it = std::find(set->vertices.begin(), set->vertices.end(), v);
  if (it != set->vertices.end()) {
    return;
  }
  ENSURES(set->size <= set->capacity);
  set->vertices.push_back(v);
  set->size = set->vertices.size();
}

void removeVertex(VertexSet *set, Vertex v)
{
  ENSURES(set->size>0);
  std::vector<Vertex>::iterator it;

  it = std::find(set->vertices.begin(), set->vertices.end(), v);

  if (it != set->vertices.end()) {
    using std::swap;

    // swap the one to be removed with the last element
    // and remove the item at the end of the container
    // to prevent moving all items after '5' by one
    swap(*it, set->vertices.back());
    set->vertices.pop_back();
  }

  set->size = set->vertices.size();

}

/**
 * Returns the union of sets u and v. Destroys u and v.
 */
VertexSet* vertexUnion(VertexSet *u, VertexSet* v)
{
  // TODO: Implement

  // STUDENTS WILL ONLY NEED TO IMPLEMENT THIS FUNCTION IN PART 3 OF
  // THE ASSIGNMENT
  
  return NULL;
}

