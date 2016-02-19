#include "vertex_set.h"

#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <stdio.h>
#include "mic.h"

/**
 * Creates an empty VertexSet with the given type and capacity.
 * numNodes is the total number of nodes in the graph.
 * 
 * Student may interpret type however they wish.  It may be helpful to
 * use different representations of VertexSets under different
 * conditions, and they different conditions can be indicated by 'type'
 */
VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes)
{
	VertexSet* vs = new VertexSet();
	vs->size = 0;
  vs->capacity = capacity;
  vs->vertices = new Vertex[capacity];
  vs->numNodes = numNodes;
  vs->type = type;
  return vs;
}

void freeVertexSet(VertexSet *set)
{
  delete set->vertices;
  delete set;
}

void addVertex(VertexSet *set, Vertex v)
{
  for (int i = 0; i < set->size; ++i)
  {
    if (set->vertices[i]==v)
    {
      return;
    }
  }
  ENSURES(set->size <= set->capacity);
  set->vertices[set->size++]=v;
}

void removeVertex(VertexSet *set, Vertex v)
{
  ENSURES(set->size>0);
  int i=0;
  for ( i = 0; i < set->size; ++i)
  {
    if (set->vertices[i]==v)
    {
      break;
    }
  }
  for (; i < set->size-1; ++i)
  {
    set->vertices[i]=set->vertices[i+1];
  }
  set->size--;

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

