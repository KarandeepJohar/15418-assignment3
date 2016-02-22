#ifndef __PARAGRAPH_H__
#define __PARAGRAPH_H__

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include "vertex_set.h"
#include "graph.h"
#include <set>
#include "mic.h"
#include "graph_internal.h"

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
                          bool removeDuplicates = true)
{
	int* degrees = new int[u->size];
	// printf("edgemap called %d\n",u->size );
	updateSparse(u, true);
	if (u->type == SPARSE)
	{	
		// printf("edgeMap\n");
		int sum_degrees = 0;
		#pragma omp parallel for reduction(+:sum_degrees)
		for (int i = 0; i < u->size; ++i)
		{
			Vertex v = u->vertices[i];
			// printf("%d %d\n",i, v );
			degrees[i] = outgoing_size(g, v);
			sum_degrees += degrees[i];
		}
		int* finalNeighbours = new int[sum_degrees];
		int* offsets = new int[u->size];

		prefix_sum(offsets, degrees, u->size);



		#pragma omp parallel for 
		for (int i = 0; i < u->size; ++i)
		{
			Vertex v = u->vertices[i];
			int offset = offsets[i];
			const Vertex* start = outgoing_begin(g, v);
			const Vertex* end = outgoing_end(g, v);
			// printf("here\n" );
			int j=0;
			for (const Vertex* neigh = start; neigh != end; neigh++, j++) {
				if (f.cond(*neigh) && f.update(v, *neigh))
				{
					finalNeighbours[offset + j] = *neigh;
				}
				else {
					finalNeighbours[offset + j] = -1;
				}
			}
		}

		// printf("99\n");
		if (removeDuplicates) {
			remDuplicates(finalNeighbours, sum_degrees, u->numNodes);
		}
		// printf("85\n");
		Vertex* newSparseVertices = new Vertex[sum_degrees];

		bool* tempBoolArray = new bool[sum_degrees];
		#pragma omp parallel for
		for (int i = 0; i < sum_degrees; ++i)
		{
			if (finalNeighbours[i] >= 0)
			{
				tempBoolArray[i] = true;
			} else {
				tempBoolArray[i] = false;
			}
		}
		// printf("99\n");
		int new_sum = packIndices(newSparseVertices,  finalNeighbours, tempBoolArray, sum_degrees);
		VertexSet *trueResult = newVertexSet(SPARSE, sum_degrees, u->numNodes, newSparseVertices, new_sum);
		delete[] tempBoolArray;
		delete[] degrees;
		delete[] offsets;
		delete[] finalNeighbours;
		// printf("100\n");
		return trueResult;
	}
	else {
		bool* newDenseVertices = new bool[u->numNodes];
		#pragma omp parallel for
		for (int i = 0; i < u->numNodes; ++i)
		{
			newDenseVertices[i] = false;
		}
		#pragma omp parallel for
		for (int i = 0; i < u->numNodes; ++i)
		{
			if (u->denseVertices[i])
			{
				Vertex v = i;
				const Vertex* start = outgoing_begin(g, v);
				const Vertex* end = outgoing_end(g, v);
				int j = 0;
				for (const Vertex* neigh = start; neigh != end; neigh++, j++) {
					if (f.cond(*neigh) && f.update(v, *neigh))
					{
						newDenseVertices[*neigh] = true;
					}
				}

			}
		}
		int sum = 0;
		#pragma omp parallel for reduction(+:sum)
		for (int i = 0; i < u->numNodes; ++i)
		{
			sum += newDenseVertices[i];
		}
		return newVertexSet(DENSE, sum ,  u->numNodes, newDenseVertices);


	}
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
static VertexSet *vertexMap(VertexSet *u, F &f, bool returnSet = true)
{

	if (returnSet) {
		//std::cout << "Size of result" << u->size;
		updateDense(u, true);
		bool* newDenseVertices = new bool[u->numNodes];
		//for (int i = 0; i < u->numNodes; ++i)
		//{
		//	newDenseVertices[i] = false;
		//}

		#pragma omp parallel for
		for (int i = 0; i < u->numNodes; ++i)
		{
			if (u->denseVertices[i])
			{
				newDenseVertices[i] = f(i);
			} else {
				newDenseVertices[i] = false;
      }
		}
		int sum = 0;
		#pragma omp parallel for reduction(+:sum)
		for (int i = 0; i < u->numNodes; ++i)
		{
			sum += newDenseVertices[i];
		}
		return newVertexSet(DENSE, sum, u->numNodes, newDenseVertices);
	}
	else {
		if (u->type == SPARSE)
		{
			# pragma omp parallel for
			for (int i = 0; i < u->size; i++)
				f(u->vertices[i]);
		}
		else {
			# pragma omp parallel for
			for (int i = 0; i < u->numNodes; ++i)
			{
				if (u->denseVertices[i])
				{
					f(i);
				}
			}
		}

	}
	return NULL;
}

#endif /* __PARAGRAPH_H__ */
