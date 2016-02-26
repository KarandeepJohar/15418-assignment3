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
#include <omp.h>
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
static VertexSet *edgeMapBottomUp(Graph g, VertexSet *u, F &f,
                                  bool removeDuplicates = true)
{
	bool* newDenseVertices = new bool[u->numNodes]();
	int sum = 0;
	bool *ptrDenseVertices = u->denseVertices;
	int uNumNodes = u->numNodes;
	int sum_degrees = 0;
	#pragma omp parallel  reduction(+:sum,sum_degrees)
	{
		#pragma omp  for schedule(dynamic,256)
		for (int i = 0; i < uNumNodes; i++) {
			if (f.cond(i))
			{

				const Vertex* start = incoming_begin(g, i);
				const Vertex* end = incoming_end(g, i);
				for (const Vertex* v = start; v != end; v++) {
					if (ptrDenseVertices[*v] == true  && f.updateNoWorries(*v, i) && !newDenseVertices[i]) {

						sum += 1;
						sum_degrees += outgoing_size(g, i);

						newDenseVertices[i] = true;
					}
					if (!f.cond(i))
					{
						break;
					}

				}
			}

		}
	}
	return newVertexSet(DENSE, sum ,  u->numNodes, newDenseVertices, sum_degrees);
	// }
}
template <class F>
static VertexSet *edgeMapDenseTopDown(Graph g, VertexSet *u, F &f,
                                      bool removeDuplicates = true)
{
	bool* newDenseVertices = new bool[u->numNodes]();
	int sum = 0;
	bool *ptrDenseVertices = u->denseVertices;

	#pragma omp parallel for default(none) shared(g, f, newDenseVertices, ptrDenseVertices)
	for (int i = 0; i < g->num_nodes; i++) {
		Vertex v = i;
		if (ptrDenseVertices[i])
		{
			const Vertex* start = outgoing_begin(g, v);
			const Vertex* end = outgoing_end(g, v);
			for (const Vertex* neigh = start; neigh != end; neigh++) {
				if (f.cond(*neigh) && f.update(v, *neigh))
				{
					newDenseVertices[*neigh] = true;
				}
			}
		}
	}
	int uNumNodes = u->numNodes;
	int sum_degrees = 0;
	#pragma omp parallel for reduction(+:sum,sum_degrees) default(none) shared(uNumNodes, newDenseVertices,g)
	for (int i = 0; i < uNumNodes; ++i)
	{
		if (newDenseVertices[i])
		{
			sum += newDenseVertices[i];
			sum_degrees += outgoing_size(g, i);
		}

	}
	return newVertexSet(DENSE, sum ,  u->numNodes, newDenseVertices, sum_degrees);
	// }
}
template <class F>
static VertexSet *edgeMap(Graph g, VertexSet * u, F & f,
                          bool removeDuplicates = true)
{
	int threshold = num_edges(g) / 100;
	int sum_degrees;
	sum_degrees = u->sum_degrees;
	
	if (sum_degrees  > threshold || u->size > u->numNodes / 100)
	{
		updateDense(u, true);
		return edgeMapBottomUp(g, u, f);
	}
	if (u->size < 8000)
	{
		updateDense(u, true);
		return edgeMapDenseTopDown(g, u, f);
	}
	updateSparse(u, true);
	sum_degrees = 0;
	Vertex *ptrVertices = u->vertices;
	int uSize = u->size;


	if (removeDuplicates)
	{	bool* newDenseVertices = new bool[u->numNodes]();
		int uNumNodes = u->numNodes;
		int sum_degrees = 0;
		int sum = 0;
		#pragma omp parallel default(none)  shared(uSize, ptrVertices, newDenseVertices, f, g, uNumNodes)  reduction(+:sum,sum_degrees)
		{
			#pragma omp for 
			for (int i = 0; i < uSize; ++i)
			{
				Vertex v = ptrVertices[i];
				const Vertex* start = outgoing_begin(g, v);
				const Vertex* end = outgoing_end(g, v);
				for (const Vertex* neigh = start; neigh != end; neigh++) {
					if (f.cond(*neigh) && f.update(v, *neigh) && !newDenseVertices[*neigh])
					{
						sum += 1;
						sum_degrees += outgoing_size(g, *neigh);

						newDenseVertices[*neigh] = true;
					}
				}
			}
		}
		return newVertexSet(DENSE, sum ,  u->numNodes, newDenseVertices, sum_degrees);
	}
	int* offsets = new int[u->size + 1];
	int *degrees;
	degrees = offsets + 1;
	updateSparse(u, true);
	sum_degrees = 0;
	#pragma omp parallel for reduction(+:sum_degrees) default(none) shared(ptrVertices, degrees, g, uSize)
	for (int i = 0; i < uSize; ++i)
	{
		Vertex v = ptrVertices[i];
		degrees[i] = outgoing_size(g, v);
		sum_degrees += degrees[i];
	}
	if (sum_degrees  > threshold || u->size > u->numNodes / 100)
	{

		updateDense(u, true);
		delete[] offsets;
		return edgeMapBottomUp(g, u, f);
	}
	int* finalNeighbours = new int[sum_degrees];
	prefix_sum(offsets, degrees, u->size);
	offsets[0] = 0;
	// #pragma omp parallel for default(none) shared(g, uSize, ptrVertices, f, offsets, finalNeighbours) 
	// for (int i = 0; i < uSize; ++i)
	// {
	// 	Vertex v = ptrVertices[i];
	// 	int offset = offsets[i];
	// 	const Vertex* start = outgoing_begin(g, v);
	// 	const Vertex* end = outgoing_end(g, v);
	// 	int j = 0;
	// 	for (const Vertex* neigh = start; neigh != end; neigh++, j++) {
	// 		if (f.cond(*neigh) && f.update(v, *neigh))
	// 		{
	// 			finalNeighbours[offset + j] = *neigh;
	// 		}
	// 		else {
	// 			finalNeighbours[offset + j] = -1;
	// 		}
	// 	}
	// }

	if (removeDuplicates) {
		remDuplicates(finalNeighbours, sum_degrees, u->numNodes);
	}
	Vertex* newSparseVertices = new Vertex[sum_degrees];
	int new_sum_degrees = 0;
	static bool* tempBoolArray = new bool[threshold]();
	#pragma omp parallel for default(none) shared(sum_degrees, finalNeighbours, tempBoolArray,g) reduction(+:new_sum_degrees)
	for (int i = 0; i < sum_degrees; ++i)
	{
		if (finalNeighbours[i] >= 0)
		{
			new_sum_degrees += outgoing_size(g, i);
			tempBoolArray[i] = true;
		} else {
			tempBoolArray[i] = false;
		}
	}
	int new_sum = packIndices(newSparseVertices,  finalNeighbours, tempBoolArray, sum_degrees);
	VertexSet *trueResult = newVertexSet(SPARSE, sum_degrees, u->numNodes, newSparseVertices, new_sum);
	trueResult->sum_degrees = new_sum_degrees;
	delete[] offsets;
	delete[] finalNeighbours;
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
static VertexSet *vertexMap(VertexSet * u, F & f, bool returnSet = true)
{

	if (returnSet) {
		if (u->type == SPARSE)
		{
			Vertex* newVertices = new Vertex[u->size];
			int sum = 0;
			int uSize = u->size;
			Vertex* ptrVertices = u->vertices;
			int j = -1;
			#pragma omp parallel for default (none) shared(f,uSize, ptrVertices, newVertices,j) reduction(+:sum) schedule(dynamic,256)
			for (int i = 0; i < uSize; ++i)
			{

				if (f(ptrVertices[i])) {
					#pragma omp atomic
					j++;
					newVertices[j] = ptrVertices[i];
					sum += 1;
				}
			}

			VertexSet* result =  newVertexSet(SPARSE, sum, u->numNodes, ptrVertices, sum);
			result->sum_degrees = u->sum_degrees;
			return result;
		} else
		{
			updateDense(u, true);
			bool* newDenseVertices = new bool[u->numNodes]();
			int sum = 0;
			int uNumNodes = u->numNodes;
			bool* ptrDenseVertices = u->denseVertices;
			#pragma omp parallel for default (none) shared(f,uNumNodes, ptrDenseVertices, newDenseVertices) reduction(+:sum)
			for (int i = 0; i < uNumNodes; ++i)
			{
				if (ptrDenseVertices[i])
				{
					if (f(i))
					{

						newDenseVertices[i] = 1;
						sum += 1;
					}
				}
			}

			return newVertexSet(DENSE, sum, u->numNodes, newDenseVertices, u->sum_degrees);
		}
	}
	else {
		if (u->type == SPARSE)
		{
			int uSize = u->size;
			Vertex *ptrVertices = u->vertices;
			# pragma omp parallel for default(none) shared(uSize, f, ptrVertices)
			for (int i = 0; i < uSize; i++)
				f(ptrVertices[i]);
		}
		else {
			int uNumNodes = u->numNodes;
			bool *ptrDenseVertices = u->denseVertices;
			# pragma omp parallel for default(none) shared(uNumNodes, f, ptrDenseVertices) schedule(dynamic,256)
			for (int i = 0; i < uNumNodes; ++i)
			{
				ptrDenseVertices[i]&&f(i);

			}
		}

	}
	return NULL;
}

#endif /* __PARAGRAPH_H__ */
