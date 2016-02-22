#ifndef __VERTEX_SET__
#define __VERTEX_SET__

#include "graph.h"
#include <vector>
typedef enum {
	SPARSE,
	DENSE,
} VertexSetType;

typedef struct {
	int size;     // Number of nodes in the set
	int numNodes; // Number of nodes in the graph
	VertexSetType type;
	// std::vector<Vertex> vertices;
	bool* denseVertices;
	bool denseUpToDate, sparseUpToDate;
	Vertex* vertices;
	int capacity;
} VertexSet;

VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes);
VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes, bool* denseVertices);
VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes, Vertex* vertices);
int packIndices(Vertex* output, Vertex* input, bool* boolArray, int n);

void remDuplicates(Vertex* input, int size, int numNodes);

void freeVertexSet(VertexSet *set);

void addVertex(VertexSet *set, Vertex v);
void removeVertex(VertexSet *set, Vertex v);

VertexSet*  vertexUnion(VertexSet *u, VertexSet* v);
void updateDense(VertexSet *set, bool convert);
void prefix_sum(Vertex* output, int* boolArray, int N);
void updateSparse(VertexSet *set, bool convert );
#endif // __VERTEX_SET__
