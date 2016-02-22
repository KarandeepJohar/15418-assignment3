#include "vertex_set.h"

#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <stdio.h>
#include <omp.h>
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
bool ispowerof2( int x) {
     return x && !(x & (x - 1));
}

unsigned long upper_power_of_two(unsigned long v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}


void prefix_sum(Vertex* output, bool* boolArray, int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
      if (boolArray[i]) {
        output[i] = 1;
      } else {
        output[i] = 0;
      }
    }
    if (ispowerof2(N)) {
        // upsweep phase.
        for (int twod = 1; twod < N; twod *= 2) {
            int twod1 = twod * 2;
            #pragma omp parallel for
            for (int i = 0; i < N; i += twod1) {
                output[i + twod1 - 1] += output[i + twod - 1];
            }
        }
        output[N - 1] = 0;
        // for (int i = 0; i<N; i++) {
        //     printf("Bool: %d ",boolArray[i]);
        //     printf("Upsweep: %d\n",output[i]);
        // }
        // downsweep phase.
        for (int twod = N / 2; twod >= 1; twod /= 2) {
            int twod1 = twod * 2;
            #pragma omp parallel for
            for (int i = 0; i < N; i += twod1) {
                int t = output[i + twod - 1];
                output[i + twod - 1] = output[i + twod1 - 1];
                output[i + twod1 - 1] += t;
            }
        }
        // for (int i = 0; i<N; i++) {
        //     printf("Bool: %d ",boolArray[i]);
        //     printf("Sum: %d\n",output[i]);
        // }
    } else {
        // printf("Going up\n");
        int powN = upper_power_of_two(N);
        //Vertex* opN = (Vertex *)calloc(powN*sizeof(int),0);
        Vertex* opN = new Vertex[powN]();
        memcpy(opN, output, N*sizeof(int));
        // upsweep phase.
        for (int twod = 1; twod < powN; twod *= 2) {
            int twod1 = twod * 2;
            #pragma omp parallel for
            for (int i = 0; i < powN; i += twod1) {
                opN[i + twod1 - 1] += opN[i + twod - 1];
            }
        }
        opN[powN - 1] = 0;
        // for (int i = 0; i<N; i++) {
        //     printf("Bool: %d ",boolArray[i]);
        //     printf("Upsweep: %d\n",opN[i]);
        // }
        // downsweep phase.
        for (int twod = powN / 2; twod >= 1; twod /= 2) {
            int twod1 = twod * 2;
            #pragma omp parallel for
            for (int i = 0; i < powN; i += twod1) {
                int t = opN[i + twod - 1];
                opN[i + twod - 1] = opN[i + twod1 - 1];
                opN[i + twod1 - 1] += t;
            }
        }
        // for (int i = 0; i<N; i++) {
        //     printf("Bool: %d ",boolArray[i]);
        //     printf("Sum: %d\n",opN[i]);
        // }

        memcpy(output, opN, N*sizeof(int));
        //free(opN);
        delete[] opN;
    }

}

void prefix_sum(Vertex* output, int* boolArray, int N) {
    memcpy(output, boolArray, N*sizeof(int));
    if (ispowerof2(N)) {
        // upsweep phase.
        for (int twod = 1; twod < N; twod *= 2) {
            int twod1 = twod * 2;
            //#pragma omp parallel for
            for (int i = 0; i < N; i += twod1) {
                output[i + twod1 - 1] += output[i + twod - 1];
            }
        }
        output[N - 1] = 0;
        // for (int i = 0; i<N; i++) {
        //     printf("Bool: %d ",boolArray[i]);
        //     printf("Upsweep: %d\n",output[i]);
        // }
        // downsweep phase.
        for (int twod = N / 2; twod >= 1; twod /= 2) {
            int twod1 = twod * 2;
            //#pragma omp parallel for
            for (int i = 0; i < N; i += twod1) {
                int t = output[i + twod - 1];
                output[i + twod - 1] = output[i + twod1 - 1];
                output[i + twod1 - 1] += t;
            }
        }
        // for (int i = 0; i<N; i++) {
        //     printf("Bool: %d ",boolArray[i]);
        //     printf("Sum: %d\n",output[i]);
        // }
    } else {
        // printf("Going up\n");
        int powN = upper_power_of_two(N);
        //Vertex* opN = (Vertex *)calloc(powN*sizeof(int),0);
        Vertex* opN = new Vertex[powN]();
        memcpy(opN, output, N*sizeof(int));
        // upsweep phase.
        for (int twod = 1; twod < powN; twod *= 2) {
            int twod1 = twod * 2;
            //#pragma omp parallel for
            for (int i = 0; i < powN; i += twod1) {
                opN[i + twod1 - 1] += opN[i + twod - 1];
            }
        }
        opN[powN - 1] = 0;
        // for (int i = 0; i<N; i++) {
        //     printf("Bool: %d ",boolArray[i]);
        //     printf("Upsweep: %d\n",opN[i]);
        // }
        // downsweep phase.
        for (int twod = powN / 2; twod >= 1; twod /= 2) {
            int twod1 = twod * 2;
            //#pragma omp parallel for
            for (int i = 0; i < powN; i += twod1) {
                int t = opN[i + twod - 1];
                opN[i + twod - 1] = opN[i + twod1 - 1];
                opN[i + twod1 - 1] += t;
            }
        }
        // for (int i = 0; i<N; i++) {
        //     printf("Bool: %d ",boolArray[i]);
        //     printf("Sum: %d\n",opN[i]);
        // }

        memcpy(output, opN, N*sizeof(int));
        //free(opN);
        delete[] opN;
    }

}

int packIndices(Vertex* output, Vertex* input, bool* boolArray, int n) {
    Vertex* sums = new Vertex[n];
    prefix_sum(sums, boolArray, n);
    int sum=0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < n; ++i)
    {
        if (boolArray[i])
        {
            output[sums[i]] = input[i];
            sum++;
        }

    }
    delete[] sums;
    return sum;
}
void remDuplicates(Vertex* input, int size, int numNodes) {
    

    Vertex* flags = new Vertex[numNodes];
    #pragma omp parallel for
    for (int i = 0; i < numNodes; ++i)
    {
        flags[i] = -1;
    }
    #pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        if (input[i] != -1 && flags[input[i]] == -1)
        {
            __sync_bool_compare_and_swap(&flags[input[i]], -1, i);
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < size; ++i)
    {
        if (input[i] != -1)
        {
            if (flags[input[i]] == i)
            {
                flags[input[i]] = -1;
            } else {
                input[i] = -1;
            }
        }
    }
    delete[] flags;
}
VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes)
{
    VertexSet* vs = new VertexSet;
    vs->numNodes = numNodes;
    vs->size = 0;
    vs->type = type;
    vs->sparseUpToDate = false;
    vs->denseUpToDate = false;
    vs->capacity = capacity;
    if (type == SPARSE)
    {
        vs->vertices = new Vertex[capacity];
        vs->denseVertices = NULL;
        vs->sparseUpToDate = true;

    } else {
        vs->denseVertices = new bool[numNodes];
        vs->vertices = NULL;
        vs->denseUpToDate = true;
    }
    return vs;
}

VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes, bool* denseVertices) {
    VertexSet* vs = new VertexSet;
    vs->numNodes = numNodes;
    vs->size = capacity;
    vs->capacity = capacity;
    vs->numNodes = numNodes;
    vs->type = DENSE;
    vs->denseUpToDate = true;
    vs->sparseUpToDate = false;
    vs->denseVertices = denseVertices;
    vs->vertices = NULL;
    return vs;
}

VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes, Vertex* vertices, int size) {
    VertexSet* vs = new VertexSet;
    vs->numNodes = numNodes;
    vs->size = size;
    vs->capacity = capacity;
    vs->numNodes = numNodes;
    vs->type = SPARSE;
    vs->denseUpToDate = false;
    vs->sparseUpToDate = true;
    vs->vertices = vertices;
    vs->denseVertices = NULL;
    return vs;
}

void freeVertexSet(VertexSet *set)
{

    // delete set->vertices;
    delete[] set->denseVertices;
    delete[] set->vertices;
    delete set;
}

void parallel_update_dense( bool* dense, Vertex* sparse,  int size, int numNodes) {

    // #pragma omp parallel for
    for (int i = 0; i < numNodes; ++i)
    {
        dense[i] = false;
    }
    #pragma omp parallel for
    for (int i = 0; i < size ; ++i)
    {
        dense[sparse[i]] = true;
    }
}

void parallel_pack_scan(Vertex* sparse, bool* dense, int size, int numNodes) {

    Vertex* sums = new Vertex[numNodes];
    prefix_sum(sums, dense, numNodes);
    #pragma omp parallel for
    for (int i = 0; i < numNodes; ++i)
    {
        if (dense[i])
        {
            sparse[sums[i]] =i;
        }
    }
    delete[] sums;
}


void updateDense(VertexSet *set, bool convert = false) {
    if (set->denseVertices == NULL)
    {
        set->denseVertices = new bool[set->numNodes];

    }
    if (!(set->denseUpToDate))
        parallel_update_dense(set->denseVertices, set->vertices, set->size, set->numNodes);
    // set->denseVertices[0];

    if (convert)
        set->type = DENSE;
    set->denseUpToDate = true;
}

void updateSparse(VertexSet *set, bool convert = false ) {
    if (set->vertices == NULL)
    {
        set->vertices = new Vertex[set->capacity];

    }
    if (!(set->sparseUpToDate))
        parallel_pack_scan(set->vertices, set->denseVertices,  set->size, set->numNodes);

    if (convert)
        set->type = SPARSE;

    set->sparseUpToDate = true;
}

void addVertex(VertexSet *set, Vertex v)
{   //this will convert it into dense

    // printf("addVertex\n");
    updateDense(set, true);
    if (set->type == DENSE)
    {
        if (set->denseVertices[v] == false)
        {
            set->size++;
            ENSURES(set->size <= set->capacity);
            set->denseVertices[v] = true;

            set->denseUpToDate = true;
            set->sparseUpToDate = false;
        }
        return;
    }
    else {
        updateDense(set);
        if (set->denseVertices[v] == false) {
            set->size++;
            ENSURES(set->size <= set->capacity);
            set->denseVertices[v] = true;
            set->denseUpToDate = true;
            set->sparseUpToDate = false;
            updateSparse(set);
        }

        return;
    }
}

void removeVertex(VertexSet *set, Vertex v)
{
    ENSURES(set->size > 0);
    //this will convert it into dense
    updateDense(set, true);
    if (set->type == DENSE) {

        if (set->denseVertices[v] == true) {

            set->size--;
            set->denseVertices[v] = false;
            set->denseUpToDate = true;
            set->sparseUpToDate = false;
        }

        return;
    } else {
        updateDense(set);
        if (set->denseVertices[v] == true) {
            set->size--;
            set->denseVertices[v] = false;
            set->denseUpToDate = true;
            set->sparseUpToDate = false;
            updateSparse(set);
        }
        return;
    }

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

