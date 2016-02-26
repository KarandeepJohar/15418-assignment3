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
    #pragma omp parallel for default(none) shared(N, output, boolArray)
    for (int i = 0; i < N; i++) {
        if (boolArray[i]) {
            output[i] = 1;
        } else {
            output[i] = 0;
        }
    }
    int *arr, *partial, *temp;
    int num_threads, work, n;
    int i, mynum, last;
    arr = output;
    n = N;
    #pragma omp parallel default(none) private(i, mynum, last) shared(arr, partial, temp, num_threads, work, n)
    {
        #pragma omp single
        {
            num_threads = omp_get_num_threads();
            partial = (int *) malloc (sizeof (int) * num_threads);
            temp = (int *) malloc (sizeof (int) * num_threads);
            work = n / num_threads + 1; /*sets length of sub-arrays*/
        }
        mynum = omp_get_thread_num();
        /*calculate prefix-sum for each subarray*/
        for (i = work * mynum + 1; i < work * mynum + work && i < n; i++)
            arr[i] += arr[i - 1];
        partial[mynum] = arr[i - 1];
        #pragma omp barrier
        /*calculate prefix sum for the array that was made from last elements of each of the previous sub-arrays*/
        for (i = 1; i < num_threads; i <<= 1) {
            if (mynum >= i)
                temp[mynum] = partial[mynum] + partial[mynum - i];
            #pragma omp barrier
            #pragma omp single
            memcpy(partial + 1, temp + 1, sizeof(int) * (num_threads - 1));
        }
        for (i = work * mynum; i < (last = work * mynum + work < n ? work * mynum + work : n); i++)
            arr[i] += partial[mynum] - arr[last - 1];
    }

}

void prefix_sum(Vertex* output, int* boolArray, int N) {
    int *arr, *partial, *temp;
    int num_threads, work, n;
    int i, mynum, last;
    arr = boolArray;
    n = N;
    #pragma omp parallel default(none) private(i, mynum, last) shared(arr, partial, temp, num_threads, work, n)
    {
        #pragma omp single
        {
            num_threads = omp_get_num_threads();
            partial = (int *) malloc (sizeof (int) * num_threads);
            temp = (int *) malloc (sizeof (int) * num_threads);
            work = n / num_threads + 1; /*sets length of sub-arrays*/
        }
        mynum = omp_get_thread_num();
        /*calculate prefix-sum for each subarray*/
        for (i = work * mynum + 1; i < work * mynum + work && i < n; i++)
            arr[i] += arr[i - 1];
        partial[mynum] = arr[i - 1];
        #pragma omp barrier
        /*calculate prefix sum for the array that was made from last elements of each of the previous sub-arrays*/
        for (i = 1; i < num_threads; i <<= 1) {
            if (mynum >= i)
                temp[mynum] = partial[mynum] + partial[mynum - i];
            #pragma omp barrier
            #pragma omp single
            memcpy(partial + 1, temp + 1, sizeof(int) * (num_threads - 1));
        }
        for (i = work * mynum; i < (last = work * mynum + work < n ? work * mynum + work : n); i++)
            arr[i] += partial[mynum] - arr[last - 1];
    }

}

int packIndices(Vertex* output, Vertex* input, bool* boolArray, int n) {
    Vertex* sums = new Vertex[n + 1];
    sums[0] = 0;
    prefix_sum(sums + 1, boolArray, n);
    int sum = 0;
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


    static Vertex* flags = NULL;
    if (!flags) {
        flags = new Vertex[numNodes];
        #pragma omp parallel for
        for (int i = 0; i < numNodes; ++i)
        {
            flags[i] = -1;
        }
    }
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for (int i = 0; i < size; ++i)
        {
            if (input[i] != -1 && flags[input[i]] == -1)
            {
                __sync_bool_compare_and_swap(&flags[input[i]], -1, i);
            }
        }
        #pragma omp for schedule(static)
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
    }
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
    vs->denseVertices = NULL;
    vs->sparseUpToDate = false;
    vs->denseUpToDate = false;

    if (type == SPARSE)
    {
        vs->sparseUpToDate = true;
    } else {
        vs->denseUpToDate = true;
    }
    vs->vertices = NULL;
    vs->sum_degrees = 0;
    return vs;
}

VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes, bool* denseVertices, int sum_degrees=0) {
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
    vs->sum_degrees = sum_degrees;
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
    vs->sum_degrees = 0;
    return vs;
}

void freeVertexSet(VertexSet *set)
{
    delete[] set->denseVertices;
    delete[] set->vertices;
    delete set;
}

void parallel_update_dense( bool* dense, Vertex* sparse,  int size, int numNodes) {

    #pragma omp parallel for
    for (int i = 0; i < size ; ++i)
    {
        dense[sparse[i]] = true;
    }
}

void parallel_pack_scan(Vertex* sparse, bool* dense, int size, int numNodes) {

    static Vertex* sums = new Vertex[numNodes + 1];
    sums[0] = 0;
    prefix_sum(sums + 1, dense, numNodes);
    #pragma omp parallel for
    for (int i = 0; i < numNodes; ++i)
    {
        if (dense[i])
        {
            sparse[sums[i]] = i;
        }
    }
}


void updateDense(VertexSet *set, bool convert = false) {
    if (set->denseVertices == NULL)
    {
        set->denseVertices = new bool[set->numNodes]();

    }
    if (!(set->denseUpToDate)) {
        parallel_update_dense(set->denseVertices, set->vertices, set->size, set->numNodes);
    }

    if (convert)
        set->type = DENSE;
    set->denseUpToDate = true;
}

void updateSparse(VertexSet *set, bool convert = false ) {
    if (set->vertices == NULL)
    {
        set->vertices = new Vertex[set->capacity];

    }

    if (!(set->sparseUpToDate)) {
        parallel_pack_scan(set->vertices, set->denseVertices,  set->size, set->numNodes);
    }
    if (convert)
        set->type = SPARSE;

    set->sparseUpToDate = true;
}

void addVertex(VertexSet *set, Vertex v)
{   //this will convert it into dense

    updateDense(set, true);
    if (set->capacity>set->numNodes/100 || set->type ==DENSE)
    {
        if (set->denseVertices[v] == false)
        {
            set->size++;
            set->denseVertices[v] = true;
            set->denseUpToDate = true;
            set->sparseUpToDate = false;
        }
        return;
    }
    else {
        updateSparse(set, true);
        ENSURES(set->size <= set->capacity);
        set->vertices[set->size++] = v;
        set->denseUpToDate = false;
        set->sparseUpToDate = true;
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
    updateDense(u,true);
    updateDense(v,true);
    bool* newDense = new bool[u->numNodes]();
    int size = 0;
    #pragma omp parallel for reduction(+:size)
    for (int i=0; i < u->numNodes; i++) {
        if (u->denseVertices[i] || v->denseVertices[i]) {
            newDense[i] = true;
            size += 1;
        }
    }
    VertexSet* vs = newVertexSet(DENSE, size, u->numNodes, newDense, 0);

    freeVertexSet(u);
    freeVertexSet(v);
    return vs;
}

