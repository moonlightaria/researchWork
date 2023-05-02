#ifndef HEAP_H 
#define HEAP_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "crack.h"
/*
#ifndef XD
#define XD 300
#define YD 300
#endif
*/


#define BLEH 7

#define BIGNUM 320001

typedef struct {
	long int distance;
	long int temp;
	int id;
} strippednode;


typedef struct {
	strippednode entry[bXD*bYD*BLEH];
	int numentries;
	int cutoff;
} heaptype;

heaptype heap, focal;
//
//int testheap(void);
//void emptyheap(void);
//strippednode maxfromheap(void);
//void printheap(void);
//void insertintoheap(strippednode innode);

/*----------------------------------------------------------------------*/
int testheap(void)
{
	int i;
	int parent;
	for (i = 0; i != heap.numentries; i++)
	{
		if (i > 0)
			parent = (i - 1) / 2;
		else
			parent = 0;

		if (heap.entry[i].distance < heap.entry[parent].distance)
		{
			fprintf(stderr, "Inconsistent heap, %i > %i (%i tot)\n",
				parent, i, heap.numentries);
			return(1);
		}
	}
	return(0);
}
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
void emptyheap(void)
{
	// throw out contents of heap, if any

	heap.numentries = 0;
}
/*----------------------------------------------------------------------*/
strippednode maxfromheap(void)
{
	// remove the max from the heap and fix the heap
	strippednode retval;
	int i;
	strippednode top;
	int curindex;
	int rightchild, leftchild, largerchild;

	retval = heap.entry[0];

	// put last entry at top of heap:
	heap.numentries--;
	if (heap.numentries < 0)
	{
		fprintf(stderr, "Tried to pop from empty heap, crash\n");
		exit(0);
	}
	heap.entry[0] = heap.entry[heap.numentries];
	top = heap.entry[0];

	// trickle down:
	curindex = 0;
	while (curindex < heap.numentries / 2)
	{
		leftchild = 2 * curindex + 1;
		rightchild = leftchild + 1;

		if ((rightchild < heap.numentries) &&
			(heap.entry[leftchild].distance > // hi pri to lower cost
			heap.entry[rightchild].distance))
		{
			largerchild = rightchild;
		}
		else
		{
			largerchild = leftchild;
		}

		if (top.distance <= heap.entry[largerchild].distance)
			break; // done
		else
		{
			heap.entry[curindex] = heap.entry[largerchild];
			curindex = largerchild;
		}

	} // end while

	heap.entry[curindex] = top;

	return(retval);

}
/*----------------------------------------------------------------------*/
void printheap(void)
{
	int i;
	for (i = 0; i != heap.numentries; i++)
	{
		fprintf(stderr, "%i(%li)..", heap.entry[i].id, heap.entry[i].distance);
	}
	fprintf(stderr, "\n");
}
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
void insertintoheap(strippednode innode)
{

	// insert innode into heap. First place at end, then trickle up.

	int curindex;
	int parent;
	strippednode bottom;
	//int found;

	if (heap.numentries == XD*YD*BLEH)
	{
		// heap filled!
		fprintf(stderr, "Heap full!\n");
		exit(0);
	}

	heap.entry[heap.numentries] = innode;

	bottom = innode;

	curindex = heap.numentries;
	heap.numentries++;
	parent = (curindex - 1) / 2;

	while ((curindex > 0) &&
		(heap.entry[parent].distance >
		bottom.distance))
	{
		heap.entry[curindex] = heap.entry[parent];
		curindex = parent;
		parent = (parent - 1) / 2;
	}

	heap.entry[curindex] = bottom;

}

#endif
