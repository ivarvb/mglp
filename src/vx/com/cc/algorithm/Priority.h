

#include <stdio.h>
#include <stdlib.h>

typedef struct HeapPair HeapPair;
typedef struct Heap Heap;
Heap* new_heap(int isminheap, int initial_size);
void free_Heap(Heap* h);
void max_heapify(HeapPair *data, int loc, int count);
void heap_push(Heap *h, HeapPair e);
void heap_display(Heap *h);
HeapPair heap_pop(Heap *h);


typedef struct HeapPair {
	double weight;
	int i;
	int k;
}HeapPair;


typedef struct Heap {
	int size;
	int count;
    int initial_size;
    double factor;
    int isminheap;
	HeapPair *heaparr;
}Heap;

//int *heap, size, count;
Heap* new_heap(int isminheap, int initial_size){
    Heap *h = (Heap*)malloc(sizeof(Heap));
    h->initial_size = initial_size;
    h->factor = 1.0;
    h->isminheap = isminheap;
    if (h->isminheap==1){
        h->factor = -1.0;
    }
	h->count = 0;
	h->size = h->initial_size;
	h->heaparr = (HeapPair *) malloc(sizeof(HeapPair) * h->initial_size);
	if(!h->heaparr) {printf("Error allocatinga memory...\n");exit(-1);}
    return h;
}

void free_Heap(Heap* h){
    free(h->heaparr);
    free(h);
}

void max_heapify(HeapPair *data, int loc, int count) {
    HeapPair temp;
	int left, right, largest;
	left = 2*(loc) + 1;
	right = left + 1;
	largest = loc;
	

	if (left <= count && data[left].weight > data[largest].weight) {
		largest = left;
	} 
	if (right <= count && data[right].weight > data[largest].weight) {
		largest = right;
	} 
	
	if(largest != loc) {
		temp = data[loc];
		data[loc] = data[largest];
		data[largest] = temp;
		max_heapify(data, largest, count);
	}

}

void heap_push(Heap *h, HeapPair e){

	//double value, int id
    e.weight = e.weight*h->factor;

	int index, parent;
 
	// Resize the heap if it is too small to hold all the data
	if (h->count == h->size){
		h->size += 1;
		h->heaparr = (HeapPair*) realloc(h->heaparr, sizeof(HeapPair) * h->size);
		if (!h->heaparr) exit(-1); // Exit if the memory allocation fails
	}
 	
 	index = h->count++; // First insert at last of array

 	// Find out where to put the element and put it
	for(;index; index = parent){
		parent = (index - 1) / 2;
		if (h->heaparr[parent].weight >= e.weight) break;
		h->heaparr[index] = h->heaparr[parent];
	}
	h->heaparr[index].weight = e.weight;
	h->heaparr[index].i = e.i;
	h->heaparr[index].k = e.k;
}

void heap_display(Heap *h) {
	int i;
	for(i=0; i<h->count; ++i) {
		printf("|%f, %d|", h->heaparr[i].weight, h->heaparr[i].i );
	}
	printf("\n");
}

HeapPair heap_pop(Heap *h){
	HeapPair removed;
	HeapPair temp = h->heaparr[--h->count];
 	
	
	if ((h->count <= (h->size + 2)) && (h->size > h->initial_size)){
		h->size -= 1;
		h->heaparr = (HeapPair*) realloc(h->heaparr, sizeof(HeapPair) * h->size);
		if (!h->heaparr) exit(-1); // Exit if the memory allocation fails
	}
 	removed = h->heaparr[0];
 	h->heaparr[0] = temp;
 	max_heapify(h->heaparr, 0, h->count);
 	return removed;
}


// int emptyPQ(Heap *pq) {
// 	int i;
//     HeapPair aux;
// 	while(pq->count != 0) {
//         aux = heap_pop(pq);
// 		printf("%f, %d, %d \n", aux.weight, aux.i, aux.k);
// 	}
// }

