
/* block.h */
/*
	Template classes Block and DBlock
	Implement adding and deleting items of the same type in blocks.

	If there there are many items then using Block or DBlock
	is more efficient than using 'new' and 'delete' both in terms
	of memory and time since
	(1) On some systems there is some minimum amount of memory
	    that 'new' can allocate (e.g., 64), so if items are
	    small that a lot of memory is wasted.
	(2) 'new' and 'delete' are designed for items of varying size.
	    If all items has the same size, then an algorithm for
	    adding and deleting can be made more efficient.
	(3) All Block and DBlock functions are inline, so there are
	    no extra function calls.

	Differences between Block and DBlock:
	(1) DBlock allows both adding and deleting items,
	    whereas Block allows only adding items.
	(2) Block has an additional operation of scanning
	    items added so far (in the order in which they were added).
	(3) Block allows to allocate several consecutive
	    items at a time, whereas DBlock can add only a single item.

	Note that no constructors or destructors are called for items.

	Example usage for items of type 'MyType':

	///////////////////////////////////////////////////
	#include "block.h"
	#define BLOCK_SIZE 1024
	typedef struct { int a, b; } MyType;
	MyType *ptr, *array[10000];

	...

	Block<MyType> *block = new Block<MyType>(BLOCK_SIZE);

	// adding items
	for (int i=0; i<sizeof(array); i++)
	{
		ptr = block -> New();
		ptr -> a = ptr -> b = rand();
	}

	// reading items
	for (ptr=block->ScanFirst(); ptr; ptr=block->ScanNext())
	{
		printf("%d %d\n", ptr->a, ptr->b);
	}

	delete block;

	...

	DBlock<MyType> *dblock = new DBlock<MyType>(BLOCK_SIZE);
	
	// adding items
	for (int i=0; i<sizeof(array); i++)
	{