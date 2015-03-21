#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

/************************** allocMultidim ***************************
* Argument  1:  Integer size of elements in the array
* Argument  2:  Integer number of dimensions
* Argument 3+:  Integers specifying array size in each dimension
* Allocate memory for a n-dimensional array, with dimensions given in the
* variable list of arguments. Allocate all the memory in one large chunk to
* ensure (in the 2d case) the last element of a row is contigous with the
* first of the next.
*
* Author: Julian Panetta
********************************************************************/
void **allocMultidim(unsigned int type_size, unsigned int nDims, ...)
{
    assert(nDims > 0);

    va_list arguments;

    unsigned int i, j, p = 0;
	unsigned long int nAbove = 1, indexEntries = 0, columnVectors = 1;
    unsigned long int blockSize;
	void **anArray; 
    
    unsigned int *dims = (unsigned int *)malloc(nDims * sizeof(unsigned int));
    va_start(arguments, nDims);

	for(i = 0; i < nDims; i++)	{
        /* read in the array sizes */
        dims[i] = va_arg(arguments, int);
        /* The first (nDims - 1) dimensions determine the number of column vectors */
        if (i < nDims - 1)  {
            columnVectors *= dims[i];
            indexEntries += columnVectors;
        }
	}

    va_end(arguments);

	blockSize = indexEntries * sizeof(void *) +
                columnVectors * dims[nDims - 1] * type_size;

	anArray = (void **)malloc(blockSize);

    if (!anArray)   {
        printf("Allocation of block of size %lu failed!\n", blockSize);
        exit(-1);
    }
	
	for(i = 0; i < nDims - 1; i++)	{	/* build index */
		anArray[p] = (void *)((char *)&anArray[p] + nAbove * dims[i] * sizeof(char *));
		for(j = p + 1; j < p + dims[i] * nAbove; j++)
			anArray[j] = (void *)((char *)anArray[j - 1] +
                dims[i + 1] * ((i < nDims - 2) ? sizeof(char *) : type_size));
		nAbove *= dims[i];
		p += nAbove;
	}

    free(dims);

	return anArray;
}
