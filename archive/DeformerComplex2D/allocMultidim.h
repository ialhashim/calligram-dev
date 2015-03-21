#ifndef allocMultidim_included
#define allocMultidim_included

/************************** allocMultidim ***************************
* Allocate memory for a n-dimensional array, with dimensions given in
* dims. Allocate all the memory in one large chunk to ensure (in the 
* 2d case) the last element of a row is contigous with the first of
* the next.
********************************************************************/
void **allocMultidim(unsigned int type_size, unsigned int n, ...);

#endif
