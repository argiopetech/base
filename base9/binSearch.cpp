#include <stdio.h>
#include <stdlib.h>

/*
// last update, 18sep05
//
// performs a binary search for interpolation
// returns the lower index of search
// higher index is simply (lower+1)
// some code taken from http://download.cdsoft.co.uk/tutorials/binarysearch/
*/

int binarySearch (double *searchArray, int size, double searchItem)
{
    // lower bounds check
    if (searchItem <= searchArray[0])
    {
        return 0;
    }
    // upper bounds check
    else if (searchItem >= searchArray[size - 1])
    {
        return size - 2;
    }

    int lo = 0;
    int hi = size;
    int mid;

    while (1)
    {
        /* Using a binary shift right because it's faster than divide */
        mid = ((lo + hi) >> 1);

        if (searchArray[mid] <= searchItem && searchItem <= searchArray[mid + 1])
            return mid;

        // not found, but should never reach here
        if (lo >= hi)
        {
            printf ("ERROR: BINARY SEARCH FAILURE \n");
            return -1;
        }

        if (searchItem > searchArray[mid])
            lo = mid + 1;
        if (searchItem < searchArray[mid])
            hi = mid - 1;
    }

    printf ("ERROR: BINARY SEARCH FAILURE \n");
    return -1;
}

int reverseBinarySearch (double *searchArray, int size, double searchItem)
{
    // lower bounds check
    if (searchItem >= searchArray[0])
    {
        return 0;
    }
    // upper bounds check
    else if (searchItem <= searchArray[size - 1])
    {
        return size - 2;
    }


    int lo = 0;
    int hi = size;
    int mid;

    while (1)
    {
        /* Using a binary shift right because it's faster than divide */
        mid = ((lo + hi) >> 1);

        if (searchArray[mid] >= searchItem && searchItem >= searchArray[mid + 1])
            return mid;

        // not found, but should never reach here
        if (lo >= hi)
        {
            printf ("ERROR: BINARY SEARCH FAILURE \n");
            return -1;
        }

        if (searchItem < searchArray[mid])
            lo = mid + 1;
        if (searchItem > searchArray[mid])
            hi = mid - 1;
    }

    printf ("ERROR: BINARY SEARCH FAILURE \n");
    return -1;
}
