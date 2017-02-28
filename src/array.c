/*==============================================================================
array.c : array routines
(C) 2010 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "array.h"

/*___________________________________________________________________________*/
/** INTEGER */
/*____________________________________________________________________________*/
/** initialise integer array */
void init_array_int(int *array, int length, int val)
{
    unsigned int i;

    for (i = 0; i < length; ++ i)
        array[i] = val;
}

/*____________________________________________________________________________*/
/** print integer array */
void print_array_int(char *outFileName, int *array, int length)
{
    FILE *outFile = 0;
    unsigned int i;

    outFile = safe_open(outFileName, "w");

    /* print array */
    for (i = 0; i < length; ++ i)
        fprintf(outFile, "%d\t%d\n", i, array[i]);
    
    fclose(outFile);
}

/*____________________________________________________________________________*/
/** print two integer arrays of same length */
void print_arrays_int(char *outFileName, int *array0, int *array1, int length)
{
    FILE *outFile = 0;
    unsigned int i;

    outFile = safe_open(outFileName, "w");

    /* print array */
    for (i = 0; i < length; ++ i)
        fprintf(outFile, "%d\t%d\t%d\n", i, array0[i], array1[i]);
    
    fclose(outFile);
}

/*____________________________________________________________________________*/
/** copy array of same length */
void copy_array_int(int *array_to, int *array_from, int length)
{
    unsigned int i;

    /* copy array */
    for (i = 0; i < length; ++ i)
        array_to[i] = array_from[i];
}

/*____________________________________________________________________________*/
/** add array */
void add_array_int(int *array_sum, int *array_A, int *array_B, int length)
{
    unsigned int i;

    for (i = 0; i < length; ++ i)
        array_sum[i] = array_A[i] + array_B[i];
}

/*___________________________________________________________________________*/
/** FLOAT */
/*____________________________________________________________________________*/
/** initialise float array */
void init_array_float(float *array, int length, float val)
{
    unsigned int i;

    for (i = 0; i < length; ++ i)
        array[i] = val;
}

/*____________________________________________________________________________*/
/** print float array */
void print_array_float(char *outFileName, float *array, int length)
{
    FILE *outFile = 0;
    unsigned int i;

    outFile = safe_open(outFileName, "w");

    /* print array */
    for (i = 0; i < length; ++ i)
        fprintf(outFile, "%d\t%f\n", i, array[i]);
    
    fclose(outFile);
}

/*____________________________________________________________________________*/
/** print two arrays of same length */
void print_arrays_float(char *outFileName, float *array0, float *array1, int length)
{
    FILE *outFile = 0;
    unsigned int i;

    outFile = safe_open(outFileName, "w");

    /* print array */
    for (i = 0; i < length; ++ i)
        fprintf(outFile, "%d\t%f\t%f\n", i, array0[i], array1[i]);
    
    fclose(outFile);
}

/*____________________________________________________________________________*/
/** copy array of same length */
void copy_array_float(float *array_to, float *array_from, int length)
{
    unsigned int i;

    /* copy array */
    for (i = 0; i < length; ++ i)
        array_to[i] = array_from[i];

}

/*____________________________________________________________________________*/
/** add array of same length */
void add_array_float(float *array_sum, float *array_A, float *array_B, int length)
{
    unsigned int i;

    for (i = 0; i < length; ++ i)
        array_sum[i] = array_A[i] + array_B[i];
}

/*___________________________________________________________________________*/
/** DOUBLE */
/*____________________________________________________________________________*/
/** initialise double array */
void init_array_double(double *array, int length, double val)
{
    unsigned int i;

    for (i = 0; i < length; ++ i)
        array[i] = val;
}

/*____________________________________________________________________________*/
/** print double array */
void print_array_double(char *outFileName, double *array, int length)
{
    FILE *outFile = 0;
    unsigned int i;

    outFile = safe_open(outFileName, "w");

    /* print array */
    for (i = 0; i < length; ++ i)
        fprintf(outFile, "%d\t%16.8e\n", i, array[i]);
    
    fclose(outFile);
}

/*____________________________________________________________________________*/
/** print two arrays of same length */
void print_arrays_double(char *outFileName, double *array0, double *array1, int length)
{
    FILE *outFile = 0;
    unsigned int i;

    outFile = safe_open(outFileName, "w");

    /* print array */
    for (i = 0; i < length; ++ i)
        fprintf(outFile, "%d\t%lf\t%lf\n", i, array0[i], array1[i]);
    
    fclose(outFile);
}


/*____________________________________________________________________________*/
/** copy array of same length */
void copy_arrays_double(double *array_to, double *array_from, int length)
{
    unsigned int i;

    /* copy array */
    for (i = 0; i < length; ++ i)
        array_to[i] = array_from[i];

}

/*____________________________________________________________________________*/
/** add array of same length */
void add_array_double(double *array_sum, double *array_A, double *array_B, int length)
{
    unsigned int i;

    for (i = 0; i < length; ++ i)
        array_sum[i] = array_A[i] + array_B[i];
}

