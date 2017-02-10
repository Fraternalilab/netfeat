/*==============================================================================
array.h : array routines
Copyright (C) 2009 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARRAY_H
#define ARRAY_H

#include "safe.h"

/*___________________________________________________________________________*/
/* prototypes */
/* INTEGER */
void init_array_int(int *array, int length, int val);
void print_array_int(char *outFileName, int *array, int length);
void print_arrays_int(char *outFileName, int *array0, int *array1, int length);
void copy_arrays_int(int *array_to, int *array_from, int length);
void add_array_int(int *array_sum, int *array_A, int *array_B, int length);
/* FLOAT */
void init_array_float(float *array, int length, float val);
void print_array_float(char *outFileName, float *array, int length);
void print_arrays_float(char *outFileName, float *array0, float *array1, int length);
void copy_arrays_float(float *array_to, float *array_from, int length);
void add_array_float(float *array_sum, float *array_A, float *array_B, int length);
/* DOUBLE */
void init_array_double(double *array, int length, double val);
void print_array_double(char *outFileName, double *array, int length);
void print_arrays_double(char *outFileName, double *array0, double *array1, int length);
void copy_arrays_double(double *array_to, double *array_from, int length);
void add_array_double(double *array_sum, double *array_A, double *array_B, int length);
#endif
