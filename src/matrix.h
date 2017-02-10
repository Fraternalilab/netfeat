/*=============================================================================
matrix.h : allocate matrices (with dicontinuous memory allocation)
(C) Jens Kleinjung 2007
Read the COPYING file for license information.
=============================================================================*/

#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>

#include "safe.h"
#include "vector.h"

/*___________________________________________________________________________*/
/* prototypes */
/* integer */
int **alloc_mat2D_int(int **mat2D_int, int x, int y);
void init_mat2D_int(int **mat2D_int, int x, int y, int val);
void free_mat2D_int(int **mat2D_int, int x);
void print_mat2D_int(char *outFileName, int **mat2D_int, int x, int y);
void print_mat2D_int_lowlim(char *outFileName, int **mat2D_int, int lowx, int x, int lowy, int y);
int symmetry_mat2D_int(int **mat2D_int, int x, int y);
int binary_mat2D_int(int **mat2D_int, int x, int y);
void multiply_mat2D_int(int **mat2D_int_C, int **mat2D_int_A, int xA, int yA, int **mat2D_int_B, int xB, int yB);
int trace_mat2D_int(int **mat2D_int, int x);
void add_mat2D_int(int **mat2D_int_sum, int **mat2D_int_A, int **mat2D_int_B, int x, int y);
void copy_mat2D_int(int **mat2D_int_B, int **mat2D_int_A, int x, int y);
void regularise_mat2D_int(int **mat2D_int, int x, int y);

int ***alloc_mat3D_int(int ***mat3D_int, int x, int y, int z);
void init_mat3D_int(int ***mat3D_int, int x, int y, int z, int val);
void free_mat3D_int(int ***mat3D_int, int x, int y);

int ****alloc_mat4D_int(int ****mat4D_int, int w, int x, int y, int z);
void init_mat4D_int(int ****mat4D_int, int w, int x, int y, int z, int val);
void free_mat4D_int(int ****mat4D_int, int w, int x, int y);

/* float */
float **alloc_mat2D_float(float **mat2D_float, int x, int y);
void init_mat2D_float(float **mat2D_float, int x, int y, float val);
void free_mat2D_float(float **mat2D_float, int x);
void print_mat2D_float(char *outFileName, float **mat2D_float, int x, int y);
void print_mat2D_float_lowlim(char *outFileName, float **mat2D_float, int lowx, int x, int lowy, int y);
void print_mat2D_floate_lowlim(char *outFileName, float **mat2D_float, int lowx, int x, int lowy, int y);
void div_mat2D_float(float **mat2D_float, int x, int y, float a);
int symmetry_mat2D_float(float **mat2D_float, int x, int y);
void add_mat2D_float(float **mat2D_int_sum, float **mat2D_int_A, float **mat2D_int_B, int x, int y);

float ***alloc_mat3D_float(float ***mat3D_float, int x, int y, int z);
void init_mat3D_float(float ***mat3D_float, int x, int y, int z, float val);
void free_mat3D_float(float ***mat3D_float, int x, int y);

float ****alloc_mat4D_float(float ****mat4D_float, int w, int x, int y, int z);
void init_mat4D_float(float ****mat4D_float, int w, int x, int y, int z, float val);
void free_mat4D_float(float ****mat4D_float, int w, int x, int y);

/* double */
double **alloc_mat2D_double(double **mat2D_double, int x, int y);
void init_mat2D_double(double **mat2D_double, int x, int y, double val);
void free_mat2D_double(double **mat2D_double, int x);
void print_mat2D_double(char *outFileName, double **mat2D_double, int x, int y);
void print_mat2D_double_lowlim(char *outFileName, double **mat2D_double, int lowx, int x, int lowy, int y);
void print_mat2D_doublee_lowlim(char *outFileName, double **mat2D_double, int lowx, int x, int lowy, int y);
void div_mat2D_double(double **mat2D_double, int x, int y, double a);
int symmetry_mat2D_double(double **mat2D_double, int x, int y);
void copy_mat2D_double(double **mat2D_double_to, double **mat2D_double_from, int x, int y);
void add_mat2D_double(double **mat2D_int_sum, double **mat2D_int_A, double **mat2D_int_B, int x, int y);

double ***alloc_mat3D_double(double ***mat3D_double, int x, int y, int z);
void init_mat3D_double(double ***mat3D_double, int x, int y, int z, double val);
void free_mat3D_double(double ***mat3D_double, int x, int y);

double ****alloc_mat4D_double(double ****mat4D_double, int w, int x, int y, int z);
void init_mat4D_double(double ****mat4D_double, int w, int x, int y, int z, double val);
void free_mat4D_double(double ****mat4D_double, int w, int x, int y);

/* vector */
Vec **alloc_mat2D_vec(Vec **mat2D_vec, int x, int y);
void init_mat2D_vec(Vec **mat2D_vec, int x, int y, Vec val);
void free_mat2D_vec(Vec **mat2D_vec, int x);

Vec ***alloc_mat3D_vec(Vec ***mat3D_vec, int x, int y, int z);
void init_mat3D_vec(Vec ***mat3D_vec, int x, int y, int z, Vec val);
void free_mat3D_vec(Vec ***mat3D_vec, int x, int y);

Vec ****alloc_mat4D_vec(Vec ****mat4D_vec, int w, int x, int y, int z);
void init__mat4D_vec(Vec ****mat4D_vec, int w, int x, int y, int z, Vec val);
void free_mat4D_vec(Vec ****mat4D_vec, int w, int x, int y);

#endif

