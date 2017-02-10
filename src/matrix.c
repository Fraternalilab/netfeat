/*=============================================================================
matrix.c : allocate matrices with dicontinuous memory allocation
(C) Jens Kleinjung 2007
Read the COPYING file for license information.
==============================================================================*/

#include "matrix.h"

/*___________________________________________________________________________*/
/** INTEGER */
/*___________________________________________________________________________*/
/** 2D integer matrix */
/** allocate */
int **alloc_mat2D_int(int **mat2D_int, int x, int y)
{
    unsigned int i;

	mat2D_int = (int **)safe_malloc(x * sizeof(int *));
	for (i = 0; i < x; ++ i)
		mat2D_int[i] = (int *)safe_malloc(y * sizeof(int));

	return mat2D_int;
}

/** initialise */
void init_mat2D_int(int **mat2D_int, int x, int y, int val)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			mat2D_int[i][j] = val;
}

/** free */
void free_mat2D_int(int **mat2D_int, int x)
{
    unsigned int i;

	for (i = 0; i < x; ++ i)
		free(mat2D_int[i]);

    free(mat2D_int);
}

/** print */
void print_mat2D_int(char *outFileName, int **mat2D_int, int x, int y)
{
    unsigned int i, j;
    FILE *outFile = 0;

    outFile = safe_open(outFileName, "w");

	for (i = 0; i < x; ++ i) {
		if (i > 0) fprintf(outFile, "\n");
		for (j = 0; j < y; ++ j) {
			fprintf(outFile, "%2d ", mat2D_int[i][j]);
		}
	}
	fprintf(outFile, "\n");

	fclose(outFile);
}

/** print with lower limit specification */
void print_mat2D_int_lowlim(char *outFileName, int **mat2D_int, int lowx, int x, int lowy, int y)
{
    unsigned int i, j;
    FILE *outFile = 0;

    outFile = safe_open(outFileName, "w");

	for (i = lowx; i < x; ++ i) {
		if (i > lowx) fprintf(outFile, "\n");
		for (j = lowy; j < y; ++ j) {
			fprintf(outFile, "%2d ", mat2D_int[i][j]);
		}
	}
	fprintf(outFile, "\n");

	fclose(outFile);
}

/** symmetry check */
int symmetry_mat2D_int(int **mat2D_int, int x, int y)
{
    unsigned int i, j, s;

	for (i = 0, s = 0; i < x; ++ i)
		for (j = 0; j < i; ++ j)
			s += (mat2D_int[i][j] == mat2D_int[j][i] ? 0 : 1);

	return s;
}

/** binarity check */
int binarity_mat2D_int(int **mat2D_int, int x, int y)
{
    unsigned int i, j, s;

	for (i = 0, s = 0; i < x; ++ i)
		for (j = 0; j < i; ++ j)
			s += (((mat2D_int[i][j] == 0) || (mat2D_int[i][j] == 1)) ? 0 : 1);

	return s;
}

/* multiply matrix */
void multiply_mat2D_int(int **mat2D_int_C, int **mat2D_int_A, int xA, int yA, \
				   int **mat2D_int_B, int xB, int yB)
{
	unsigned int i, j, k;

	assert(yA == xB);

	for (i = 0; i < xA; ++ i) 
		for (j = 0; j < yB; ++ j) {
				mat2D_int_C[i][j] = 0;
				for (k = 0; k < yA; ++ k) 
					mat2D_int_C[i][j] += mat2D_int_A[i][k] * mat2D_int_B[k][j];
		}
}

/* matrix trace */
int trace_mat2D_int(int **mat2D_int, int x)
{
	unsigned int i;
	int trace = 0;

    for (i = 0; i < x; ++ i) 
		trace += mat2D_int[i][i];

	return trace;	
}

/** add matrix */
void add_mat2D_int(int **mat2D_int_sum, int **mat2D_int_A, int **mat2D_int_B, int x, int y)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			mat2D_int_sum[i][j] = mat2D_int_A[i][j] + mat2D_int_B[i][j];
}

/* copy matrix */
void copy_mat2D_int(int **mat2D_int_B, int **mat2D_int_A, int x, int y)
{
	unsigned int i, j;

    for (i = 0; i < x; ++ i) 
		for (j = 0; j < y; ++ j)
			mat2D_int_B[i][j] = mat2D_int_A[i][j];
}

/* regularise matrix */
/* this is a special function applied after interaction matrix self-multiplication :
 * set all values >0 to 1 and all diagonal values to zero */
void regularise_mat2D_int(int **mat2D_int, int x, int y)
{
	unsigned int i, j;

	for (i = 0; i < x; ++ i) {
		mat2D_int[i][i] = 0;
		for (j = 0; j < y; ++ j) {
			if (mat2D_int[i][j] > 0)
				mat2D_int[i][j] = 1;
		}
	}
}

/*___________________________________________________________________________*/
/** 3D integer matrix */
/** allocate */
int ***alloc_mat3D_int(int ***mat3D_int, int x, int y, int z)
{
	unsigned int i, j;

	mat3D_int = (int ***)safe_malloc(x * sizeof(int **));
	for (i = 0; i < x; ++ i) {
		mat3D_int[i] = (int **)safe_malloc(y * sizeof(int *));
		for (j = 0; j < y; ++ j)
			mat3D_int[i][j] = (int *)safe_malloc(z * sizeof(int));
	}

	return mat3D_int;
}

/** initialise */
void init_mat3D_int(int ***mat3D_int, int x, int y, int z, int val)
{
    unsigned int i, j, k;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			for (k = 0; k < z; ++ k)
				mat3D_int[i][j][z] = val;
}

/** free */
void free_mat3D_int(int ***mat3D_int, int x, int y)
{
    unsigned int i, j;

    for (i = 0; i < x; ++ i) {
        for (j = 0; j < y; ++ j)
            free(mat3D_int[i][j]);
        free(mat3D_int[i]);
    }
    free(mat3D_int);
}

/*___________________________________________________________________________*/
/** 4D integer matrix */
/** allocate */
int ****alloc_mat4D_int(int ****mat4D_int, int w, int x, int y, int z)
{
	unsigned int h, i, j;

	mat4D_int = (int ****)safe_malloc(w * sizeof(int ***));
	for (h = 0; h < w; ++ h) {
		mat4D_int[h] = (int ***)safe_malloc(x * sizeof(int **));
		for (i = 0; i < x; ++ i) {
			mat4D_int[h][i] = (int **)safe_malloc(y * sizeof(int *));
			for (j = 0; j < y; ++ j)
				mat4D_int[h][i][j] = (int *)safe_malloc(z * sizeof(int));
		}
	}

	return mat4D_int;
}

/** initialise */
void init_mat4D_int(int ****mat4D_int, int w, int x, int y, int z, int val)
{
    unsigned int h, i, j, k;

	for (h = 0; h < w; ++ h)
		for (i = 0; i < x; ++ i)
			for (j = 0; j < y; ++ j)
				for (k = 0; k < z; ++ k)
					mat4D_int[h][i][j][k] = val;
}

/** free */
void free_mat4D_int(int ****mat4D_int, int w, int x, int y)
{
    unsigned int h, i, j;

    for (h = 0; h < w; ++ h) {
		for (i = 0; i < x; ++ i) {
			for (j = 0; j < y; ++ j)
				free(mat4D_int[h][i][j]);
			free(mat4D_int[h][i]);
	    }
		free(mat4D_int[h]);
	}
	free(mat4D_int);
}

/*___________________________________________________________________________*/
/** FLOAT */
/*___________________________________________________________________________*/
/** 2D float matrix */
/** allocate */
float **alloc_mat2D_float(float **mat2D_float, int x, int y)
{
    unsigned int i;

	mat2D_float = (float **)safe_malloc(x * sizeof(float *));
	for (i = 0; i < x; ++ i)
		mat2D_float[i] = (float *)safe_malloc(y * sizeof(float));

	return mat2D_float;
}

/** initialise */
void init_mat2D_float(float **mat2D_float, int x, int y, float val)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			mat2D_float[i][j] = val;
}

/* free */
void free_mat2D_float(float **mat2D_float, int x)
{
    unsigned int i;

	for (i = 0; i < x; ++ i)
		free(mat2D_float[i]);
    free(mat2D_float);
}

/** print */
void print_mat2D_float(char *outFileName, float **mat2D_float, int x, int y)
{
    unsigned int i, j;
    FILE *outFile = 0;

    outFile = safe_open(outFileName, "w");

	for (i = 0; i < x; ++ i) {
		if (i > 0) fprintf(outFile, "\n");
		for (j = 0; j < y; ++ j) {
			fprintf(outFile, "%3.2f ", mat2D_float[i][j]);
		}
	}
	fprintf(outFile, "\n");

	fclose(outFile);
}

/** print with lower limit specification */
void print_mat2D_float_lowlim(char *outFileName, float **mat2D_float, int lowx, int x, int lowy, int y)
{
    unsigned int i, j;
    FILE *outFile = 0;

    outFile = safe_open(outFileName, "w");

	for (i = lowx; i < x; ++ i) {
		if (i > lowx) fprintf(outFile, "\n");
		for (j = lowy; j < y; ++ j) {
			fprintf(outFile, "%3.2f ", mat2D_float[i][j]);
		}
	}
	fprintf(outFile, "\n");

	fclose(outFile);
}

/** print with lower limit specification */
void print_mat2D_floate_lowlim(char *outFileName, float **mat2D_float, int lowx, int x, int lowy, int y)
{
    unsigned int i, j;
    FILE *outFile = 0;

    outFile = safe_open(outFileName, "w");

	for (i = lowx; i < x; ++ i) {
		if (i > lowx) fprintf(outFile, "\n");
		for (j = lowy; j < y; ++ j) {
			fprintf(outFile, "%10.5e ", mat2D_float[i][j]);
		}
	}
	fprintf(outFile, "\n");

	fclose(outFile);
}

/** divide by factor */
void div_mat2D_float(float **mat2D_float, int x, int y, float a)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i) {
		for (j = 0; j < y; ++ j) {
			mat2D_float[i][j] /= a;
		}
	}
}

/** symmetry check */
int symmetry_mat2D_float(float **mat2D_float, int x, int y)
{
    unsigned int i, j, s;

	for (i = 0, s = 0; i < x; ++ i)
		for (j = 0; j < i; ++ j)
			s += (mat2D_float[i][j] == mat2D_float[j][i] ? 0 : 1);

	return s;
}

/** add matrix */
void add_mat2D_float(float **mat2D_float_sum, float **mat2D_float_A, float **mat2D_float_B, int x, int y)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			mat2D_float_sum[i][j] = mat2D_float_A[i][j] + mat2D_float_B[i][j];
}

/*___________________________________________________________________________*/
/** 3D float matrix */
/** allocate */
float ***alloc_mat3D_float(float ***mat3D_float, int x, int y, int z)
{
	unsigned int i, j;

	mat3D_float = (float ***)safe_malloc(x * sizeof(float **));
	for (i = 0; i < x; ++ i) {
		mat3D_float[i] = (float **)safe_malloc(y * sizeof(float *));
		for (j = 0; j < y; ++ j)
			mat3D_float[i][j] = (float *)safe_malloc(z * sizeof(float));
	}

	return mat3D_float;
}

/** initialise */
void init_mat3D_float(float ***mat3D_float, int x, int y, int z, float val)
{
    unsigned int i, j, k;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			for (k = 0; k < z; ++ k)
				mat3D_float[i][j][z] = val;
}


/** free */
void free_mat3D_float(float ***mat3D_float, int x, int y)
{
    unsigned int i, j;

    for (i = 0; i < x; ++ i) {
        for (j = 0; j < y; ++ j)
            free(mat3D_float[i][j]);
        free(mat3D_float[i]);
    }
    free(mat3D_float);
}

/*___________________________________________________________________________*/
/** 4D float matrix */
/** allocate */
float ****alloc_mat4D_float(float ****mat4D_float, int w, int x, int y, int z)
{
	unsigned int h, i, j;

	mat4D_float = (float ****)safe_malloc(w * sizeof(float ***));
	for (h = 0; h < w; ++ h) {
		mat4D_float[h] = (float ***)safe_malloc(x * sizeof(float **));
		for (i = 0; i < x; ++ i) {
			mat4D_float[h][i] = (float **)safe_malloc(y * sizeof(float *));
			for (j = 0; j < y; ++ j)
				mat4D_float[h][i][j] = (float *)safe_malloc(z * sizeof(float));
		}
	}

	return mat4D_float;
}

/** initialise */
void init_mat4D_float(float ****mat4D_float, int w, int x, int y, int z, float val)
{
    unsigned int h, i, j, k;

	for (h = 0; h < w; ++ h)
		for (i = 0; i < x; ++ i)
			for (j = 0; j < y; ++ j)
				for (k = 0; k < z; ++ k)
					mat4D_float[h][i][j][k] = val;
}


/** free */
void free_mat4D_float(float ****mat4D_float, int w, int x, int y)
{
    unsigned int h, i, j;

    for (h = 0; h < w; ++ h) {
		for (i = 0; i < x; ++ i) {
			for (j = 0; j < y; ++ j)
				free(mat4D_float[h][i][j]);
			free(mat4D_float[h][i]);
	    }
		free(mat4D_float[h]);
	}
	free(mat4D_float);
}

/*___________________________________________________________________________*/
/** DOUBLE */
/*___________________________________________________________________________*/
/** 2D double matrix */
/** allocate */
double **alloc_mat2D_double(double **mat2D_double, int x, int y)
{
    unsigned int i;

	mat2D_double = (double **)safe_malloc(x * sizeof(double *));
	for (i = 0; i < x; ++ i)
		mat2D_double[i] = (double *)safe_malloc(y * sizeof(double));

	return mat2D_double;
}

/** initialise */
void init_mat2D_double(double **mat2D_double, int x, int y, double val)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			mat2D_double[i][j] = val;
}

/* free */
void free_mat2D_double(double **mat2D_double, int x)
{
    unsigned int i;

	for (i = 0; i < x; ++ i)
		free(mat2D_double[i]);
    free(mat2D_double);
}

/** print */
void print_mat2D_double(char *outFileName, double **mat2D_double, int x, int y)
{
    unsigned int i, j;
    FILE *outFile = 0;

    outFile = safe_open(outFileName, "w");

	for (i = 0; i < x; ++ i) {
		if (i > 0) fprintf(outFile, "\n");
		for (j = 0; j < y; ++ j) {
			fprintf(outFile, "%4.2f ", mat2D_double[i][j]);
		}
	}
	fprintf(outFile, "\n");

	fclose(outFile);
}

/** print with lower limit specification */
void print_mat2D_double_lowlim(char *outFileName, double **mat2D_double, int lowx, int x, int lowy, int y)
{
    unsigned int i, j;
    FILE *outFile = 0;

    outFile = safe_open(outFileName, "w");

	for (i = lowx; i < x; ++ i) {
		if (i > lowx) fprintf(outFile, "\n");
		for (j = lowy; j < y; ++ j) {
			fprintf(outFile, "%4.2f ", mat2D_double[i][j]);
		}
	}
	fprintf(outFile, "\n");

	fclose(outFile);
}

/** print with lower limit specification */
void print_mat2D_doublee_lowlim(char *outFileName, double **mat2D_double, int lowx, int x, int lowy, int y)
{
    unsigned int i, j;
    FILE *outFile = 0;

    outFile = safe_open(outFileName, "w");

	for (i = lowx; i < x; ++ i) {
		if (i > lowx) fprintf(outFile, "\n");
		for (j = lowy; j < y; ++ j) {
			fprintf(outFile, "%10.5e ", mat2D_double[i][j]);
		}
	}
	fprintf(outFile, "\n");

	fclose(outFile);
}

/** divide by factor */
void div_mat2D_double(double **mat2D_double, int x, int y, double a)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i) {
		for (j = 0; j < y; ++ j) {
			mat2D_double[i][j] /= a;
		}
	}
}

/** symmetry check */
int symmetry_mat2D_double(double **mat2D_double, int x, int y)
{
    unsigned int i, j, s;

	for (i = 0, s = 0; i < x; ++ i)
		for (j = 0; j < i; ++ j)
			s += (mat2D_double[i][j] == mat2D_double[j][i] ? 0 : 1);
			if (s)
				fprintf(stderr, "asymmetry at [%d %d] = %lf, [%d %d] = %lf\n", 
					i, j, mat2D_double[i][j], j, i, mat2D_double[j][i]);

	return s;
}

/* copy matrix */
void copy_mat2D_double(double **mat2D_double_to, double **mat2D_double_from, int x, int y)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			mat2D_double_to[i][j] = mat2D_double_from[i][j];
}

/** add matrix */
void add_mat2D_double(double **mat2D_double_sum, double **mat2D_double_A, double **mat2D_double_B, int x, int y)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			mat2D_double_sum[i][j] = mat2D_double_A[i][j] + mat2D_double_B[i][j];
}

/*___________________________________________________________________________*/
/** 3D double matrix */
/** allocate */
double ***alloc_mat3D_double(double ***mat3D_double, int x, int y, int z)
{
	unsigned int i, j;

	mat3D_double = (double ***)safe_malloc(x * sizeof(double **));
	for (i = 0; i < x; ++ i) {
		mat3D_double[i] = (double **)safe_malloc(y * sizeof(double *));
		for (j = 0; j < y; ++ j)
			mat3D_double[i][j] = (double *)safe_malloc(z * sizeof(double));
	}

	return mat3D_double;
}

/** initialise */
void init_mat3D_double(double ***mat3D_double, int x, int y, int z, double val)
{
    unsigned int i, j, k;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			for (k = 0; k < z; ++ k)
				mat3D_double[i][j][z] = val;
}


/** free */
void free_mat3D_double(double ***mat3D_double, int x, int y)
{
    unsigned int i, j;

    for (i = 0; i < x; ++ i) {
        for (j = 0; j < y; ++ j)
            free(mat3D_double[i][j]);
        free(mat3D_double[i]);
    }
    free(mat3D_double);
}

/*___________________________________________________________________________*/
/** 4D double matrix */
/** allocate */
double ****alloc_mat4D_double(double ****mat4D_double, int w, int x, int y, int z)
{
	unsigned int h, i, j;

	mat4D_double = (double ****)safe_malloc(w * sizeof(double ***));
	for (h = 0; h < w; ++ h) {
		mat4D_double[h] = (double ***)safe_malloc(x * sizeof(double **));
		for (i = 0; i < x; ++ i) {
			mat4D_double[h][i] = (double **)safe_malloc(y * sizeof(double *));
			for (j = 0; j < y; ++ j)
				mat4D_double[h][i][j] = (double *)safe_malloc(z * sizeof(double));
		}
	}

	return mat4D_double;
}

/** initialise */
void init_mat4D_double(double ****mat4D_double, int w, int x, int y, int z, double val)
{
    unsigned int h, i, j, k;

	for (h = 0; h < w; ++ h)
		for (i = 0; i < x; ++ i)
			for (j = 0; j < y; ++ j)
				for (k = 0; k < z; ++ k)
					mat4D_double[h][i][j][k] = val;
}


/** free */
void free_mat4D_double(double ****mat4D_double, int w, int x, int y)
{
    unsigned int h, i, j;

    for (h = 0; h < w; ++ h) {
		for (i = 0; i < x; ++ i) {
			for (j = 0; j < y; ++ j)
				free(mat4D_double[h][i][j]);
			free(mat4D_double[h][i]);
	    }
		free(mat4D_double[h]);
	}
	free(mat4D_double);
}

/*___________________________________________________________________________*/
/** VECTOR */
/*___________________________________________________________________________*/
/** 2D vector matrix */
/** allocate */
Vec **alloc_mat2D_vec(Vec **mat2D_vec, int x, int y)
{
    unsigned int i;

	mat2D_vec = (Vec **)safe_malloc(x * sizeof(Vec *));
	for (i = 0; i < x; ++ i)
		mat2D_vec[i] = (Vec *)safe_malloc(y * sizeof(Vec));

	return mat2D_vec;
}

/** initialise */
void init_mat2D_vec(Vec **mat2D_vec, int x, int y, Vec val)
{
    unsigned int i, j;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			v_copy(&(mat2D_vec[i][j]), &(val));
}

/** free */
void free_mat2D_vec(Vec **mat2D_vec, int x)
{
    unsigned int i;

	for (i = 0; i < x; ++ i)
		free(mat2D_vec[i]);
    free(mat2D_vec);
}

/*___________________________________________________________________________*/
/** 3D vector matrix */
/** allocate */
Vec ***alloc_mat3D_vec(Vec ***mat3D_vec, int x, int y, int z)
{
	unsigned int i, j;

	mat3D_vec = (Vec ***)safe_malloc(x * sizeof(Vec **));
	for (i = 0; i < x; ++ i) {
		mat3D_vec[i] = (Vec **)safe_malloc(y * sizeof(Vec *));
		for (j = 0; j < y; ++ j)
			mat3D_vec[i][j] = (Vec *)safe_malloc(z * sizeof(Vec));
	}

	return mat3D_vec;
}

/** initialise */
void init_mat3D_vec(Vec ***mat3D_vec, int x, int y, int z, Vec val)
{
    unsigned int i, j, k;

	for (i = 0; i < x; ++ i)
		for (j = 0; j < y; ++ j)
			for (k = 0; k < z; ++ k)
				v_copy(&(mat3D_vec[i][j][z]), &(val));
}

/** free */
void free_mat3D_vec(Vec ***mat3D_vec, int x, int y)
{
    unsigned int i, j;

    for (i = 0; i < x; ++ i) {
        for (j = 0; j < y; ++ j)
            free(mat3D_vec[i][j]);
        free(mat3D_vec[i]);
    }
    free(mat3D_vec);
}

/*___________________________________________________________________________*/
/** 4D vector matrix */
/** allocate */
Vec ****alloc_mat4D_vec(Vec ****mat4D_vec, int w, int x, int y, int z)
{
	unsigned int h, i, j;

	mat4D_vec = (Vec ****)safe_malloc(w * sizeof(Vec ***));
	for (h = 0; h < w; ++ h) {
		mat4D_vec[h] = (Vec ***)safe_malloc(x * sizeof(Vec **));
		for (i = 0; i < x; ++ i) {
			mat4D_vec[h][i] = (Vec **)safe_malloc(y * sizeof(Vec *));
			for (j = 0; j < y; ++ j)
				mat4D_vec[h][i][j] = (Vec *)safe_malloc(z * sizeof(Vec));
		}
	}

	return mat4D_vec;
}

/** initialise */
void init_mat4D_vec(Vec ****mat4D_vec, int w, int x, int y, int z, Vec val)
{
    unsigned int h, i, j, k;

	for (h = 0; h < w; ++ h)
		for (i = 0; i < x; ++ i)
			for (j = 0; j < y; ++ j)
				for (k = 0; k < z; ++ k)
					v_copy(&(mat4D_vec[h][i][j][k]), &(val));
}

/** free */
void free_mat4D_vec(Vec ****mat4D_vec, int w, int x, int y)
{
    unsigned int h, i, j;

    for (h = 0; h < w; ++ h) {
		for (i = 0; i < x; ++ i) {
			for (j = 0; j < y; ++ j)
				free(mat4D_vec[h][i][j]);
			free(mat4D_vec[h][i]);
	    }
		free(mat4D_vec[h]);
	}
	free(mat4D_vec);
}

