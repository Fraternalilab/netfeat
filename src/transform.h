/*==============================================================================
transform.h : transform network matricess
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "config.h"
#include "error.h"
#include "float.h"
#include "ints.h"
#include "matrix.h"
#include "putprop.h"

/*____________________________________________________________________________*/
/* prototypes */
void smoothen_Pi(Ints *ints, Arg *arg);
void smoothen_p(Ints *ints, Arg *arg);
void smoothen_W(Ints *ints, Arg *arg);
void recompute_w(Ints *ints, Dist *dist, Arg *arg);
void norm_corr_mat(double **matNorm, double **matIn, int x);

#endif
