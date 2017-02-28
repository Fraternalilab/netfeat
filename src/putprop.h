/*==============================================================================
putprop.h : print properties
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#ifndef PUTPROP_H
#define PUTPROP_H

#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "array.h"
#include "config.h"
#include "ints.h"
#include "safe.h"
#include "matrix.h"

/*____________________________________________________________________________*/
/* prototypes */
void print_k_list(int n, Ints *ints);
void print_k(int n, Ints *ints);
void print_p(int n, Ints *ints);
void print_P_conn(int n, Ints *ints);
void print_Pi_corr(int n, int level, Ints *ints);
void print_etc(int n, Ints *ints);
void print_entropy(int n, Ints *ints);

#endif
