/*==============================================================================
getinteraction.h : interaction input as list or matrix
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#ifndef GETINTERACTION_H
#define GETINTERACTION_H

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "config.h"
#include "ints.h"
#include "matrix.h"
#include "safe.h"
#include "vector.h"

void read_node_list(char *protFileName, Ints *ints, Arg *arg);
void read_interaction_list(char *intsFileName, Ints *ints, Arg *arg);
void assign_interaction_matrix(int n, Ints *ints, Arg *arg);
int read_interaction_matrix(char *matFileName, Ints *ints, Arg *arg, int alldata);

#endif
