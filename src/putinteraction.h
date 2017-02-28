/*==============================================================================
putinteraction.h : interaction output as list or matrix
(C) 2008-2010 A Annibale, L Fernandes, ACC Coolen, J Kleinjung, F Fraternali
Read the COPYING file for license information.
==============================================================================*/

#ifndef PUTINTERACTION_H
#define PUTINTERACTION_H

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>

#include "arg.h"
#include "config.h"
#include "ints.h"
#include "safe.h"

void print_interaction_matrix_ascii(char *outFileName, Ints *ints, Arg *arg);
void print_interaction_matrix_binary(char *outFileName, Ints *ints);
void print_interaction_list_ascii(char *outFileName, Ints *ints, Arg *arg);

#endif
