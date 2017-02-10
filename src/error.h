/*==============================================================================
error.h : error message routines
(C) 2008 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
==============================================================================*/

#ifndef ERROR_H
#define ERROR_H

#include <stdio.h>
#include <stdlib.h>

/*____________________________________________________________________________*/
/* protypes */
void Warning(char *message);
void WarningSpec(char *message, char *spec);
void Error(char *message);
void ErrorSpec(char *message, char *spec); 
void ErrorSpecNoexit(char *message, char *spec); 
void ErrorLoc(char *message, char *file, int line);

#endif
