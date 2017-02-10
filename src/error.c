/*==============================================================================
error.c : error message routines
(C) 2008 Jens Kleinjung and Alessandro Pandini
Read the COPYING file for license information.
==============================================================================*/

#include "error.h"

/*____________________________________________________________________________*/
/** print warning message */
void Warning(char *message) {
    fprintf(stderr, "Warning: %s\n", message);
}

/*____________________________________________________________________________*/
/** print warning message with specification */
void WarningSpec(char *message, char *spec) {
    fprintf(stderr, "Warning: %s (%s)\n", message, spec);
}

/*____________________________________________________________________________*/
/** print error message and exit */
void Error(char *message) {
    fprintf(stderr, "Error: %s\n", message);
    exit(1);
}

/*____________________________________________________________________________*/
/** print error message with specification and exit */
void ErrorSpec(char *message, char *spec) {
    fprintf(stderr, "Error: %s (%s)\n", message, spec);
    exit(1);
}

/*____________________________________________________________________________*/
/** print error message with specification, don't exit */
void ErrorSpecNoexit(char *message, char *spec) {
    fprintf(stderr, "Error: %s (%s)\n", message, spec);
}

/*____________________________________________________________________________*/
/** print error message and exit:
	call with arguments __FILE__ and __LINE__ */
void ErrorLoc(char *message, char *file, int line) {
    fprintf(stderr, "Error: %s (%s:%d)\n", message, file, line);
    exit(1);
}
