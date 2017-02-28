#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(netFeat, 1.4.2, jens@jkleinj.eu)

AC_PREREQ(2.61)
AC_CONFIG_SRCDIR([src/safe.c])
AC_CONFIG_HEADER([src/config.h])

AM_INIT_AUTOMAKE

AC_SUBST(INTI_CFLAGS)
AC_SUBST(INTI_LIBS)

# Checks for programs.
AC_PROG_CC
AC_PROG_LN_S

# Checks for libraries.
AC_CHECK_LIB([m], [cos])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# enable debug
AC_ARG_ENABLE(debug, [--enable-debug  Turn on debugging], [ac_use_debug="yes"], [ac_use_debug="no"])
if test "$ac_use_debug" = "yes"; then
    CFLAGS="-ggdb"
    AC_DEFINE(DEBUG,1,[Define to compile with DEBUG settings.])
fi

# enable profiling 
AC_ARG_ENABLE(profiling, [--enable-profiling  Turn on profiling], [ac_use_profiling="yes"], [ac_use_profiling="no"])
if test "$ac_use_profiling" = "yes"; then
    CFLAGS="-pg -O2"
    AC_DEFINE(PROFILING,1,[Define to compile with PROFILING settings.])
fi

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile]
				[src/Makefile]
				[tests/Makefile])
AC_OUTPUT

