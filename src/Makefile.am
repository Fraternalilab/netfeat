#===============================================================================
# Makefile.am : automake input
#===============================================================================

bin_PROGRAMS = netfeat

AM_CPPFLAGS = $(INTI_CFLAGS)

netfeat_SOURCES	= \
	arg.c arg.h array.c array.h config.h error.c error.h \
	getinteraction.c getinteraction.h ints.h loops.c loops.h matrix.c matrix.h \
	netcons.c netcons.h netdist.c netdist.h netfeat.c netfeat.h netprop.c netprop.h \
	putprop.c putprop.h putinteraction.c putinteraction.h netmult.c netmult.h \
	safe.c safe.h transform.c transform.h vector.c vector.h

netfeat_LDADD = $(INTI_LIBS)

EXTRA_DIST = doxygen.cfg

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

