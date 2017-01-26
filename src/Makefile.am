#===============================================================================
# Makefile.am : for automake
# Copyright (C) 2016 Jens Kleinjung
#===============================================================================

AM_CPPFLAGS = $(INTI_CFLAGS)
AM_CFLAGS = -Wall

#_______________________________________________________________________________
bin_PROGRAMS = ga readpoolbin

#_______________________________________________________________________________
ga_SOURCES = \
freqs.c freqs_data.h freqs.h ga.c ga.h getfreqs.c getfreqs.h \
getseqs.c getseqs.h huffman.c huffman.h lz.c lz.h \
safe.c safe.h suffix_tree.c suffix_tree.h

ga_LDADD = $(INTI_LIBS)

#_______________________________________________________________________________
readpoolbin_SOURCES = \
read_poolbin.c

readpoolbin_LDADD = $(INTI_LIBS)

#_______________________________________________________________________________
EXTRA_DIST = doxygen.cfg

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

