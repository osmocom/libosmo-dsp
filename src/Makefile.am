# This is _NOT_ the library release version, it's an API version.
# Please read Chapter 6 "Library interface versions" of the libtool documentation before making any modification
LIBVERSION=0:0:0

INCLUDES = $(all_includes) -I$(top_srcdir)/include
AM_CFLAGS = -fPIC -Wall -ffast-math

lib_LTLIBRARIES = libosmodsp.la

libosmodsp_la_SOURCES = cfile.c cxvec.c cxvec_math.c
libosmodsp_la_LDFALGS = -version-info $(LIBVERSION)
