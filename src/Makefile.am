# This is _NOT_ the library release version, it's an API version.
# Please read Chapter 6 "Library interface versions" of the libtool documentation before making any modification
LIBVERSION=0:1:0

AM_CPPFLAGS = $(all_includes) -I$(top_srcdir)/include
AM_CFLAGS = -Wall -ffast-math $(FFTW3F_CFLAGS)

lib_LTLIBRARIES = libosmodsp.la

libosmodsp_la_SOURCES = cfile.c cxvec.c cxvec_math.c iqbal.c
libosmodsp_la_LDFLAGS = -version-info $(LIBVERSION) $(FFTW3F_LIBS)
