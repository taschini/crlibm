## Process this file with automake to produce Makefile.in


lib_LIBRARIES = libscs.a


# shared libraries need libtool. We'll do it later.
#lib_LTLIBRARIES = libscs.la

libscs_a_SOURCES = scs.h scs_private.h scs_private.c \
	double2scs.c print_scs.c \
	addition_scs.c multiplication_scs.c division_scs.c\
	rand_scs.c poly_fct.c\
	scs2double.c zero_scs.c \
	scs2mpf.c scs2mpfr.c \
	wrapper_scs.h


include_HEADERS = scs.h scs_config.h

SUBDIRS = . tests 

AM_CFLAGS =  -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline 

EXTRA_DIST = VERSION  Doxyfile.user Doxyfile.dev 
