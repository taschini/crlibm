## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libcrlibm.a

# shared libraries need libtool. We'll do it later.
#lib_LTLIBRARIES = libcrlibm.la

libcrlibm_a_SOURCES =  \
	crlibm.h crlibm_private.h crlibm_private.c\
	exp.h exp.c exp_fast.h exp_fast.c\
	log.h log.c log_fast.c log_fast.h\
	log10.h	log10.c\
	log2.h  log2.c \
	coefpi2.h rem_pio2.c\
	sine.h sine.c\
	cosine.h cosine.c\
	tan.h tan.c cotan.c\
	atan.h atan.c 

# add all the scs_lib useful functions
libcrlibm_a_LIBADD = scs_lib/addition_scs.o \
	scs_lib/division_scs.o scs_lib/print_scs.o\
	scs_lib/double2scs.o scs_lib/zero_scs.o\
	scs_lib/multiplication_scs.o scs_lib/scs2double.o

include_HEADERS = crlibm.h

SUBDIRS = scs_lib . tests

AM_CFLAGS =  -fPIC -Wall -Wshadow -Wpointer-arith  -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline 
EXTRA_DIST = VERSION  
